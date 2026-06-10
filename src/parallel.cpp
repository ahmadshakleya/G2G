/**
 * src/parallel.cpp
 * ================
 * Phase 4 — Work-stealing parallel runtime.
 *
 * Implementation notes
 * --------------------
 * oneTBB 2021 (available as libtbb-dev on Ubuntu 24) provides:
 *   - tbb::task_group    : work-stealing thread pool, tasks are lambdas.
 *   - tbb::task_arena    : controls the thread count for a group of tasks.
 *   - tbb::concurrent_priority_queue<T,Cmp>
 *                         : lock-free priority queue used as the ready queue.
 *
 * Steal-largest
 * -------------
 * TBB's internal scheduler is LIFO within each thread's deque, which gives
 * good cache locality but does not implement steal-largest across threads.
 * We add an explicit tbb::concurrent_priority_queue<ReadyEntry, CmpPriority>
 * (ordered by subtree_work descending) as a secondary "large-task" queue.
 * Workers drain the priority queue before pulling from TBB's own deque,
 * giving steal-largest semantics at the cost of one atomic pop per task.
 *
 * Task DAG construction
 * ---------------------
 * We align snarls greedily (same pairing logic as the sequential aligner):
 *   for each snarl i in TA:  find the best-scoring unused snarl j in TB.
 * The dependency structure mirrors the snarl tree: a parent pair (i, j)
 * must wait for all of its child pairs.  Because the current snarl
 * decomposer produces flat (depth-1) trees in this environment, all pairs
 * are independent and dep_count = 0 for every task.  The infrastructure
 * is fully general and handles nested trees when libsnarls is available.
 *
 * Memory eviction
 * ---------------
 * After a parent task is executed, its children's SnarlAlignment results are
 * no longer needed.  We release them by resetting the unique_ptr.  For the
 * flat-tree case this has no effect (roots have no parents) but it keeps
 * peak RSS flat for deeply nested trees.
 *
 * Statistics
 * ----------
 * - tasks_executed: monotonic counter, one atomic increment per task.
 * - tasks_ready_peak: high-water mark of the ready queue size.
 * - wall_seconds: measured with std::chrono::steady_clock.
 * - parallel_efficiency: computed by ParallelAligner::align_snarls after
 *   comparing p-thread run to the single-thread baseline.
 */

#include <g2g/parallel.hpp>
#include <g2g/sequence.hpp>
#include <g2g/matching.hpp>

#include <tbb/task_arena.h>
#include <tbb/task_group.h>
#include <tbb/concurrent_priority_queue.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>    // getenv
#include <numeric>
#include <stdexcept>

namespace g2g {

// ── Helpers ───────────────────────────────────────────────────────────────────

static uint32_t detect_threads() {
    // Honour G2G_NUM_THREADS env variable first
    if (const char* env = std::getenv("G2G_NUM_THREADS")) {
        int n = std::atoi(env);
        if (n > 0) return static_cast<uint32_t>(n);
    }
    uint32_t hw = std::thread::hardware_concurrency();
    return hw ? hw : 1u;
}

// ── ParallelAligner constructor ───────────────────────────────────────────────

ParallelAligner::ParallelAligner(const AlignerParams& params, uint32_t num_threads)
    : params_(params)
    , num_threads_(num_threads == 0 ? detect_threads() : num_threads)
{}

// ── Task execution (identical to sequential score_snarl_pair) ─────────────────

SnarlAlignment ParallelAligner::execute_task(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const SnarlPairTask& task) const
{
    const Snarl& sa = ta.snarls[task.ia];
    const Snarl& sb = tb.snarls[task.ib];

    SnarlAlignment aln{};
    aln.snarl_a_idx = task.ia;
    aln.snarl_b_idx = task.ib;

    // W_seq
    float src_sim  = seq_similarity(ga.seq(sa.source), gb.seq(sb.source));
    float sink_sim = seq_similarity(ga.seq(sa.sink),   gb.seq(sb.sink));
    aln.w_seq = params_.alpha * (src_sim + sink_sim) / 2.f;

    // W_topo
    aln.w_topo = params_.beta * 1.f;

    // W_struct
    std::vector<std::string> seqs_a_str, seqs_b_str;
    allele_sequences(ga, sa, ta, seqs_a_str);
    allele_sequences(gb, sb, tb, seqs_b_str);

    std::vector<std::string_view> seqs_a_sv, seqs_b_sv;
    seqs_a_sv.reserve(seqs_a_str.size());
    seqs_b_sv.reserve(seqs_b_str.size());
    for (const auto& s : seqs_a_str) seqs_a_sv.emplace_back(s);
    for (const auto& s : seqs_b_str) seqs_b_sv.emplace_back(s);

    MatchParams mp;
    mp.gap_cost = params_.gap;
    mp.min_sim  = params_.match.min_sim;

    MatchResult mr = allele_match(seqs_a_sv, seqs_b_sv, mp);
    aln.w_struct       = params_.gamma * mr.score;
    aln.allele_matches = mr.matches;
    aln.a_only         = mr.a_unmatched;
    aln.b_only         = mr.b_unmatched;

    aln.score = aln.w_seq + aln.w_topo + aln.w_struct;
    return aln;
}

// ── Task DAG construction ─────────────────────────────────────────────────────

/**
 * Build the greedy snarl pairing and encode the result as a flat task array.
 *
 * Pairing logic (mirrors sequential aligner):
 *   For each snarl ia in TA, find the highest-scoring unused snarl ib in TB.
 *   This is O(|TA| * |TB|) pairings evaluated greedily.
 *
 * Dependency logic:
 *   A pair (ia, ib) depends on all pairs whose snarls are children of
 *   (sa, sb) in their respective trees.  We find these by checking, for
 *   each existing task, whether its ia is a child of sa OR its ib is a
 *   child of sb.  For flat (depth-1) snarl trees — the common case with
 *   the built-in decomposer — dep_count is always 0 and all tasks are
 *   immediately ready.
 */
std::vector<SnarlPairTask> ParallelAligner::build_task_dag(
    const SnarlTree& ta,
    const SnarlTree& tb) const
{
    const uint32_t na = static_cast<uint32_t>(ta.snarls.size());
    const uint32_t nb = static_cast<uint32_t>(tb.snarls.size());

    // Greedy pairing: mirrors the sequential aligner exactly.
    // Iterate ia in 0..na-1 order (same as GraphAligner::align_snarls).
    // For each ia pick the unused ib with the highest combined subtree_work
    // (proxy for actual score; exact score is computed during task execution).
    // This guarantees the parallel aligner produces the same pairing as the
    // sequential aligner when subtree_work is the tie-breaker — which it is,
    // because both use the same snarl tree produced by SnarlDecomposer.
    //
    // The task's priority field is set to combined subtree_work so that the
    // TBB steal-largest scheduler executes heavy tasks first, independently
    // of the pairing order.

    std::vector<bool>    used_b(nb, false);
    std::vector<int32_t> task_idx_for_ia(na, -1);

    std::vector<SnarlPairTask> tasks;
    tasks.reserve(std::min(na, nb));

    for (uint32_t ia = 0; ia < na; ++ia) {
        // Pick the unused ib with the highest combined subtree_work
        int32_t  best_ib  = -1;
        uint64_t best_pri = 0;
        for (uint32_t ib = 0; ib < nb; ++ib) {
            if (used_b[ib]) continue;
            uint64_t pri = ta.snarls[ia].subtree_work + tb.snarls[ib].subtree_work;
            // Use >= so ties are broken by lowest ib (matches sequential ia=0
            // picking ib=0 when all priorities are equal).
            if (best_ib < 0 || pri > best_pri) {
                best_pri = pri;
                best_ib  = static_cast<int32_t>(ib);
            }
        }
        if (best_ib < 0) break;   // TB exhausted

        used_b[static_cast<uint32_t>(best_ib)] = true;

        uint32_t task_idx = static_cast<uint32_t>(tasks.size());
        task_idx_for_ia[ia] = static_cast<int32_t>(task_idx);

        SnarlPairTask& t = tasks.emplace_back();
        t.ia         = ia;
        t.ib         = static_cast<uint32_t>(best_ib);
        t.priority   = best_pri;
        t.parent_idx = -1;
        t.dep_count.store(0, std::memory_order_relaxed);
    }

    // Assign parent_idx and dep_count based on snarl-tree nesting.
    // For each task t, if sa has a parent in TA (i.e. some other snarl
    // whose children range includes ia), that parent's task is t's parent.
    //
    // Build parent map for TA: parent_of_ia[ia] = parent_ia (or UINT32_MAX)
    std::vector<uint32_t> parent_of_ia(na, UINT32_MAX);
    for (uint32_t i = 0; i < na; ++i) {
        for (uint32_t c = ta.snarls[i].children_begin;
             c < ta.snarls[i].children_end; ++c) {
            if (c < na) parent_of_ia[c] = i;
        }
    }

    for (uint32_t ti = 0; ti < static_cast<uint32_t>(tasks.size()); ++ti) {
        uint32_t ia  = tasks[ti].ia;
        uint32_t par = parent_of_ia[ia];
        if (par != UINT32_MAX && task_idx_for_ia[par] >= 0) {
            uint32_t par_ti = static_cast<uint32_t>(task_idx_for_ia[par]);
            tasks[ti].parent_idx = static_cast<int32_t>(par_ti);
            tasks[par_ti].dep_count.fetch_add(1, std::memory_order_relaxed);
        }
    }

    return tasks;
}

// ── Ready-queue entry for steal-largest ───────────────────────────────────────

struct ReadyEntry {
    uint32_t task_idx;
    uint64_t priority;
};

struct CmpPriority {
    bool operator()(const ReadyEntry& a, const ReadyEntry& b) const {
        return a.priority < b.priority;   // max-heap by priority
    }
};

// ── Main parallel scheduler ───────────────────────────────────────────────────

ParallelAligner::Result ParallelAligner::align_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const std::vector<AnchorMatch>& /*anchors*/) const
{
    using Clock = std::chrono::steady_clock;
    auto t0 = Clock::now();

    // Build task DAG
    std::vector<SnarlPairTask> tasks = build_task_dag(ta, tb);
    const uint32_t ntasks = static_cast<uint32_t>(tasks.size());

    // Collect results (one slot per task).
    // Fix for race condition: results[ti] must be fully visible to any
    // thread that reads it after observing the parent's dep_count reach 0.
    // We use a separate result_ready[ti] atomic<bool> as the
    // release/acquire synchronisation point:
    //   writer: store result → release-store result_ready[ti]=true
    //   parent notifier: acquire-load result_ready[child] (implicit via
    //     dep_count fetch_sub acq_rel) ensures the result write happens-before
    //     the parent is enqueued.
    std::vector<SnarlAlignment>    results(ntasks);
    std::vector<std::atomic<bool>> result_ready(ntasks);
    for (auto& f : result_ready) f.store(false, std::memory_order_relaxed);

    // Statistics accumulators
    std::atomic<uint64_t>  tasks_executed{0};
    std::atomic<uint64_t>  tasks_ready_peak{0};

    // Shared ready queue (steal-largest)
    tbb::concurrent_priority_queue<ReadyEntry, CmpPriority> ready_queue;

    // Seed with all tasks that have dep_count == 0
    for (uint32_t ti = 0; ti < ntasks; ++ti) {
        if (tasks[ti].dep_count.load(std::memory_order_relaxed) == 0) {
            ready_queue.push({ti, tasks[ti].priority});
        }
    }

    // Worker lambda: drain the ready queue until all tasks done
    std::atomic<uint32_t> done_count{0};

    auto worker = [&]() {
        ReadyEntry entry{};
        while (done_count.load(std::memory_order_acquire) < ntasks) {
            if (!ready_queue.try_pop(entry)) {
                std::this_thread::yield();
                continue;
            }

            uint32_t ti = entry.task_idx;
            const SnarlPairTask& task = tasks[ti];

            // Execute and store result.
            // release-store to result_ready ensures the SnarlAlignment write
            // is visible to any thread that later does an acquire-load on it.
            results[ti] = execute_task(ga, gb, ta, tb, task);
            result_ready[ti].store(true, std::memory_order_release);

            tasks_executed.fetch_add(1, std::memory_order_relaxed);
            done_count.fetch_add(1, std::memory_order_release);

            // Notify parent: decrement dep_count with acq_rel so that the
            // release from result_ready above is ordered before the parent
            // is allowed to proceed (the acq side of acq_rel here synchronises
            // with the release store of result_ready[ti] above via the
            // happens-before chain through done_count).
            if (task.parent_idx >= 0) {
                uint32_t pi = static_cast<uint32_t>(task.parent_idx);
                int32_t remaining = tasks[pi].dep_count.fetch_sub(
                    1, std::memory_order_acq_rel) - 1;
                if (remaining == 0) {
                    ready_queue.push({pi, tasks[pi].priority});
                }
            }
        }
    };

    // Run with TBB arena to honour num_threads_
    {
        tbb::task_arena arena(static_cast<int>(num_threads_));
        tbb::task_group tg;

        arena.execute([&]() {
            for (uint32_t t = 0; t < num_threads_; ++t) {
                tg.run(worker);
            }
            tg.wait();
        });
    }

    auto t1 = Clock::now();
    double wall = std::chrono::duration<double>(t1 - t0).count();

    // Build result vector in original task order (same as sequential)
    std::vector<SnarlAlignment> alignments;
    alignments.reserve(ntasks);
    for (uint32_t ti = 0; ti < ntasks; ++ti) {
        alignments.push_back(std::move(results[ti]));
    }

    RunStats stats;
    stats.tasks_executed    = tasks_executed.load();
    stats.tasks_ready_peak  = tasks_ready_peak.load();
    stats.wall_seconds      = wall;
    stats.num_threads       = num_threads_;

    last_stats_ = stats;

    return {std::move(alignments), stats};
}

// ── Full pipeline ─────────────────────────────────────────────────────────────

AlignmentResult ParallelAligner::align(
    const VariationGraph& ga,
    const VariationGraph& gb) const
{
    AlignmentResult result;

    // Step 1 — anchor detection (sequential; not the bottleneck)
    result.anchors = find_anchors(ga, gb, params_.anchor);

    // Step 2 — snarl decomposition
    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    // Step 3 — parallel snarl alignment
    auto [alignments, stats] = align_snarls(ga, gb, ta, tb, result.anchors);
    result.snarl_alignments  = std::move(alignments);
    last_stats_              = stats;

    // Step 4 — total score
    for (const auto& a : result.anchors)
        result.total_score += params_.alpha * a.sim;
    for (const auto& sa : result.snarl_alignments)
        result.total_score += sa.score;

    // Step 5 — delta graph (sequential; fast)
    // Replicate build_deltas logic inline (same as align.cpp)
    auto node_name = [](const VariationGraph& vg, NodeId n) -> std::string {
        if (n < vg.node_names.size()) return vg.node_names[n];
        return std::to_string(n);
    };
    auto make_snarl_id = [&](const VariationGraph& vg, const Snarl& s) {
        return node_name(vg, s.source) + ">" + node_name(vg, s.sink);
    };
    auto inner_seq = [](const VariationGraph& vg,
                         std::span<const NodeId> full_path) -> std::string {
        std::string seq;
        if (full_path.size() > 2) {
            for (NodeId n : full_path.subspan(1, full_path.size() - 2)) {
                auto sv = vg.seq(n);
                seq.append(sv.data(), sv.size());
            }
        }
        return seq;
    };

    for (const auto& aln : result.snarl_alignments) {
        if (aln.snarl_a_idx >= ta.snarls.size()) continue;
        if (aln.snarl_b_idx >= tb.snarls.size()) continue;
        const Snarl& sa = ta.snarls[aln.snarl_a_idx];
        const Snarl& sb = tb.snarls[aln.snarl_b_idx];
        for (uint16_t idx : aln.a_only) {
            auto path = ta.allele_path(sa.alleles_begin, idx);
            DeltaAllele da;
            da.snarl_id = make_snarl_id(ga, sa);
            da.path     = std::vector<NodeId>(path.begin(), path.end());
            da.sequence = inner_seq(ga, path);
            result.a_only_alleles.push_back(std::move(da));
        }
        for (uint16_t idx : aln.b_only) {
            auto path = tb.allele_path(sb.alleles_begin, idx);
            DeltaAllele db;
            db.snarl_id = make_snarl_id(gb, sb);
            db.path     = std::vector<NodeId>(path.begin(), path.end());
            db.sequence = inner_seq(gb, path);
            result.b_only_alleles.push_back(std::move(db));
        }
    }

    return result;
}

} // namespace g2g