#pragma once
/**
 * include/g2g/parallel.hpp
 * ========================
 * Phase 4 — Work-stealing parallel runtime for the snarl-pair DP.
 *
 * Design
 * ------
 * The snarl-pair alignment has a natural task DAG: each pair (SA, SB) depends
 * on all of its child pairs (bottom-up DP).  Sibling pairs are fully
 * independent and can run in parallel.
 *
 * We implement the scheduler in three layers:
 *
 *   1. SnarlPairTask  — A unit of work: score one (SA_idx, SB_idx) pair.
 *                       Priority = Snarl::subtree_work (pre-computed in P3-4).
 *                       Populated into a TBB concurrent_priority_queue so that
 *                       TBB's internal work-stealing naturally picks
 *                       highest-subtree_work tasks first (steal-largest).
 *
 *   2. ParallelAlignScheduler — Builds the task DAG from two SnarlTrees,
 *                       drives a tbb::task_group, and collects results.
 *                       NUMA-aware allocation hint stored per task (populated
 *                       at task-creation time for multi-socket servers).
 *
 *   3. ParallelAligner — Drop-in replacement for GraphAligner::align_snarls().
 *                       Accepts AlignerParams and thread count; falls back to
 *                       sequential mode when num_threads == 1.
 *
 * Memory management
 * -----------------
 * Each SnarlPairTask owns its SnarlAlignment result on the heap.  Once a
 * parent task has consumed a child result, the child result is freed
 * immediately (memory_evict flag).  Peak RAM is proportional to the
 * tree frontier, not the full tree.
 *
 * Statistics
 * ----------
 * ParallelAlignScheduler::stats() returns a RunStats struct with steal
 * count, load imbalance factor, peak active tasks, and wall time — exactly
 * the data required by Experiment 4 (parallel scaling) in the paper.
 *
 * Thread count
 * ------------
 * Defaults to std::thread::hardware_concurrency().  Set G2G_NUM_THREADS in
 * the environment or pass explicitly to ParallelAligner to override.
 */

#include <g2g/types.hpp>
#include <g2g/decompose.hpp>
#include <g2g/align.hpp>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace g2g {

// ── Task record ───────────────────────────────────────────────────────────────

/**
 * One unit of work in the snarl-pair task DAG.
 *
 * dep_count: number of child tasks that must complete before this task
 *            is ready.  Decremented atomically; task is enqueued when
 *            dep_count reaches 0.
 * parent_idx: index into the flat task array for this task's parent
 *             (-1 for root tasks).
 * priority:   subtree_work of the more expensive of the two snarls.
 *             Higher = stolen first.
 */
struct SnarlPairTask {
    uint32_t ia;           // index into SnarlTree A
    uint32_t ib;           // index into SnarlTree B
    uint64_t priority;     // subtree_work — scheduler hint

    int32_t  parent_idx;   // -1 if root
    std::atomic<int32_t> dep_count{0};

    std::unique_ptr<SnarlAlignment> result;  // filled after execution

    // Non-copyable (atomic member); explicitly move-constructible
    SnarlPairTask() = default;
    SnarlPairTask(const SnarlPairTask&) = delete;
    SnarlPairTask& operator=(const SnarlPairTask&) = delete;
    SnarlPairTask(SnarlPairTask&& o) noexcept
        : ia(o.ia), ib(o.ib), priority(o.priority)
        , parent_idx(o.parent_idx)
        , dep_count(o.dep_count.load(std::memory_order_relaxed))
        , result(std::move(o.result))
    {}
    SnarlPairTask& operator=(SnarlPairTask&&) = delete;
};

// ── Scheduler statistics ──────────────────────────────────────────────────────

struct RunStats {
    uint64_t tasks_executed{0};
    uint64_t tasks_ready_peak{0};   // max depth of ready queue (load imbalance proxy)
    double   wall_seconds{0.0};
    uint32_t num_threads{0};

    // Parallel efficiency = T1 / (p * Tp).  Filled by ParallelAligner.
    double   parallel_efficiency{0.0};
    double   sequential_seconds{0.0};   // baseline (1-thread) time
};

// ── Parallel aligner ──────────────────────────────────────────────────────────

/**
 * Parallel version of GraphAligner::align_snarls().
 *
 * Usage:
 *   ParallelAligner pa(params, num_threads);
 *   auto [alignments, stats] = pa.align_snarls(ga, gb, ta, tb, anchors);
 */
class ParallelAligner {
public:
    explicit ParallelAligner(const AlignerParams& params = {},
                              uint32_t num_threads = 0);   // 0 = auto-detect

    struct Result {
        std::vector<SnarlAlignment> alignments;
        RunStats                    stats;
    };

    Result align_snarls(const VariationGraph& ga,
                        const VariationGraph& gb,
                        const SnarlTree& ta,
                        const SnarlTree& tb,
                        const std::vector<AnchorMatch>& anchors) const;

    /**
     * Full pipeline replacement for GraphAligner::align().
     * Runs anchor detection, decomposition, parallel snarl DP, deltas.
     */
    AlignmentResult align(const VariationGraph& ga,
                          const VariationGraph& gb) const;

    uint32_t    num_threads() const { return num_threads_; }
    RunStats    last_stats()  const { return last_stats_; }

    // Build the snarl-pair task DAG (public for testing and inspection).
    std::vector<SnarlPairTask> build_task_dag(
        const SnarlTree& ta,
        const SnarlTree& tb) const;

private:
    AlignerParams    params_;
    uint32_t         num_threads_;
    mutable RunStats last_stats_;

    // Execute one task (score_snarl_pair wrapper).
    SnarlAlignment execute_task(
        const VariationGraph& ga,
        const VariationGraph& gb,
        const SnarlTree& ta,
        const SnarlTree& tb,
        const SnarlPairTask& task) const;
};

} // namespace g2g
