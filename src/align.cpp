#include <g2g/align.hpp>
#include <g2g/anchor.hpp>
#include <g2g/sequence.hpp>

namespace g2g {

AlignmentResult GraphAligner::align(const VariationGraph& ga,
                                     const VariationGraph& gb) const {
    AlignmentResult result;

    // Step 1 — anchor detection
    result.anchors = find_anchors(ga, gb, params_.anchor);

    // Step 2 — snarl decomposition (stub returns empty trees for now)
    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    // Step 3 — snarl alignment (stub returns empty for now)
    result.snarl_alignments = align_snarls(ga, gb, ta, tb, result.anchors);

    // Step 4 — total score: anchor contributions + snarl scores
    for (const auto& a : result.anchors)
        result.total_score += params_.alpha * a.sim;
    for (const auto& sa : result.snarl_alignments)
        result.total_score += sa.score;

    // Step 5 — delta graph (stub for now)
    build_deltas(ga, gb, result, ta, tb);

    return result;
}

std::vector<SnarlAlignment> GraphAligner::align_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const std::vector<AnchorMatch>&) const {

    std::vector<SnarlAlignment> results;
    std::vector<bool> used_b(tb.snarls.size(), false);

    for (uint32_t ia = 0; ia < ta.snarls.size(); ++ia) {
        float best_score = 0.f;
        int   best_ib    = -1;

        for (uint32_t ib = 0; ib < tb.snarls.size(); ++ib) {
            if (used_b[ib]) continue;
            SnarlAlignment aln = score_snarl_pair(ga, gb,
                                                   ta.snarls[ia],
                                                   tb.snarls[ib]);
            if (aln.score > best_score) {
                best_score = aln.score;
                best_ib    = static_cast<int>(ib);
            }
        }
        if (best_ib >= 0) {
            results.push_back(score_snarl_pair(ga, gb,
                                                ta.snarls[ia],
                                                tb.snarls[best_ib]));
            used_b[best_ib] = true;
        }
    }
    return results;
}

SnarlAlignment GraphAligner::score_snarl_pair(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const Snarl& sa,
    const Snarl& sb) const {

    SnarlAlignment aln{};
    aln.snarl_a_idx = 0;  // index lookup not needed for scoring
    aln.snarl_b_idx = 0;

    // W_seq: average similarity of source and sink boundary nodes
    float src_sim  = seq_similarity(ga.seq(sa.source), gb.seq(sb.source));
    float sink_sim = seq_similarity(ga.seq(sa.sink),   gb.seq(sb.sink));
    aln.w_seq  = params_.alpha * (src_sim + sink_sim) / 2.f;

    // W_topo: 1.0 — all detected snarls are simple bubbles
    aln.w_topo = params_.beta * 1.f;

    // W_struct: skipped until allele path storage is wired up (Phase 3 P3-6)
    aln.w_struct = 0.f;

    aln.score = aln.w_seq + aln.w_topo;
    return aln;
}

void GraphAligner::build_deltas(
    const VariationGraph&,
    const VariationGraph&,
    AlignmentResult&,
    const SnarlTree&,
    const SnarlTree&) const {
    // TODO
}

} // namespace g2g