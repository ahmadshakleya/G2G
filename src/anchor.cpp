#include <g2g/anchor.hpp>
#include <g2g/sequence.hpp>

namespace g2g {

// ── MinimizerIndex stubs (Phase 3b — implement after naïve version works) ──

MinimizerIndex::MinimizerIndex(const VariationGraph& vg, const AnchorParams& p)
    : vg_(vg), params_(p) { build(); }

std::vector<NodeId> MinimizerIndex::query(std::string_view) const { return {}; }
void MinimizerIndex::build() {}
std::vector<uint32_t> MinimizerIndex::minimizers(std::string_view) const { return {}; }
uint32_t MinimizerIndex::hash_kmer(std::string_view) { return 0; }

// ── find_anchors — naïve O(|VA|·|VB|) implementation ──────────────────────
// Direct port of _find_anchors() from align.py.
// Backbone = nodes with in_degree ≤ 1 AND out_degree ≤ 1 (chain nodes).
// Greedy best-first: each node used at most once.
// Replace inner loop with MinimizerIndex once this passes all tests.

std::vector<AnchorMatch> find_anchors(const VariationGraph& ga,
                                       const VariationGraph& gb,
                                       const AnchorParams& params) {
    // Collect backbone nodes from each graph
    std::vector<NodeId> backbone_a, backbone_b;
    for (NodeId n = 0; n < ga.num_nodes(); ++n)
        if (ga.in_degree(n) <= 1 && ga.out_degree(n) <= 1)
            backbone_a.push_back(n);
    for (NodeId n = 0; n < gb.num_nodes(); ++n)
        if (gb.in_degree(n) <= 1 && gb.out_degree(n) <= 1)
            backbone_b.push_back(n);

    std::vector<AnchorMatch> anchors;
    std::vector<bool> used_b(gb.num_nodes(), false);

    for (NodeId na : backbone_a) {
        float   best_sim = params.min_sim - 1e-9f;
        NodeId  best_nb  = kInvalidNode;

        for (NodeId nb : backbone_b) {
            if (used_b[nb]) continue;
            float s = seq_similarity(ga.seq(na), gb.seq(nb));
            if (s > best_sim) { best_sim = s; best_nb = nb; }
        }

        if (best_nb != kInvalidNode) {
            anchors.push_back({na, best_nb, best_sim});
            used_b[best_nb] = true;
        }
    }
    return anchors;
}

} // namespace g2g