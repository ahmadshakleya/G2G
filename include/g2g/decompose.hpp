#pragma once
/**
 * include/g2g/decompose.hpp
 * =========================
 * Snarl decomposition and snarl-tree layout.
 *
 * SnarlTree stores snarls in van Emde Boas (vEB) recursive order
 * (§2.5.2 of the paper) so that bottom-up DP traversal exhibits
 * cache-oblivious spatial locality.
 *
 * Two back-ends are provided at compile time:
 *   G2G_HAVE_HANDLEGRAPH  → delegate to libsnarls from vg toolkit.
 *                            Handles nested snarls, unary snarls,
 *                            cycles (via DAG projection).
 *   (default)             → built-in dominator approximation
 *                            (port of decompose.py), suitable for
 *                            DAGs only.  Used for testing.
 *
 * The vEB layout is applied by SnarlTree::build_veb_order() regardless
 * of which back-end produced the snarls.
 */

#include <g2g/types.hpp>
#include <g2g/graph.hpp>
#include <optional>   // add this line
#include <cstdint>
#include <vector>

namespace g2g {

// ── Snarl tree ────────────────────────────────────────────────

/**
 * Flat array of Snarl objects in van Emde Boas order.
 * The tree structure is encoded in Snarl::children_begin/end —
 * all children of snarl[i] are stored at snarls[children_begin..children_end).
 *
 * vEB ordering guarantees that for any subtree of height h,
 * all its nodes fit in O(N^{1/2}) cache lines, giving an
 * O(N log_B N / B) cache complexity for bottom-up traversal
 * (optimal for comparison-based tree algorithms).
 */
struct SnarlTree {
    std::vector<Snarl>    snarls;     // vEB-ordered
    std::vector<uint32_t> roots;      // indices of top-level snarls
    uint32_t              depth{0};   // tree height

    uint32_t size() const { return static_cast<uint32_t>(snarls.size()); }
};

// ── Decomposer ────────────────────────────────────────────────

class SnarlDecomposer {
public:
    explicit SnarlDecomposer(const VariationGraph& vg) : vg_(vg) {}

    /**
     * Decompose the graph into a SnarlTree in vEB order.
     * Raises std::runtime_error if the graph contains cycles
     * and G2G_HAVE_HANDLEGRAPH is not defined.
     */
    SnarlTree decompose() const;

private:
    const VariationGraph& vg_;

    // Convert a list of snarls (DFS order) to vEB order in-place.
    static void apply_veb_order(SnarlTree& tree);

    // Compute subtree_work for each snarl (sum of k^3 over all descendants).
    static void annotate_work(SnarlTree& tree);
};

// ── Topology helpers ──────────────────────────────────────────

/** Topological sort of a DAG; throws if cycles detected. */
std::vector<NodeId> topological_sort(const VariationGraph& vg);

/**
 * Cache the topological order in the graph for reuse.
 * Phase 3 fix: Python prototype called topological_sort 1.8M times;
 * C++ implementation computes it once per graph and stores it here.
 */
struct TopoCache {
    const VariationGraph*  vg{nullptr};
    std::vector<NodeId>    order;

    void build(const VariationGraph& g);
    const std::vector<NodeId>& get(const VariationGraph& g);
};

} // namespace g2g
