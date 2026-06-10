#include <g2g/decompose.hpp>
#include <algorithm>
#include <functional>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>

namespace g2g {

// ── Topological sort (Kahn's algorithm) ──────────────────────────────────────

std::vector<NodeId> topological_sort(const VariationGraph& vg) {
    uint32_t n = vg.num_nodes();
    std::vector<uint32_t> in_deg(n);
    for (NodeId u = 0; u < n; ++u)
        for (NodeId v : vg.successors(u))
            in_deg[v]++;

    std::queue<NodeId> q;
    for (NodeId u = 0; u < n; ++u)
        if (in_deg[u] == 0) q.push(u);

    std::vector<NodeId> order;
    order.reserve(n);
    while (!q.empty()) {
        NodeId u = q.front(); q.pop();
        order.push_back(u);
        for (NodeId v : vg.successors(u))
            if (--in_deg[v] == 0) q.push(v);
    }
    if (order.size() != n)
        throw std::runtime_error("Graph contains cycles — DAG required");
    return order;
}

void TopoCache::build(const VariationGraph& g) {
    vg    = &g;
    order = topological_sort(g);
}

const std::vector<NodeId>& TopoCache::get(const VariationGraph& g) {
    if (vg != &g) build(g);
    return order;
}

// ── find_sink ─────────────────────────────────────────────────────────────────

static std::optional<NodeId> find_sink(
        NodeId source,
        const VariationGraph& vg,
        const std::vector<NodeId>& topo,
        const std::vector<uint32_t>& /*topo_pos*/) {

    auto succs = vg.successors(source);
    if (succs.size() < 2) return std::nullopt;

    std::vector<std::vector<bool>> reachable(succs.size(),
                                              std::vector<bool>(vg.num_nodes(), false));
    for (size_t si = 0; si < succs.size(); ++si) {
        std::queue<NodeId> q;
        q.push(succs[si]);
        reachable[si][succs[si]] = true;
        while (!q.empty()) {
            NodeId u = q.front(); q.pop();
            for (NodeId v : vg.successors(u))
                if (!reachable[si][v]) { reachable[si][v] = true; q.push(v); }
        }
    }

    for (NodeId n : topo) {
        if (n == source) continue;
        bool in_all = true;
        for (size_t si = 0; si < succs.size(); ++si)
            if (!reachable[si][n]) { in_all = false; break; }
        if (in_all) return n;
    }
    return std::nullopt;
}

// ── SnarlDecomposer::decompose() ──────────────────────────────────────────────

SnarlTree SnarlDecomposer::decompose() const {
    std::vector<NodeId> topo = topological_sort(vg_);
    std::vector<uint32_t> topo_pos(vg_.num_nodes());
    for (uint32_t i = 0; i < topo.size(); ++i)
        topo_pos[topo[i]] = i;

    SnarlTree tree;
    std::vector<bool> visited(vg_.num_nodes(), false);

    // Detect all candidate snarls and store their allele paths temporarily
    struct RawSnarl {
        Snarl snarl;
        std::vector<std::vector<NodeId>> alleles;
    };
    std::vector<RawSnarl> raw;

    for (NodeId s : topo) {
        if (visited[s]) continue;
        if (vg_.out_degree(s) < 2) continue;

        auto sink_opt = find_sink(s, vg_, topo, topo_pos);
        if (!sink_opt) continue;
        NodeId t = *sink_opt;

        std::vector<std::vector<NodeId>> alleles;
        std::function<void(NodeId, std::vector<NodeId>&)> dfs =
            [&](NodeId u, std::vector<NodeId>& path) {
                if (u == t) { alleles.push_back(path); return; }
                for (NodeId v : vg_.successors(u)) {
                    if (topo_pos[v] > topo_pos[t]) continue;
                    path.push_back(v);
                    dfs(v, path);
                    path.pop_back();
                }
            };
        std::vector<NodeId> path = {s};
        dfs(s, path);

        if (alleles.size() < 2) continue;

        Snarl snarl{};
        snarl.source       = s;
        snarl.sink         = t;
        snarl.alleles_begin = 0;
        snarl.alleles_end   = static_cast<uint32_t>(alleles.size());
        snarl.children_begin = 0;
        snarl.children_end   = 0;
        snarl.subtree_work   = 0;

        raw.push_back({snarl, alleles});
        visited[s] = true;
    }

    // Collect interior nodes for every raw snarl
    // Interior = nodes strictly between source and sink (excluding both)
    std::vector<std::set<NodeId>> interiors(raw.size());
    for (size_t i = 0; i < raw.size(); ++i) {
        for (const auto& allele_path : raw[i].alleles)
            for (size_t k = 1; k + 1 < allele_path.size(); ++k)
                interiors[i].insert(allele_path[k]);
    }

    // A snarl is top-level if its source is NOT an interior node of any other snarl
    std::set<NodeId> all_interiors;
    for (const auto& iset : interiors)
        all_interiors.insert(iset.begin(), iset.end());

    // Build flat allele-path store for top-level snarls (P3-6).
    // For each accepted snarl, append its allele paths into the flat arrays
    // and set alleles_begin/end on the Snarl to index into those arrays.
    uint32_t flat_allele_idx = 0;
    for (auto& r : raw) {
        if (all_interiors.find(r.snarl.source) != all_interiors.end()) continue;

        r.snarl.alleles_begin = flat_allele_idx;
        r.snarl.alleles_end   = flat_allele_idx + static_cast<uint32_t>(r.alleles.size());

        for (const auto& path : r.alleles) {
            tree.allele_path_offsets.push_back(
                static_cast<uint32_t>(tree.allele_paths_data.size()));
            tree.allele_path_lengths.push_back(
                static_cast<uint32_t>(path.size()));
            tree.allele_paths_data.insert(tree.allele_paths_data.end(),
                                          path.begin(), path.end());
        }
        flat_allele_idx += static_cast<uint32_t>(r.alleles.size());
        tree.snarls.push_back(r.snarl);
    }

    tree.roots.resize(tree.snarls.size());
    std::iota(tree.roots.begin(), tree.roots.end(), 0u);
    tree.depth = 1;

    return tree;
}

void SnarlDecomposer::apply_veb_order(SnarlTree&) {}
void SnarlDecomposer::annotate_work(SnarlTree&) {}

} // namespace g2g