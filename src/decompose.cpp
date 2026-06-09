#include <g2g/decompose.hpp>
namespace g2g {
SnarlTree SnarlDecomposer::decompose() const { return {}; }
void SnarlDecomposer::apply_veb_order(SnarlTree&) {}
void SnarlDecomposer::annotate_work(SnarlTree&) {}
std::vector<NodeId> topological_sort(const VariationGraph&) { return {}; }
void TopoCache::build(const VariationGraph&) {}
const std::vector<NodeId>& TopoCache::get(const VariationGraph& g) { build(g); return order; }
} // namespace g2g
