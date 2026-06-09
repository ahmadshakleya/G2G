#include <g2g/anchor.hpp>
namespace g2g {
MinimizerIndex::MinimizerIndex(const VariationGraph& vg, const AnchorParams& p) : vg_(vg), params_(p) { build(); }
std::vector<NodeId> MinimizerIndex::query(std::string_view) const { return {}; }
void MinimizerIndex::build() {}
std::vector<uint32_t> MinimizerIndex::minimizers(std::string_view) const { return {}; }
uint32_t MinimizerIndex::hash_kmer(std::string_view) { return 0; }
std::vector<AnchorMatch> find_anchors(const VariationGraph&, const VariationGraph&, const AnchorParams&) { return {}; }
} // namespace g2g
