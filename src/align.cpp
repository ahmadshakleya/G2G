#include <g2g/align.hpp>
namespace g2g {
AlignmentResult GraphAligner::align(const VariationGraph&, const VariationGraph&) const { return {}; }
std::vector<SnarlAlignment> GraphAligner::align_snarls(const VariationGraph&, const VariationGraph&, const SnarlTree&, const SnarlTree&, const std::vector<AnchorMatch>&) const { return {}; }
SnarlAlignment GraphAligner::score_snarl_pair(const VariationGraph&, const VariationGraph&, const Snarl&, const Snarl&) const { return {}; }
void GraphAligner::build_deltas(const VariationGraph&, const VariationGraph&, AlignmentResult&, const SnarlTree&, const SnarlTree&) const {}
} // namespace g2g
