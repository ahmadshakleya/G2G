#include <g2g/io.hpp>
#include <stdexcept>
namespace g2g {
VariationGraph read_gfa(const std::filesystem::path&) { throw std::runtime_error("GFA reader not yet implemented"); }
VariationGraph read_gfa(std::istream&, const std::string&) { throw std::runtime_error("GFA reader not yet implemented"); }
void write_alignment_gaf(std::ostream&, const AlignmentResult&, const VariationGraph&, const VariationGraph&, const SnarlTree&, const SnarlTree&) {}
void write_delta_vcf(std::ostream&, const AlignmentResult&, const VariationGraph&, const SnarlTree&) {}
} // namespace g2g
