/**
 * src/sequence.cpp
 * ================
 * C++ port of sequence.py.
 *
 * edit_distance: Wagner-Fischer DP, space O(n).
 * TODO Phase 3 SIMD: for sequences ≤ 16 bytes (most inner-node seqs),
 *   use SSE4.2 PCMPEQB to compare 16 chars in one instruction.
 *   Speedup estimate: 8–16× on x86 (see §2.5.1 finding 1 in the paper).
 *   Gate behind #ifdef __SSE4_2__ so the scalar path remains the fallback.
 */

#include <g2g/sequence.hpp>

#include <algorithm>
#include <cassert>
#include <vector>

namespace g2g {

int edit_distance(std::string_view a, std::string_view b) {
    const size_t m = a.size(), n = b.size();
    if (m == 0) return static_cast<int>(n);
    if (n == 0) return static_cast<int>(m);

    // Rolling single-row DP (space O(n))
    std::vector<int> dp(n + 1);
    for (size_t j = 0; j <= n; ++j) dp[j] = static_cast<int>(j);

    for (size_t i = 1; i <= m; ++i) {
        int prev = static_cast<int>(i) - 1;
        dp[0] = static_cast<int>(i);
        for (size_t j = 1; j <= n; ++j) {
            int cur = std::min({
                prev + (a[i-1] == b[j-1] ? 0 : 1),
                dp[j]   + 1,
                dp[j-1] + 1
            });
            prev  = dp[j];
            dp[j] = cur;
        }
    }
    return dp[n];
}

void path_sequence(const VariationGraph& vg,
                   std::span<const NodeId> path,
                   std::string& out) {
    for (NodeId n : path) {
        auto sv = vg.seq(n);
        out.append(sv.data(), sv.size());
    }
}

void allele_sequences(const VariationGraph& vg,
                      const Snarl& s,
                      std::vector<std::string>& out) {
    // Each allele is a slice [alleles_begin + i*stride, alleles_begin + i*stride + stride)
    // of path_node_ids; inner nodes only (exclude first and last = source/sink).
    // TODO: this placeholder assumes allele paths are stored contiguously in
    //       VariationGraph::path_node_ids.  Fill in once VariationGraph::add_path
    //       is complete and snarl indexing is wired up.
    out.clear();
    out.resize(s.allele_count());
    // Real implementation: iterate s.alleles_begin..s.alleles_end,
    // skip first and last NodeId of each allele path, concatenate the rest.
}

} // namespace g2g
