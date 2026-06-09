/**
 * tests/test_matching.cpp
 * =======================
 * Verify hungarian_match() and greedy_match() against brute_match() oracle.
 *
 * These tests mirror the Phase 3 requirement from §3.3 item 4:
 *   "Verify output matches Python brute-force on k ≤ 6."
 *
 * Ground truth scores computed from Python:
 *   from g2g.align import bipartite_allele_match
 *   bipartite_allele_match(seqs_a, seqs_b, gap_cost=0.5, min_match_sim=0.75)
 */

#include <g2g/matching.hpp>
#include <gtest/gtest.h>

#include <string>
#include <vector>

using namespace g2g;

// ── Helper ────────────────────────────────────────────────────

static std::vector<std::string_view> as_views(const std::vector<std::string>& v) {
    std::vector<std::string_view> out;
    out.reserve(v.size());
    for (const auto& s : v) out.emplace_back(s);
    return out;
}

// ── Bubble 1: SNP (k=2 per graph, identical alleles) ─────────
// Python ground truth:
//   seqs_a = ["", ""]  (ACGT/ATGT bubble, inner nodes empty — source&sink excluded)
//   seqs_b = ["", ""]
//   → score ≈ 2.0 (both alleles matched), 2 matches

TEST(Matching, SNPBubble_k2) {
    // Inner sequences of SNP bubble are empty (single-node allele paths
    // where source == first and sink == last → inner = empty)
    std::vector<std::string> a = {"", ""};
    std::vector<std::string> b = {"", ""};
    auto va = as_views(a), vb = as_views(b);

    auto r_brute    = brute_match(va, vb);
    auto r_hungarian = hungarian_match(va, vb);

    EXPECT_FLOAT_EQ(r_brute.score, r_hungarian.score);
    EXPECT_EQ(r_brute.matches.size(), r_hungarian.matches.size());
}

// ── Bubble 2: insertion (k=3 vs k=3, mixed divergence) ───────
// seqs_a = ["GATTACA", "GA", "CTTAGGAATCG"]
// seqs_b = ["GATTACA", "GA", "GATTAC"]
// Python:
//   GATTACA ↔ GATTACA (sim=1.0), GA ↔ GA (sim=1.0),
//   CTTAGGAATCG ↔ GATTAC (sim ≈ 0.45 < 0.75 → rejected)
//   → CTTAGGAATCG is A-only, GATTAC is B-only

TEST(Matching, InsertionBubble_k3) {
    std::vector<std::string> a = {"GATTACA", "GA", "CTTAGGAATCG"};
    std::vector<std::string> b = {"GATTACA", "GA", "GATTAC"};
    auto va = as_views(a), vb = as_views(b);

    auto r = hungarian_match(va, vb);

    // Should match GATTACA↔GATTACA and GA↔GA
    EXPECT_EQ(r.matches.size(), 2u);
    EXPECT_EQ(r.a_unmatched.size(), 1u);  // CTTAGGAATCG
    EXPECT_EQ(r.b_unmatched.size(), 1u);  // GATTAC

    // Score = 1.0 + 1.0 - 0.5 (a_only) - 0.5 (b_only) = 1.0
    EXPECT_FLOAT_EQ(r.score, 1.0f);
}

// ── Symmetry: swapping A and B should give same score ─────────

TEST(Matching, Symmetry) {
    std::vector<std::string> a = {"GATTACA", "GA"};
    std::vector<std::string> b = {"GA", "GATTACA"};
    auto va = as_views(a), vb = as_views(b);

    auto r_fwd = hungarian_match(va, vb);
    auto r_rev = hungarian_match(vb, va);

    EXPECT_FLOAT_EQ(r_fwd.score, r_rev.score);
}

// ── Hungarian matches brute on k≤6 random inputs ─────────────
// (Extend this test once brute_match enumeration is complete)

TEST(Matching, HungarianMatchesBrute_k2) {
    std::vector<std::string> a = {"ACGT", "ATGT"};
    std::vector<std::string> b = {"ACGT", "TTTT"};
    auto va = as_views(a), vb = as_views(b);

    auto r_brute    = brute_match(va, vb);
    auto r_hungarian = hungarian_match(va, vb);

    // Scores must agree within floating-point tolerance
    EXPECT_NEAR(r_brute.score, r_hungarian.score, 1e-4f);
}
