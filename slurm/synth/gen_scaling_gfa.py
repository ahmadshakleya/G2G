#!/usr/bin/env python3
"""
gen_scaling_gfa.py — Generate a pair of synthetic GFA files for G2G scaling experiments.

Produces GA and GB as variation graphs with a linear backbone of n_snarls snarls.
Each snarl has k alleles drawn from a Zipf-like distribution (heavy-tailed, matching
the HPRC snarl statistics from Phase 2 profiling).

Usage:
    python gen_scaling_gfa.py --n_snarls 1000 --out_dir /path/to/out --seed 42

Outputs:
    <out_dir>/ga_<n_snarls>.gfa
    <out_dir>/gb_<n_snarls>.gfa
    <out_dir>/truth_<n_snarls>.tsv   # ground-truth shared/a_only/b_only labels
"""

import argparse
import math
import os
import random
import sys

BASES = "ACGT"

def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))

def zipf_k(rng: random.Random, k_min: int = 2, k_max: int = 50, alpha: float = 2.0) -> int:
    """
    Sample allele count k from a Zipf distribution truncated to [k_min, k_max].
    alpha=2.0 gives a realistic heavy tail matching HPRC snarl size data.
    """
    weights = [1.0 / (i ** alpha) for i in range(k_min, k_max + 1)]
    total = sum(weights)
    r = rng.random() * total
    cumulative = 0.0
    for i, w in enumerate(weights):
        cumulative += w
        if r <= cumulative:
            return k_min + i
    return k_max

def write_gfa(path: str, nodes: dict, edges: list, paths: list):
    """Write a minimal GFA v1.1 file."""
    with open(path, "w") as f:
        f.write("H\tVN:Z:1.1\n")
        for name, seq in nodes.items():
            f.write(f"S\t{name}\t{seq}\n")
        for src, dst in edges:
            f.write(f"L\t{src}\t+\t{dst}\t+\t0M\n")
        for path_name, node_list in paths:
            cigar_list = ",".join(f"{n}+" for n in node_list)
            f.write(f"P\t{path_name}\t{cigar_list}\t{'*,' * (len(node_list) - 1)}*\n")

def build_graphs(n_snarls: int, rng: random.Random, anchor_seq_len: int = 20,
                 inner_seq_len: int = 50):
    """
    Build a pair of variation graphs GA and GB with planted variants.

    Returns:
        ga_nodes, ga_edges, ga_paths,
        gb_nodes, gb_edges, gb_paths,
        truth_rows  (list of dicts for TSV)
    """
    ga_nodes, gb_nodes = {}, {}
    ga_edges, gb_edges = [], []
    ga_path_nodes, gb_path_nodes = [], []
    truth_rows = []

    # ── Backbone anchors ────────────────────────────────────────
    anc_seqs = [random_dna(anchor_seq_len, rng) for _ in range(n_snarls + 1)]
    for i, seq in enumerate(anc_seqs):
        name = f"anc_{i}"
        ga_nodes[name] = seq
        gb_nodes[name] = seq

    ga_path_nodes.append("anc_0")
    gb_path_nodes.append("anc_0")

    # ── Snarl bubbles ────────────────────────────────────────────
    for s in range(n_snarls):
        src = f"anc_{s}"
        snk = f"anc_{s + 1}"
        k = zipf_k(rng)

        n_shared = max(1, k // 2)
        n_a_only = (k - n_shared) // 2 + ((k - n_shared) % 2)
        n_b_only = (k - n_shared) // 2

        # Shared alleles
        for j in range(n_shared):
            seq = random_dna(inner_seq_len, rng)
            na = f"s{s}_sh{j}_a"
            nb = f"s{s}_sh{j}_b"
            ga_nodes[na] = seq
            gb_nodes[nb] = seq
            ga_edges += [(src, na), (na, snk)]
            gb_edges += [(src, nb), (nb, snk)]
            truth_rows.append({"snarl": s, "type": "shared", "seq": seq,
                                "node_a": na, "node_b": nb})

        # A-only alleles
        for j in range(n_a_only):
            seq = random_dna(inner_seq_len, rng)
            na = f"s{s}_ao{j}"
            ga_nodes[na] = seq
            ga_edges += [(src, na), (na, snk)]
            truth_rows.append({"snarl": s, "type": "a_only", "seq": seq,
                                "node_a": na, "node_b": ""})

        # B-only alleles
        for j in range(n_b_only):
            seq = random_dna(inner_seq_len, rng)
            nb = f"s{s}_bo{j}"
            gb_nodes[nb] = seq
            gb_edges += [(src, nb), (nb, snk)]
            truth_rows.append({"snarl": s, "type": "b_only", "seq": seq,
                                "node_a": "", "node_b": nb})

        # Haplotype paths: use first shared allele
        ga_path_nodes += [f"s{s}_sh0_a", snk]
        gb_path_nodes += [f"s{s}_sh0_b", snk]

    ga_paths = [("hap_A0", ga_path_nodes)]
    gb_paths = [("hap_B0", gb_path_nodes)]

    return (ga_nodes, ga_edges, ga_paths,
            gb_nodes, gb_edges, gb_paths,
            truth_rows)


def main():
    ap = argparse.ArgumentParser(description="Generate GFA pair for G2G scaling experiments")
    ap.add_argument("--n_snarls", type=int, required=True,
                    help="Number of snarls in each graph")
    ap.add_argument("--out_dir", default=".", help="Output directory")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--anchor_len", type=int, default=20)
    ap.add_argument("--inner_len", type=int, default=50)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)

    print(f"[gen] building n_snarls={args.n_snarls} seed={args.seed}", flush=True)
    (ga_nodes, ga_edges, ga_paths,
     gb_nodes, gb_edges, gb_paths,
     truth_rows) = build_graphs(args.n_snarls, rng,
                                args.anchor_len, args.inner_len)

    n = args.n_snarls
    gfa_a = os.path.join(args.out_dir, f"ga_{n}.gfa")
    gfa_b = os.path.join(args.out_dir, f"gb_{n}.gfa")
    truth  = os.path.join(args.out_dir, f"truth_{n}.tsv")

    write_gfa(gfa_a, ga_nodes, ga_edges, ga_paths)
    write_gfa(gfa_b, gb_nodes, gb_edges, gb_paths)

    with open(truth, "w") as f:
        f.write("snarl\ttype\tseq\tnode_a\tnode_b\n")
        for row in truth_rows:
            f.write(f"{row['snarl']}\t{row['type']}\t{row['seq']}\t"
                    f"{row['node_a']}\t{row['node_b']}\n")

    total_nodes_a = len(ga_nodes)
    total_nodes_b = len(gb_nodes)
    print(f"[gen] GA nodes={total_nodes_a}  GB nodes={total_nodes_b}")
    print(f"[gen] truth rows={len(truth_rows)}")
    print(f"[gen] wrote {gfa_a}")
    print(f"[gen] wrote {gfa_b}")
    print(f"[gen] wrote {truth}")


if __name__ == "__main__":
    main()
