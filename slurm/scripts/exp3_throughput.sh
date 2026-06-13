#!/bin/bash --login
# =============================================================================
# Experiment 3 — G2G Throughput Scaling (Q2: Is G2G fast?)
# =============================================================================
# Roadmap §2.7.2 / Phase 5, Experiment 3
#
# What this measures
# ------------------
# Wall-clock time and peak RSS as graph size grows:
#   n_snarls ∈ {1K, 10K, 100K, 1M}
#
# Baselines compared
#   B0: G2G sequential (same algorithm, 1 thread)  → Baseline 3 in paper
#   B1: G2G parallel   (all cores)                 → G2G system
#   B2: linearise + minimap2 asm5                  → Baseline 1 in paper
#       (extract all haplotype sequences, all-vs-all minimap2)
#
# Outputs (written to $OUTDIR)
#   throughput_scaling.tsv        — one row per (n_snarls, baseline, run)
#   throughput_scaling_summary.tsv — mean ± stderr per (n_snarls, baseline)
#
# Usage (from mcmgt01 head node)
#   sbatch slurm/scripts/exp3_throughput.sh
# =============================================================================

#SBATCH --job-name=g2g_exp3_throughput
#SBATCH --output=slurm/logs/exp3-%j.out
#SBATCH --error=slurm/logs/exp3-%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=256GB
#SBATCH --constraint=cpu_amd_epyc_7763

# ── Environment ───────────────────────────────────────────────────────────────
set -euo pipefail

REPO="$HOME/Research/G2G/G2G_Project_Summer_2026"
BUILD="$REPO/build"
G2G="$BUILD/g2g"
SYNTH_GEN="$REPO/slurm/synth/gen_scaling_gfa.py"
OUTDIR="$REPO/results/exp3"
GRAPHDIR="$OUTDIR/graphs"
LOGDIR="$REPO/slurm/logs"
RUNS=3          # repeat each measurement this many times for stability

mkdir -p "$OUTDIR" "$GRAPHDIR" "$LOGDIR"

# Activate the conda environment that has cmake, gcc, tbb, minimap2
source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate phase4

# Rebuild in Release mode with the cluster compiler
cd "$REPO"
cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER="$(which g++)" \
      -DCMAKE_C_COMPILER="$(which gcc)" -DCMAKE_VERBOSE_MAKEFILE=OFF \
      -DCMAKE_CXX_FLAGS="-O3 -march=native" > /dev/null
cmake --build "$BUILD" -j"$(nproc)" --target g2g > /dev/null
echo "[exp3] binary: $G2G"

# ── Helper: time a command, capture wall seconds and peak RSS (kB) ────────────
# Writes to a temp file because bash can't easily capture /usr/bin/time stderr.
measure() {
    local label="$1"; shift
    local tmp
    tmp=$(mktemp)
    /usr/bin/time -v "$@" 2>"$tmp"
    local wall
    wall=$(grep "Elapsed (wall clock)" "$tmp" | awk '{print $NF}' | \
           awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2}')
    local rss
    rss=$(grep "Maximum resident set size" "$tmp" | awk '{print $NF}')
    rm -f "$tmp"
    echo "$label $wall $rss"
}

# ── Output TSV header ─────────────────────────────────────────────────────────
RESULTS="$OUTDIR/throughput_scaling.tsv"
echo -e "n_snarls\tbaseline\trun\twall_sec\tpeak_rss_kb" > "$RESULTS"

# ── Sizes to sweep ────────────────────────────────────────────────────────────
SIZES=(1000 10000 100000 1000000)

for N in "${SIZES[@]}"; do
    echo "[exp3] ===== n_snarls=$N ====="

    # Generate GFA pair once
    GFA_A="$GRAPHDIR/ga_${N}.gfa"
    GFA_B="$GRAPHDIR/gb_${N}.gfa"
    if [[ ! -f "$GFA_A" ]]; then
        echo "[exp3]   generating graphs..."
        python3 "$SYNTH_GEN" --n_snarls "$N" --out_dir "$GRAPHDIR" --seed 42
    else
        echo "[exp3]   graphs exist, skipping generation"
    fi

    for RUN in $(seq 1 "$RUNS"); do
        echo "[exp3]   run=$RUN"

        # ── Baseline 0: G2G sequential ────────────────────────────
        echo "[exp3]     B0: sequential"
        read -r _ WALL RSS < <(measure "b0" "$G2G" "$GFA_A" "$GFA_B" \
            --sequential --gaf /dev/null)
        echo -e "${N}\tsequential\t${RUN}\t${WALL}\t${RSS}" >> "$RESULTS"

        # ── Baseline 1: G2G parallel (all cores) ─────────────────
        echo "[exp3]     B1: parallel ($(nproc) threads)"
        read -r _ WALL RSS < <(measure "b1" "$G2G" "$GFA_A" "$GFA_B" \
            --threads "$(nproc)" --gaf /dev/null)
        echo -e "${N}\tparallel_$(nproc)t\t${RUN}\t${WALL}\t${RSS}" >> "$RESULTS"

        # ── Baseline 2: linearise + minimap2 asm5 ────────────────
        # Extract haplotype sequences as FASTA from both graphs,
        # then run minimap2 all-vs-all. Skip for n>100K (infeasible).
        if [[ $N -le 100000 ]]; then
            echo "[exp3]     B2: minimap2 asm5"
            FASTA_A="$GRAPHDIR/seqs_a_${N}.fa"
            FASTA_B="$GRAPHDIR/seqs_b_${N}.fa"

            # Extract path sequences using awk on GFA P-lines
            if [[ ! -f "$FASTA_A" ]]; then
                python3 - "$GFA_A" "$FASTA_A" <<'EOF'
import sys, re

def rc(s):
    t = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(t.get(c,c) for c in reversed(s))

gfa, out = sys.argv[1], sys.argv[2]
seqs = {}
paths = {}
with open(gfa) as f:
    for line in f:
        if line.startswith('S\t'):
            parts = line.rstrip().split('\t')
            seqs[parts[1]] = parts[2]
        elif line.startswith('P\t'):
            parts = line.rstrip().split('\t')
            name = parts[1]
            nodes = [n.rstrip('+-') for n in parts[2].split(',')]
            orients = [n[-1] for n in parts[2].split(',')]
            seq = ''
            for n, o in zip(nodes, orients):
                s = seqs.get(n, '')
                seq += s if o == '+' else rc(s)
            paths[name] = seq

with open(out, 'w') as f:
    for name, seq in paths.items():
        f.write(f'>{name}\n{seq}\n')
EOF
                python3 - "$GFA_B" "$FASTA_B" <<'EOF'
import sys, re

def rc(s):
    t = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(t.get(c,c) for c in reversed(s))

gfa, out = sys.argv[1], sys.argv[2]
seqs = {}
paths = {}
with open(gfa) as f:
    for line in f:
        if line.startswith('S\t'):
            parts = line.rstrip().split('\t')
            seqs[parts[1]] = parts[2]
        elif line.startswith('P\t'):
            parts = line.rstrip().split('\t')
            name = parts[1]
            nodes = [n.rstrip('+-') for n in parts[2].split(',')]
            orients = [n[-1] for n in parts[2].split(',')]
            seq = ''
            for n, o in zip(nodes, orients):
                s = seqs.get(n, '')
                seq += s if o == '+' else rc(s)
            paths[name] = seq

with open(out, 'w') as f:
    for name, seq in paths.items():
        f.write(f'>{name}\n{seq}\n')
EOF
            fi

            read -r _ WALL RSS < <(measure "b2" minimap2 -x asm5 \
                -t "$(nproc)" --secondary=no \
                "$FASTA_A" "$FASTA_B" -o /dev/null)
            echo -e "${N}\tminimap2_asm5\t${RUN}\t${WALL}\t${RSS}" >> "$RESULTS"
        else
            echo "[exp3]     B2: minimap2 SKIPPED for n=$N (infeasible)"
        fi
    done
done

# ── Compute mean ± stderr per (n_snarls, baseline) ───────────────────────────
SUMMARY="$OUTDIR/throughput_scaling_summary.tsv"
python3 - "$RESULTS" "$SUMMARY" <<'EOF'
import sys, csv, math
from collections import defaultdict

rows = defaultdict(list)
with open(sys.argv[1]) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        key = (int(r['n_snarls']), r['baseline'])
        rows[key].append(float(r['wall_sec']))

with open(sys.argv[2], 'w') as f:
    f.write("n_snarls\tbaseline\tmean_wall_sec\tstderr_wall_sec\tn_runs\n")
    for (n, bl), vals in sorted(rows.items()):
        mean = sum(vals) / len(vals)
        if len(vals) > 1:
            var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
            stderr = math.sqrt(var / len(vals))
        else:
            stderr = 0.0
        f.write(f"{n}\t{bl}\t{mean:.4f}\t{stderr:.4f}\t{len(vals)}\n")
EOF

echo "[exp3] Results written to $RESULTS"
echo "[exp3] Summary written to $SUMMARY"
echo "[exp3] DONE"
