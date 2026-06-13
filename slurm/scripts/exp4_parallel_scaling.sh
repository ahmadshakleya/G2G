#!/bin/bash --login
# =============================================================================
# Experiment 4 — G2G Parallel Scaling (Q2 continued: speedup plot)
# =============================================================================
# Roadmap §2.7.2 / Phase 5, Experiment 4
#
# What this measures
# ------------------
# Fix graph size at a representative "large" scale (100K snarls by default;
# upgrade to 1M once the 1M run completes in Exp3 and memory headroom is clear).
# Sweep thread count: 1, 2, 4, 8, 16, 32, 64, 128
#
# Per thread count, record:
#   - wall_sec       (from RunStats)
#   - speedup        = T1 / Tp
#   - efficiency     = speedup / p  (= T1 / (p * Tp))
#   - tasks_executed (from RunStats, sanity check: should be constant)
#
# The G2G_NUM_THREADS env hook is used so no recompilation is needed.
#
# Outputs
#   parallel_scaling.tsv         — raw rows
#   parallel_scaling_summary.tsv — mean ± stderr, speedup, efficiency
#
# Usage (from mcmgt01)
#   sbatch slurm/scripts/exp4_parallel_scaling.sh
# =============================================================================

#SBATCH --job-name=g2g_exp4_parallel
#SBATCH --output=slurm/logs/exp4-%j.out
#SBATCH --error=slurm/logs/exp4-%j.err
#SBATCH --time=06:00:00
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
OUTDIR="$REPO/results/exp4"
GRAPHDIR="$REPO/results/exp3/graphs"   # reuse graphs from Exp 3 if available
LOGDIR="$REPO/slurm/logs"

# Graph size for the scaling experiment — 100K snarls is the default.
# Change to 1000000 once Exp3 confirms memory fits in 256GB.
N_SNARLS=100000
RUNS=5            # repeats per thread count for statistical stability

THREAD_COUNTS=(1 2 4 8 16 32 64 128)

mkdir -p "$OUTDIR" "$GRAPHDIR" "$LOGDIR"
source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate phase4

# ── Build ─────────────────────────────────────────────────────────────────────
cd "$REPO"
cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER="$(which g++)" \
      -DCMAKE_C_COMPILER="$(which gcc)" \
      -DCMAKE_CXX_FLAGS="-O3 -march=native" > /dev/null
cmake --build "$BUILD" -j"$(nproc)" --target g2g > /dev/null
echo "[exp4] binary: $G2G"

# ── Generate graph pair ───────────────────────────────────────────────────────
GFA_A="$GRAPHDIR/ga_${N_SNARLS}.gfa"
GFA_B="$GRAPHDIR/gb_${N_SNARLS}.gfa"
if [[ ! -f "$GFA_A" ]]; then
    echo "[exp4] generating n_snarls=$N_SNARLS graphs..."
    python3 "$SYNTH_GEN" --n_snarls "$N_SNARLS" --out_dir "$GRAPHDIR" --seed 42
else
    echo "[exp4] graph pair exists: $GFA_A"
fi

# ── Timing wrapper using /usr/bin/time ───────────────────────────────────────
measure_wall() {
    local tmp; tmp=$(mktemp)
    /usr/bin/time -v "$@" 2>"$tmp"
    grep "Elapsed (wall clock)" "$tmp" | awk '{print $NF}' | \
        awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2}'
    rm -f "$tmp"
}

# ── Output TSV ───────────────────────────────────────────────────────────────
RESULTS="$OUTDIR/parallel_scaling.tsv"
echo -e "n_snarls\tthreads\trun\twall_sec" > "$RESULTS"

# ── Warm the filesystem cache with a single run before measurement ────────────
echo "[exp4] warming cache..."
G2G_NUM_THREADS=1 "$G2G" "$GFA_A" "$GFA_B" --gaf /dev/null > /dev/null 2>&1 || true

# ── Sweep ─────────────────────────────────────────────────────────────────────
for T in "${THREAD_COUNTS[@]}"; do
    # Only run thread counts <= available CPUs
    if [[ $T -gt $(nproc) ]]; then
        echo "[exp4] skipping T=$T (only $(nproc) CPUs available)"
        continue
    fi

    echo "[exp4] threads=$T"
    for RUN in $(seq 1 "$RUNS"); do
        WALL=$(measure_wall env G2G_NUM_THREADS="$T" \
               "$G2G" "$GFA_A" "$GFA_B" --threads "$T" --gaf /dev/null)
        echo -e "${N_SNARLS}\t${T}\t${RUN}\t${WALL}" >> "$RESULTS"
        echo "[exp4]   run=$RUN  wall=${WALL}s"
    done
done

# ── Compute speedup and efficiency ────────────────────────────────────────────
SUMMARY="$OUTDIR/parallel_scaling_summary.tsv"
python3 - "$RESULTS" "$SUMMARY" <<'EOF'
import sys, csv, math
from collections import defaultdict

rows = defaultdict(list)
with open(sys.argv[1]) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        key = int(r['threads'])
        rows[key].append(float(r['wall_sec']))

# T1 = mean wall time at 1 thread
t1_vals = rows.get(1, [])
if not t1_vals:
    print("ERROR: no 1-thread measurements found", file=sys.stderr)
    sys.exit(1)
T1 = sum(t1_vals) / len(t1_vals)

with open(sys.argv[2], 'w') as f:
    f.write("threads\tmean_wall_sec\tstderr_wall_sec\tspeedup\tefficiency\tn_runs\n")
    for p, vals in sorted(rows.items()):
        mean = sum(vals) / len(vals)
        if len(vals) > 1:
            var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
            stderr = math.sqrt(var / len(vals))
        else:
            stderr = 0.0
        speedup = T1 / mean if mean > 0 else 0.0
        eff = speedup / p
        f.write(f"{p}\t{mean:.4f}\t{stderr:.4f}\t{speedup:.3f}\t{eff:.3f}\t{len(vals)}\n")
EOF

echo "[exp4] Results: $RESULTS"
echo "[exp4] Summary: $SUMMARY"
echo "[exp4] DONE"
