#!/bin/bash
# =============================================================================
# setup_cluster.sh — One-time setup on mcmgt01 before running Phase 5 jobs.
#
# Run this interactively (NOT via sbatch) after SSHing into mcmgt01:
#   bash ~/G2G/slurm/scripts/setup_cluster.sh
# =============================================================================

set -euo pipefail

REPO="$HOME/G2G"
cd "$HOME"

echo "======================================================="
echo "  G2G Phase 5 — Cluster Setup"
echo "  Node: $(hostname)"
echo "  Date: $(date)"
echo "======================================================="

# ── 1. Clone or update the repo ───────────────────────────────────────────────
if [[ -d "$REPO/.git" ]]; then
    echo "[setup] Pulling latest from GitHub..."
    git -C "$REPO" pull --ff-only
else
    echo "[setup] Cloning G2G repository..."
    git clone https://github.com/ahmadshakleya/G2G.git "$REPO"
fi

# ── 2. Install Mamba (if not present) ────────────────────────────────────────
if ! command -v mamba &>/dev/null; then
    echo "[setup] Installing Miniforge (mamba)..."
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash "Miniforge3-$(uname)-$(uname -m).sh" -b -p "$HOME/miniforge3"
    eval "$("$HOME/miniforge3/bin/conda" shell.bash hook)"
    rm -f "Miniforge3-$(uname)-$(uname -m).sh"
else
    echo "[setup] mamba already installed: $(which mamba)"
    eval "$(conda shell.bash hook)" 2>/dev/null || true
fi

# ── 3. Create / update the phase4 conda env ──────────────────────────────────
if conda env list | grep -q "^phase4 "; then
    echo "[setup] Updating existing 'phase4' environment..."
    mamba activate phase4
    mamba install -y -c conda-forge cmake gcc gxx libtbb-devel minimap2
else
    echo "[setup] Creating 'phase4' conda environment..."
    mamba create -y -n phase4 -c conda-forge \
        cmake gcc gxx libtbb-devel minimap2 python=3.11
fi

mamba activate phase4
echo "[setup] Compiler:  $(g++ --version | head -1)"
echo "[setup] CMake:     $(cmake --version | head -1)"
echo "[setup] minimap2:  $(minimap2 --version)"
echo "[setup] TBB:       $(find "$CONDA_PREFIX" -name 'libtbb*' 2>/dev/null | head -1)"

# ── 4. Initial Release build ─────────────────────────────────────────────────
echo "[setup] Building G2G (Release)..."
mkdir -p "$REPO/build"
cmake -B "$REPO/build" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER="$(which g++)" \
      -DCMAKE_C_COMPILER="$(which gcc)" \
      -DCMAKE_CXX_FLAGS="-O3 -march=native"
cmake --build "$REPO/build" -j"$(nproc)"

echo ""
echo "[setup] Running existing test suite to confirm baseline..."
cd "$REPO/build"
ctest --output-on-failure -j"$(nproc)" 2>&1 | tail -10

# ── 5. Create output directories ─────────────────────────────────────────────
mkdir -p "$REPO/results"/{exp3,exp4,exp5} "$REPO/slurm/logs"

# ── 6. Check perf availability ───────────────────────────────────────────────
echo ""
echo "[setup] Checking perf availability..."
if perf stat true 2>/dev/null; then
    echo "[setup] perf OK — L3 miss-rate measurements will be collected (Exp5 A1)"
else
    PARANOID=$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo "unknown")
    echo "[setup] WARNING: perf unavailable (perf_event_paranoid=$PARANOID)"
    echo "[setup]   Exp5 A1 will report wall-time only (no L3 miss rate)."
    echo "[setup]   To enable: ask admin to run: sudo sysctl kernel.perf_event_paranoid=2"
fi

# ── 7. Check storage headroom ────────────────────────────────────────────────
echo ""
echo "[setup] Storage check:"
df -h "$HOME" | tail -1
echo "[setup] Results will be written to $REPO/results/"
echo ""

# ── 8. Quick smoke test of the g2g binary ────────────────────────────────────
echo "[setup] Smoke test..."
python3 "$REPO/slurm/synth/gen_scaling_gfa.py" \
    --n_snarls 100 --out_dir /tmp/g2g_smoke --seed 1 2>/dev/null
"$REPO/build/g2g" /tmp/g2g_smoke/ga_100.gfa /tmp/g2g_smoke/gb_100.gfa \
    --gaf /dev/null --threads 4 2>&1 | tail -3
rm -rf /tmp/g2g_smoke
echo "[setup] Smoke test PASSED"

echo ""
echo "======================================================="
echo "  Setup complete.  You can now submit jobs:"
echo ""
echo "  # Submit all three experiments (can run in parallel):"
echo "  cd ~/G2G"
echo "  sbatch slurm/scripts/exp3_throughput.sh"
echo "  sbatch slurm/scripts/exp4_parallel_scaling.sh"
echo "  sbatch slurm/scripts/exp5_ablation.sh"
echo ""
echo "  # Monitor:"
echo "  squeue --me"
echo "  jobstats <job_id>"
echo "======================================================="
