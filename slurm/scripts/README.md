# G2G Phase 5 — HPC Experiment Scripts (SANDS Lab Testbed)

## Cluster overview

| Resource | Spec |
|---|---|
| Head node | `mcmgt01` (32 CPU, 128 GB RAM) — submit jobs from here |
| Compute target | `mcnode26/27/29` or `synode01-04` (128 CPU, 512 GB RAM) |
| Scheduler | Slurm (`sbatch`, `squeue`, `srun`) |
| Storage | `/home/<user>/` — 112 TB RAIDZ3, visible on all nodes |
| Env manager | Mamba (conda drop-in) |

## First-time setup

SSH into `mcmgt01`, then run the setup script **once**:

```bash
ssh mcmgt01
bash ~/G2G/slurm/scripts/setup_cluster.sh
```

This will:
1. Clone / pull the latest G2G code
2. Install Mamba if needed
3. Create the `phase4` conda environment (gcc, cmake, TBB, minimap2)
4. Build G2G in Release mode and run the full test suite
5. Check `perf` availability for L3 miss-rate measurement (Exp 5 A1)
6. Smoke-test the `g2g` binary

## Submitting experiments

All three Phase 5 experiments can be submitted simultaneously — they share
the same graph files (generated once into `results/exp3/graphs/`) but use
different compute.

```bash
cd ~/G2G
sbatch slurm/scripts/exp3_throughput.sh       # ~8 h, 128 CPU, 256 GB
sbatch slurm/scripts/exp4_parallel_scaling.sh # ~6 h, 128 CPU, 256 GB
sbatch slurm/scripts/exp5_ablation.sh         # ~4 h, 128 CPU, 128 GB

squeue --me
```

## Experiment descriptions

### Experiment 3 — Throughput scaling (`exp3_throughput.sh`)
**Paper question:** Q2 — Is G2G fast?

Sweeps graph size over `n_snarls ∈ {1K, 10K, 100K, 1M}` and measures
wall-clock time and peak RSS for three baselines:

| Baseline | Description |
|---|---|
| `sequential` | G2G, 1 thread — Baseline 3 in paper |
| `parallel_128t` | G2G, 128 threads — the G2G system |
| `minimap2_asm5` | Linearise + minimap2 (skipped for n>100K, infeasible) |

**Outputs:** `results/exp3/throughput_scaling.tsv` + `_summary.tsv`

### Experiment 4 — Parallel scaling (`exp4_parallel_scaling.sh`)
**Paper question:** Q2 continued — speedup & efficiency plot (Figure 3 in paper)

Fixed graph size (100K snarls), thread sweep: `1, 2, 4, 8, 16, 32, 64, 128`.
Uses `G2G_NUM_THREADS` env hook — no recompilation needed per thread count.
5 repeats per thread count for statistical stability.

**Outputs:** `results/exp4/parallel_scaling.tsv` + `_summary.tsv`

To change graph size to 1M (once Exp3 confirms memory fits), edit
`N_SNARLS=1000000` in the script header.

### Experiment 5 — Ablation study (`exp5_ablation.sh`)
**Paper question:** Q3 — Does the system design matter?

Five ablations, each controlled by a compile-time or env flag:

| ID | Ablation | Metric |
|---|---|---|
| A1 | vEB layout vs. BFS layout | Wall time + L3 miss rate (if `perf` available) |
| A2/A3 | Steal-largest vs. steal-random | Wall time (std graph + heavy-tailed graph) |
| A4 | MinimizerIndex vs. naive O(n²) | Anchor detection wall time |
| A5 | Hungarian (LAPJV) vs. brute-force | µs/trial at k = 2,3,5,7,10,20,50 |

**Outputs:** `results/exp5/ablation_results.tsv` + `results/exp5/matching_bench.tsv`

> **Note on perf / L3 miss rate (A1):** The script checks for `perf`
> availability and gracefully skips L3 counters if `perf_event_paranoid` is
> too restrictive. Ask the testbed admin to run
> `sudo sysctl kernel.perf_event_paranoid=2` on the compute node.
> Alternatively, request sudo access per the testbed docs and set it yourself.

## Monitoring jobs

```bash
squeue --me                          # live queue status
jobstats <job_id>                    # resource utilisation summary
jobstats <job_id> -g                 # Grafana dashboard link (rocs/rocs)
tail -f slurm/logs/exp3-<job_id>.out # live log streaming
```

## Output structure

```
results/
├── exp3/
│   ├── graphs/              ← GFA pairs (shared with Exp4 & Exp5)
│   │   ├── ga_1000.gfa
│   │   ├── gb_1000.gfa
│   │   └── ...
│   ├── throughput_scaling.tsv
│   └── throughput_scaling_summary.tsv
├── exp4/
│   ├── parallel_scaling.tsv
│   └── parallel_scaling_summary.tsv
└── exp5/
    ├── ablation_results.tsv
    └── matching_bench.tsv
```

## After experiments complete

Run the figure-generation scripts (to be added in Phase 6):

```bash
python3 python/plot_exp3.py   # throughput vs graph size (log-log)
python3 python/plot_exp4.py   # speedup curve
python3 python/plot_exp5.py   # ablation bar chart + matching bench
```

## Resource limits

The SANDS testbed uses QoS-based limits:
- `normal` QoS: max 2 A100 GPUs (not needed here), no CPU/RAM cap specified
- Interactive sessions: max 4 hours — use `sbatch` for all experiments
- Max job length: 14 days

The scripts request 128 CPUs + 256 GB RAM on `cpu_amd_epyc_7763` nodes
(`mcnode22-29`, `mcnode41-44`, `synode01-04`). If those nodes are busy,
remove the `--constraint` line and Slurm will schedule on any 128-CPU node.
