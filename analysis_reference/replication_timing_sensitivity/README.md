# Replication-timing sensitivity

These sources compare a 15-cell-line consensus with HeLa-S3 and HepG2 timing references on a shared set of supported integration loci. The same counting definitions and paired random streams support comparisons across tracks.

## Source files

| File | Purpose |
|---|---|
| `build_matched_candidate_pools.py` | Builds candidate pools and common support across timing references. |
| `baseline_api.py` | Loads the fixed coordinates, breakpoints and event geometry for timing sensitivity. |
| `null_core.py` | Reuses spatial draw/count routines and defines the composition statistic. |
| `run_replication_timing_sensitivity.py` | Computes composition and class-atlas statistics under the three timing references. |
| `stage_replication_timing_tracks.py` | Reconstructs the consensus, HeLa-S3 and HepG2 annotations from cached 50-kb values. |

## Required inputs and layout

Repli-seq component tracks and chromosome sizes; fixed non-timing covariates; coordinate and breakpoint tables; event geometry; saved candidate pools and fixed-design manifests.

The expected layout contains `tracks/`, `baseline/`, `analysis/`, `design/` and `outputs/`. Track preparation precedes candidate-pool construction and estimation. It reconstructs standardized annotations from cached 50-kb values; validation of default versus exact bigWig summaries is a separate check and is not performed by this script. `null_core.py` dynamically loads the shared spatial-null implementation, whose source path must also be configured.

See the release input inventory for data sources and `../source_map.csv` for source locations. Configure the required files, dependencies and output paths in a separate work directory before execution. Release verification checks syntax and file hashes without running these upstream analyses.
