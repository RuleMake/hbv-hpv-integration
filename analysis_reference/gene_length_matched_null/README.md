# Gene-length-matched spatial null

These sources separate the restriction to a common supported coordinate set from additional matching on gene overlap and length. Candidate pools are constructed before outcomes are scored. `gene_span_annotation.py` provides the gene-body annotation.

## Source files

| File | Purpose |
|---|---|
| `build_gene_length_matched_design.py` | Builds nested context-, gene-overlap- and gene-length-matched placement pools. |
| `gene_span_annotation.py` | Annotates distinct GENCODE gene-body overlaps and median log gene span. |
| `load_inputs_and_original_null.py` | Loads coordinate and duplication representations and exposes the original spatial estimator. |
| `matching_core.py` | Defines gene-length calipers, coupled candidate draws and statistic summaries. |
| `run_gene_length_sensitivity.py` | Scores the fixed common-support design under each matching condition. |

## Required inputs and layout

Baseline coordinate and variant-geometry tables; GENCODE v19 gene spans and chromosome sizes; genomic-context features; fixed candidate pools; support tables and the design manifest.

The expected layout separates `analysis/`, `baseline/`, `annotation/` and `design/`. Loaders read `outputs/` beneath `baseline/` and `annotation/`. Geometry and annotation are prepared before candidate pools are fixed and outcomes are scored. Hash checks preserve that ordering.

See the release input inventory for data sources and `../source_map.csv` for source locations. Configure the required files, dependencies and output paths in a separate work directory before execution. Release verification checks syntax and file hashes without running these upstream analyses.
