# Analysis source code

These 79 Python files describe the statistical analyses behind the manuscript: structural-variant counting, matched genomic comparisons, donor-level estimates and supporting molecular analyses. They are grouped by scientific purpose. Each directory lists its required inputs and dependencies.

## Using the sources

These programs require external catalogues, prepared annotation tables, candidate pools and saved comparison outputs in the layout specified by each analysis. `source_map.csv` records released paths and source checksums; the input inventory identifies available data sources.

Prepare a separate analysis work directory, then:

1. Obtain the required catalogues and annotations under their access conditions.
2. Arrange the `data/`, `outputs/`, `baseline/`, `design/`, `tracks/`, `annotation/` and `shared/` directories as described in the relevant README, or update the path constants consistently.
3. Configure the relative `source_project` or `upstream_project` input root used by the selected program, together with any dynamic module loaders. Adapt Windows separators if needed, and keep source inputs separate from generated outputs.
4. Install the required dependencies and check the input manifests and fixed-design hashes before running the analysis.

Some scripts read data or write results at import time. Release verification checks their file hashes and Python syntax without importing them. The figure-reproduction commands elsewhere in the package run separately from saved numerical results; they do not rerun these upstream analyses.

## Dependencies

The files were checked with Python 3.12. `requirements-analysis.txt` lists imported third-party packages, including pandas' parquet and Excel readers. It is a dependency list rather than a tested environment lock for every upstream analysis. `import_inventory.csv` records direct imports; dynamically loaded files are identified in the source code. Figure rendering has its own tested requirements.

## Analysis groups

- `structural_variant_catalogue/`: paired calls, breakends and donor-signature definitions.
- `donor_level_duplication_analysis/`: donor proportions, bootstrap intervals and permutation comparisons.
- `spatial_null_library/`: event counting, coordinate clustering and genomic matching.
- `structural_variant_class_composition/`: class enrichment, local composition and the four-class virus contrast.
- `gene_length_matched_null/`: gene annotation, candidate pools and length-matched sensitivity.
- `replication_timing_sensitivity/`: common-support comparisons across timing references.
- `supporting_analyses/`: null calibration, external catalogues and molecular-context comparisons.
- `supplemental_context/`: clonality, liver-cohort and spatial-sensitivity analyses.
- `independent_verification/`: separate recounting and numerical checks.

Source paths and import names have been adapted for this release. Numerical calculations, statistical parameters and random seeds are unchanged. The verification record distinguishes source-code inspection from analyses executed with the packaged inputs.
