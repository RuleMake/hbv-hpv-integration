# Statistical verification sources

These sources document separate recounting and reconciliation procedures for the catalogue, composition, enrichment and sensitivity analyses. The required inputs and saved outputs are listed below.

## Source files

| File | Purpose |
|---|---|
| `verify_class_composition_axis.py` | Recomputes local composition summaries and density-matched results. |
| `verify_class_enrichment_atlas.py` | Recounts class-specific spatial statistics from saved input and null tables. |
| `verify_donor_proportion_statistics.py` | Re-derives donor proportions and their inferential summaries from saved inputs. |
| `verify_gene_length_matching.py` | Recounts gene-length-matched results and checks candidate-pool constraints. |
| `verify_replication_timing_sensitivity.py` | Recomputes timing-reference results and checks shared random streams. |
| `verify_structural_variant_catalogue.py` | Rebuilds source-call identities and reconciles saved catalogue representations. |
| `verify_virus_class_interaction.py` | Recomputes the four-class virus contrast and checks paired random streams. |

## Required inputs and layout

Analysis-specific source tables, outputs, candidate pools, null arrays and manifests. Exact dependencies differ by script and are visible in its path constants.

Six scripts execute data loading and write check summaries at module scope. Do not import or autodiscover this directory as a test suite. Their expected directory position relative to the source inputs and saved results must be reproduced in a separate work directory. Package verification parses their syntax only.

See the release input inventory for data sources and `../source_map.csv` for source locations. Configure the required files, dependencies and output paths in a separate work directory before execution. Release verification checks syntax and file hashes without running these upstream analyses.
