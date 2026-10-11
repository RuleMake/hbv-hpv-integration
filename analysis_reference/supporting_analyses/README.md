# Supporting analyses

These sources cover null calibration, external coordinate catalogues, the HTLV-1 comparator, junction exclusion, transcript-status enrichment and amplicon slope comparisons. They also include the shared helpers required by these analyses.

## Source files

| File | Purpose |
|---|---|
| `amplification_and_histology_control.py` | Implements amplification masking and histology-related analyses. |
| `boundary_geometry_junction_exclusion.py` | Repeats boundary geometry after excluding variants close to integration junctions. |
| `build_htlv_comparator_catalog.py` | Filters and harmonizes HTLV-1 comparator coordinates. |
| `build_visdb_coordinate_set.py` | Harmonizes VISDB coordinates and removes overlap with the primary coordinate catalogue. |
| `clca_class_replication.py` | Loads CLCA structural-variant and integration records and computes class comparisons. |
| `coamplification_slope_control_and_signflip.py` | Fits matched amplicon slope estimators and a patient sign-flip comparison. |
| `corrected_estimator_baseline.py` | Loads primary catalogues and implements the original and locus-based spatial estimators. |
| `covariate_balance_check.py` | Summarizes observed and randomized genomic-covariate balance. |
| `discriminating_tests.py` | Implements boundary geometry and associated spatial controls. |
| `external_catalogue_class_composition.py` | Computes class-composition statistics for the CLCA and Sasa coordinate sets. |
| `integrase_positive_boundary_geometry.py` | Computes the HTLV-1 duplication-boundary geometry comparator. |
| `integrase_positive_class_atlas.py` | Computes the HTLV-1 class-enrichment comparator. |
| `integrase_positive_composition_axis.py` | Computes the HTLV-1 class-composition comparator. |
| `integrase_positive_covariate_balance.py` | Checks genomic-covariate balance for comparator placement pools. |
| `joint_correction.py` | Combines covariate matching with amplification-peak masking. |
| `locus_class_propensity.py` | Defines class-composition counting and summary routines used by external comparisons. |
| `minimal_bigwig_reader.py` | Reads compressed bigWig blocks without a separate genome-track library. |
| `null_calibration_baseline.py` | Simulates the baseline counting and random-placement calibration comparisons. |
| `null_calibration_recordblock_and_rotation.py` | Evaluates counting and placement combinations under block and rotation generators. |
| `replication_catalogue_power_probe.py` | Computes design-specific resolution summaries for external coordinate sets. |
| `sasa_external_replication.py` | Loads and harmonizes the Sasa coordinate catalogue and computes spatial enrichment. |
| `score_visdb_class_composition.py` | Computes class-composition statistics for the VISDB coordinate sets. |
| `transcript_status_duplication_enrichment.py` | Compares duplication enrichment across detected fusion-transcript categories. |

## Required inputs and layout

Primary coordinate and variant catalogues; genomic features; external VISDB, CLCA, Sasa and Fan source tables; genome-build chain file; GISTIC peaks; amplicon measurements; stored reference outputs.

Several scripts use relative `source_project/data` and `source_project/results` paths, sometimes with Windows separators. Configure these paths for a separate analysis work directory and keep writable outputs apart from source data.

See the release input inventory for data sources and `../source_map.csv` for source locations. Configure the required files, dependencies and output paths in a separate work directory before execution. Release verification checks syntax and file hashes without running these upstream analyses.

## Amplicon slope intervals

`coamplification_slope_control_and_signflip.py` uses normal-based intervals. The current Figure 6 `prepare_sources.py` computes the t-based intervals displayed in the manuscript from saved estimates and standard errors. Use the figure-preparation calculation to reproduce the displayed intervals.
