# Spatial-null library

These shared modules define the event counting and genomic background used by the class atlas and sensitivity analyses. Matching rules and amplification masks are kept in separate modules so their roles can be inspected.

## Source files

| File | Purpose |
|---|---|
| `amplification_peak_mask.py` | Builds the recurrent-amplification exclusion mask. |
| `bigwig_reader.py` | Reads bigWig tracks for genomic-context annotation. |
| `joint_covariate_matching.py` | Constructs matching strata using replication timing, local timing gradient and mappability. |
| `shared_analysis_context.py` | Loads consistent coordinate, variant and genomic-context tables for the spatial analyses. |
| `spatial_null_core.py` | Implements event counting, leader clustering and spatial placement routines. |

## Required inputs and layout

Coordinate catalogues; normalized event tables; 50-kb genomic features; chromosome sizes; recurrent-amplification peaks; mappability and replication-timing annotations.

These modules read sibling `data/` directories and structural-variant catalogue outputs. Import statements use the release filenames; relative input paths must be configured consistently in the analysis work directory.

See the release input inventory for data sources and `../source_map.csv` for source locations. Configure the required files, dependencies and output paths in a separate work directory before execution. Release verification checks syntax and file hashes without running these upstream analyses.
