# Structural-variant class composition

These sources estimate class-specific enrichment and the local composition contrast between HBV and HPV integration loci. `atlas_core.py` contains the event-index and summary routines used by the atlas and gene-length analyses.

## Source files

| File | Purpose |
|---|---|
| `atlas_core.py` | Provides event indexes, tracked random draws and enrichment summaries. |
| `local_sv_class_composition.py` | Computes the local class-composition axis and the HBV-HPV comparison. |
| `sv_class_enrichment_atlas.py` | Computes class-specific spatial enrichment at integration loci. |
| `virus_class_enrichment_contrast.py` | Computes the four-class contrast from shared spatial randomizations. |

## Required inputs and layout

Normalized paired calls and donor signatures; integration coordinates; genomic matching context; saved comparison summaries; saved spatial-null arrays for contrasts.

The expected `atlas/` and `axis/` directories each contain `outputs/`, with a sibling `shared/` directory. Recover this layout or explicitly reconfigure all input and output constants. The source files rely on these configured input and output relationships.

See the release input inventory for data sources and `../source_map.csv` for source locations. Configure the required files, dependencies and output paths in a separate work directory before execution. Release verification checks syntax and file hashes without running these upstream analyses.
