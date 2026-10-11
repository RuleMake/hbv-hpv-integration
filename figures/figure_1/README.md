# Figure 1: Study design and integration-site distributions

The panels introduce the study design and show how integration records are grouped into loci. Catalogue counts, within-locus distances and genomic distributions are drawn from the bundled tables.

Run from this directory:

```sh
python prepare_sources.py
python build_fig1.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 4 display tables:

- `plotted_data/Figure_1B_legend_values.csv`
- `plotted_data/Figure_1C_plotted_values.csv`
- `plotted_data/Figure_1D_plotted_values.csv`
- `plotted_data/Figure_1D_query_exclusions.csv`

The following table is retained as supporting data; the builder does not read it:

- `plotted_data/Figure_1B_plotted_values.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../README.md).
