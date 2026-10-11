# Supplemental Figure S10: Tissue-exchange and external-catalogue comparisons

The panels compare tissue-exchange results with external integration catalogues. The builder assembles the displayed estimates from three supplied tables and checks their source hashes.

Run from this directory:

```sh
python build_figS10.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 1 display table:

- `plotted_data/Figure_S10_plotted_values.csv`

The builder also reads these precomputed tables without rewriting them:

- `plotted_data/External_composition.csv`
- `plotted_data/Tissue_exchange.csv`
- `plotted_data/VISDB_composition.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../../README.md).
