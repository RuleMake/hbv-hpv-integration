# Supplemental Figure S4: Window sensitivity and cluster influence

The panels compare window sizes and summarize the influence of individual clusters. The builder reads the supplied display tables directly.

Run from this directory:

```sh
python build_figS4.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The builder reads these precomputed tables without rewriting them:

- `plotted_data/S4A_display.csv`
- `plotted_data/S4B_display.csv`
- `plotted_data/S4C_display.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../../README.md).
