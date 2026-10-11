# Supplemental Figure S5: Hotspot membership and class definitions

The panels compare class-specific estimates, hotspot membership and alternative class definitions. The builder reads the supplied display table directly.

Run from this directory:

```sh
python build_figS5.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The builder reads these precomputed tables without rewriting them:

- `plotted_data/S5_plotted_values.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../../README.md).
