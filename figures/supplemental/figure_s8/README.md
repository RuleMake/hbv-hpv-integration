# Supplemental Figure S8: HTLV-1 catalogue comparisons

The panels show HTLV-1 class estimates and sensitivity analyses. The builder reads the supplied display table and checks the hashes of its bundled sources.

Run from this directory:

```sh
python build_figS8.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The builder reads these precomputed tables without rewriting them:

- `plotted_data/S8_source_data.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../../README.md).
