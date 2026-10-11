# Supplemental Figure S2: Calibration, sampling pools and fold variation

The panels summarize null-model calibration, the availability of eligible sampling pools and variation in fold estimates. The builder reads the supplied display tables directly.

Run from this directory:

```sh
python build_figS2.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The builder reads these precomputed tables without rewriting them:

- `plotted_data/S2ABC_calibration.csv`
- `plotted_data/S2D_rejection_intervals.csv`
- `plotted_data/S2E_pool_diagnostics.csv`
- `plotted_data/S2F_fold_summaries.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../../README.md).
