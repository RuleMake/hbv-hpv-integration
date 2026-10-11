# Supplemental Figure S6: Sensitivity to complex cases and amplicon masking

The panels compare exclusion and masking settings, case ratios and the extent of filtering. Display preparation retains the saved estimates and distinguishes spatial-null reference ranges from bootstrap confidence intervals.

Run from this directory:

```sh
python prepare_sources.py
python build_figS6.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 2 display tables:

- `plotted_data/S6_exclusion_accounting.csv`
- `plotted_data/S6_plotted_values.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../../README.md).
