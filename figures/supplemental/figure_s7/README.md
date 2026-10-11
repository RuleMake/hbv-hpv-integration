# Supplemental Figure S7: Composition contrasts across matching schemes

The panels compare a fixed set of integration loci across seven matching schemes and show the resulting changes in the composition contrast. These schemes share the same observations; their estimates are not independent replications. Single-covariate schemes are parallel comparisons, followed by joint matching schemes.

Run from this directory:

```sh
python prepare_sources.py
python build_figS7.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 2 display tables:

- `plotted_data/S7_conditioning_ladder_saved_values.csv` — copied from a bundled result table.
- `plotted_data/S7_plotted_values.csv`

Copied tables retain the estimates from the underlying analysis.

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../../README.md).
