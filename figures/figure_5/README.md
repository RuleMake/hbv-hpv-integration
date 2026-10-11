# Figure 5: Genomic matching and donor-group comparisons

The panels compare the composition contrast across genomic matching schemes and donor groups. Display labels and units are prepared from the saved estimates.

Run from this directory:

```sh
python prepare_sources.py
python build_fig5.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 3 display tables:

- `plotted_data/Figure_5A_plotted_values.csv`
- `plotted_data/Figure_5B_plotted_values.csv`
- `plotted_data/Figure_5C_plotted_values.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../README.md).
