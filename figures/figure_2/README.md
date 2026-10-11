# Figure 2: Null-model calibration and counting units

The panels compare enrichment estimates, calibration rates and catalogue sizes across counting units and sampling schemes.

Run from this directory:

```sh
python prepare_sources.py
python build_fig2.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 4 display tables:

- `plotted_data/Figure_2A_plotted_values.csv`
- `plotted_data/Figure_2B_plotted_values.csv`
- `plotted_data/Figure_2C_plotted_values.csv`
- `plotted_data/Figure_2D_plotted_values.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../README.md).
