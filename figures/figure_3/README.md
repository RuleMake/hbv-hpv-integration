# Figure 3: Integration sites at duplication boundaries

The panels relate boundary and interior enrichment to site composition and exclusion-distance sensitivity. The preparation scripts copy the saved enrichment tables and calculate the display conversions used in the figure.

Run from this directory:

```sh
python prepare_sources.py
python prepare_display.py
python build_fig3.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 5 display tables:

- `plotted_data/Figure_3A_plotted_values.csv` — copied from a bundled result table.
- `plotted_data/Figure_3B_plotted_values.csv` — copied from a bundled result table.
- `plotted_data/Figure_3C_plotted_values.csv`
- `plotted_data/Figure_3E_display_values.csv`
- `plotted_data/Figure_3_composition_plotted_values.csv`

Copied tables retain the estimates from the underlying analysis.

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../README.md).
