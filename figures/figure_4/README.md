# Figure 4: Rearrangement classes and composition contrasts

The panels connect class-specific enrichment with changes in rearrangement composition and the contrast between viral catalogues. Display preparation retains the reported estimators, reference ranges and multiplicity corrections.

Run from this directory:

```sh
python prepare_sources.py
python build_fig4.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 4 display tables:

- `plotted_data/Figure_4A_plotted_values.csv`
- `plotted_data/Figure_4B_detail_plotted_values.csv`
- `plotted_data/Figure_4B_total_plotted_values.csv`
- `plotted_data/Figure_4C_plotted_values.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../README.md).
