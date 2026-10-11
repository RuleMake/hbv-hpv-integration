# Figure 6: Copy number and integration-support relationships

The panels examine chromosomal-arm, amplicon and sample-level copy-number relationships. The preparation script also reproduces the reported unweighted ordinary least-squares fits to sample means, with two-sided t intervals using n − 2 degrees of freedom. The earlier normal-approximation summaries remain separately labelled in the source tables.

Run from this directory:

```sh
python prepare_sources.py
python build_fig6.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 5 display tables:

- `plotted_data/Figure_6A_plotted_values.csv`
- `plotted_data/Figure_6B_plotted_values.csv`
- `plotted_data/Figure_6C_fits_plotted_values.csv`
- `plotted_data/Figure_6C_line_data.csv`
- `plotted_data/Figure_6C_plotted_values.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../README.md).
