# Supplemental Figure S11: HPV clonality and integration-site patterns

The panels display HPV clonality counts and the corresponding saved test summaries. The builder reads two local input tables and does not refit a biological model.

Run from this directory:

```sh
python build_figure.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The reference graphic for this figure is SVG only.

## Display tables

The commands write 1 display table:

- `displayed_presence_counts.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../../README.md).
