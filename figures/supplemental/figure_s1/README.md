# Supplemental Figure S1: Catalogue composition and counting units

The panels trace catalogue filtering and identify the input sets used in each analysis. Counts retain their original units: paired structural-variant calls, breakends, intervals or large-duplication records.

Run from this directory:

```sh
python prepare_sources.py
python build_figS1.py
```

`INPUTS.csv` records the bundled inputs and their hashes. The caption is in `legend.md`; the manuscript reference graphics are in `reference/`. New graphics are exported separately.

The exporter also writes PNG previews and layout-check reports.

## Display tables

The commands write 3 display tables:

- `plotted_data/S1A_catalogue_flow.csv`
- `plotted_data/S1B_analysis_inputs.csv`
- `plotted_data/S1C_reference_units.csv`

These commands reproduce the figure from saved results. Requirements and the scope of verification are described in the [figure reproduction guide](../../README.md).
