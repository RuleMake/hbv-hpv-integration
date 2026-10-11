# Reproducing the figures

This directory contains the six main figures and eleven supplemental figures. Each figure directory includes a plotting script, saved input tables, a caption and reference graphics. `FIGURES.json` lists the commands and comparison files; each `INPUTS.csv` records the source hashes.

Run the commands in a figure's README from that directory, or use the repository verifier. Figure 3 requires two preparation steps before plotting. Figures S8–S10 check their bundled inputs directly. Figure S11 reads two local tables and has an SVG reference only; the other sixteen figures have both SVG and PDF references.

The commands reproduce the graphics from saved results and apply the display transformations documented in the scripts. Figure 6 also reproduces the reported ordinary least-squares fits to sample means. Raw-data processing and the upstream permutation or bootstrap analyses are outside this figure-reproduction workflow. Reference graphics in `reference/` are kept separate from new exports.

Arial is required for matching text layout. Verification compares numerical tables and rendered graphics: file bytes alone are unsuitable for vector exports because timestamps and element identifiers can vary.

## Display tables

The figure directories contain 52 display CSV files. The documented commands write 34 of them, including three copies of saved result tables: the Figure 3A and 3B tables and the Supplemental Figure S7 matching-scheme table. The remaining 18 are precomputed inputs. The Figure 1B table is retained as supporting data and is not read by its builder; the other precomputed tables are read directly.

`FIGURES.json` records these roles under `generated_tables`, `input_tables` and `copied_saved_tables`. The complete list is under `tables`. Recreating a display table does not imply that its underlying statistical analysis has been repeated.

For verification, all 34 generated tables are removed from an isolated copy before the commands run. The verifier checks that they are recreated and that the 18 retained tables remain unchanged, then compares the numerical values and rendered graphics with the references.
