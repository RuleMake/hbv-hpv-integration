# Matching schemes and detectable-effect markers

Supplemental Figure S7 compares the same integration loci under seven genomic matching schemes. The single-covariate comparisons match replication timing, gene density or common fragile-site status separately; the later schemes use joint matching. They are alternative analyses of shared observations, not independent replications.

The supplied `conditioning_ladder_saved.csv` records the composition score, its null expectation and their difference:

```
excess = observed score − mean null score
z = excess / null standard deviation
```

The per-virus excess was saved to five decimal places and z to two. Its detectable-effect field was left empty. For the display, `prepare_sources.py` reconstructs the approximate two-sided threshold as:

```
threshold in percentage points = 100 × 2.80 × |excess / z|
```

The factor 2.80 corresponds to the approximate threshold for 80% power at a two-sided significance level of 0.05. It is not a confidence limit. Filled markers use the saved rule `|z| ≥ 2.80`.

The script checks the effect and z values against the saved matching-scheme table. It also evaluates rounding bounds using ±0.000005 for excess and ±0.005 for z, and confirms that rounding cannot change any marker classification. Source values and empty detectable-effect fields remain intact; reconstructed values are placed in additional display columns.

In the reduction panel, the sum of the three single-covariate reductions and the increment from adding local density are arithmetic comparisons. They are distinct from estimates obtained under a separately evaluated matching scheme. These changes describe how the contrast depends on matching choices; they do not establish whether a covariate is a confounder or mediator.
