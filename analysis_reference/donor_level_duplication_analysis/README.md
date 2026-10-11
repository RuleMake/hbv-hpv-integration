# Donor-level duplication analysis

These sources compare each donor's fraction of large duplications near a shared coordinate catalogue between detected-integration groups. Donor-bootstrap confidence intervals describe donor sampling, whereas the other spatial analyses use randomized genomic positions.

## Source files

| File | Purpose |
|---|---|
| `fig5c_donor_duplication_comparison.py` | Estimates donor-level near-duplication proportions, bootstrap intervals and permutation tests. |
| `test_fig5c_donor_duplication_comparison.py` | Checks donor proportions, event proximity and statistical edge cases. |

## Required inputs and layout

Duplication ledger; eligible viral coordinates; matching audit; 50-kb features; GISTIC peaks; histology table; input manifest and saved reference summaries.

The analysis directory requires `data/`, `provenance/` and `outputs/` subdirectories. Run `--stage reproduce` to check the reference summaries, then `--stage statistics` to compute the donor-level comparisons. Its tests also expect saved output tables; they are not part of automatic package verification.

See the release input inventory for data sources and `../source_map.csv` for source locations. Configure the required files, dependencies and output paths in a separate work directory before execution. Release verification checks syntax and file hashes without running these upstream analyses.
