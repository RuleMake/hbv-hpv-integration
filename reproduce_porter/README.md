# Reproduce the Porter event-weighting example

This example compares how breakpoint records, author-defined integration events, and equal patient contributions represent the same published integration catalogue. The comparison uses a single catalogue from Porter et al. (2025), with all three estimators applied to the same events.

The supplied analytical table contains **129 author-defined events in 69 patients with integration events, represented by 438 breakpoint records**. The published analysis cohort contained 72 patients; the three patients with no eligible event do not enter these event-composition estimators.

## Run

Use Python 3.12 with NumPy and pandas. From the repository root:

```sh
python -m pip install -r reproduce_porter/requirements.txt
python reproduce_porter/run.py --verify
```

The second command reads the bundled derived table, performs the fixed 20,000-replicate joint patient bootstrap, writes results to `reproduce_porter/outputs/`, and compares all eight numerical output files with `expected/`. It requires no network access, credentials, original article files, or external datasets.

To choose a different output directory:

```sh
python reproduce_porter/run.py --verify --output-dir reproduce_porter/outputs_rerun
```

The input and reference directories are protected against accidental overwriting. The release was tested with Python 3.12.10, NumPy 2.3.5, and pandas 2.3.3. The original saved analysis used NumPy 2.3.5 and pandas 3.0.1; `--verify` checks the numerical results directly rather than relying on matching version labels.

## Estimands and resampling

For event *e* in patient *p*, let *m_e* be the number of published breakpoint records, *n_p* the patient's number of eligible events, *N* the number of eligible events, and *P* the number of patients with an eligible event. Events receive one of three weights:

- **Breakpoint-record weighting:** `m_e / sum(m_e)`.
- **Author-event weighting:** `1 / N`.
- **Equal-patient weighting:** `1 / (P * n_p)`.

The third estimator averages within-patient event proportions, which differs from the proportion of patients positive for a category. Record counts represent published breakpoint records, not sequencing molecules or read support. The event definitions and structural labels are those of the source authors; an event need not correspond to one historical integration episode.

The structural categories are mutually exclusive author labels: multi-breakpoint, dup-like, del-like, unmatched, translocation, and repeat. Published fusion-transcript detection is a separate yes/no measure. A published “no” does not establish biological silence or verify negative-assay sensitivity. Reported ecDNA annotations in the input are retained for provenance and are not treated as a seventh exclusive structural category.

The random-number generator uses seed **20260909**. Each of **20,000** draws samples 69 patients with replacement, keeping all events of each selected patient together. The same multinomial patient multiplicities are used for all three estimators and their paired contrasts. Confidence intervals are pointwise 2.5th and 97.5th percentiles using NumPy's linear quantile method. Leave-one-patient-out contrasts describe sensitivity to individual patients. Total variation summarizes redistribution across the six structural labels; no primary null-hypothesis tests are performed.

## Files

`inputs/normalized_author_events.csv` is the fixed derived input table. The estimator consumes `patient`, `event_id`, `category`, `record_count`, and `fusion_transcript_detected`; the remaining columns retain published labels, breakpoint identifiers, source-row references, and eligibility annotations. Patient identifiers are the published catalogue labels. The input SHA-256 is:

```text
1db96456448f3d789ba703013c3df9368263c63632bdfe8efa9b2b3e9b05433e
```

`run_config.json` fixes the input checksum, categories, denominators, seed, and bootstrap size. The script verifies the checksum and the 129/69/438 denominators before estimation.

| Output | Contents |
| --- | --- |
| `weight_estimates.csv` | Three weighted estimates and their patient-bootstrap intervals for each measure |
| `weight_differences.csv` | Paired contrasts, intervals, and leave-one-patient-out extrema |
| `composition_redistribution.csv` | Descriptive total variation across the six structural labels |
| `event_weights.csv` | Derived input rows with the three event weights |
| `patient_contributions.csv` | Patient event/record totals and category contributions |
| `leave_one_patient_out.csv` | Contrasts after omitting each patient in turn |
| `patient_bootstrap_counts.npy` | The 20,000 × 69 joint patient multiplicity matrix |
| `joint_bootstrap_estimates.npz` | Replicate estimates, method names, measure names, and patient order |
| `run_metadata.json` | Current input/configuration/script hashes, versions, and denominators |
| `verification.json` | Reference comparison, produced when `--verify` is specified |

The eight numerical reference files are included in `expected/`, with checksums in `expected/manifest.json`. The verifier checks row/column and array ordering, exact integer draws and labels, and floating-point values to relative and absolute tolerances of `1e-12`. It also reports byte identity and exact value identity. Runtime metadata is excluded because the portable script, configuration path, and pandas version differ from the original run.

For a quick point-estimate check, the multi-breakpoint percentages are approximately **60.05%, 31.78%, and 27.39%**, and the dup-like percentages are **16.89%, 28.68%, and 38.06%**, in record/event/patient order. These are alternative summaries of one common event set, not measurements from three independent samples.

## Provenance

Porter VL et al. (2025). *Rearrangements of viral and human genomes at human papillomavirus integration events and their allele-specific impacts on cancer genome regulation.* Genome Research 35:653–670. [doi:10.1101/gr.279041.124](https://doi.org/10.1101/gr.279041.124).

The derived input was checked against the authors' published supplementary event and breakpoint tables, available in the [Europe PMC supplementary package](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12047271/supplementaryFiles). The normalized table retains the source row indices. The original article text and supplementary workbooks are not included here. This example starts from the derived event table; source-workbook parsing is outside its scope.

`PROVENANCE.json` records the input and source-script checksums. Aggregation, weighting, patient resampling, quantiles and contrasts retain the original algorithm; the command-line interface adds configurable output paths and comparison with the supplied reference results.
