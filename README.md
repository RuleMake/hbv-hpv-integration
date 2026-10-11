# Host rearrangement landscapes distinguish HBV and HPV integration sites

Code and numerical data accompanying the manuscript by Yan Yan, Yuhan Zhou,
Haonan Liu, Tianli Liu and Xin Zheng. Version **2026.10.10**.

The package contains the data and scripts for Figures 1–6 and Supplemental
Figures S1–S11, Supplemental Tables S1–S7, and the study's statistical analysis
programs. A worked example reproduces the comparison of record, event and
patient weighting in the Porter integration catalogue.

## Reproduce the figures

With Python 3.12 installed, run these commands from this directory:

```bash
python -m pip install -r requirements-figures.txt
python verify_release.py
python verify_release.py --rebuild
```

The first verification command checks file integrity and numerical references.
The second prepares the plotting data and redraws all 17 figures in temporary
directories. It compares the results with the supplied reference figures and
saves the new figures and a report in `verification_output/`. All inputs needed
for these commands are included.

To select figures and choose where to save the results:

```bash
python verify_release.py --figures F3 S11 --output-dir my_check
```

Use `F1`–`F6` for main figures and `S1`–`S11` for supplemental figures.
`figures/FIGURES.json` lists their inputs and commands. Each figure directory
also contains a caption and instructions.

The reference figures use **Arial**. Install this font and use the package
versions in `requirements-figures.txt` for an exact visual match; differences
in fonts or rendering libraries can affect layout. Fonts are not distributed
with the code. All figures have SVG references, and all except S11 also have
PDF references.

## Reproduce the weighting analysis

```bash
python reproduce_porter/run.py --verify --output-dir verification_output/porter
```

This example starts from the included table of 129 published events in 69
patients with integration events. It performs 20,000 joint patient-bootstrap
replicates and checks the results against the supplied outputs. The patient
statistic averages within-patient proportions. See `reproduce_porter/README.md`
for the event definitions and weighting methods.

## Find data and analysis programs

| Directory or file | Contents |
|---|---|
| `figures/` | Plotting scripts, input data, captions and reference figures |
| `tables/` | Supplemental Tables S1–S7 and column definitions |
| `data/` | Numerical results and external-catalogue comparisons |
| `reproduce_porter/` | Weighting example, input table and expected results |
| `analysis_reference/` | Statistical analysis programs and their input requirements |
| `provenance/` | Data sources and links between reported values and analysis outputs |
| `release_manifest.json` | File sizes and SHA-256 checksums |
| `VERIFICATION.json` | Tested software environment and verification results |

The analysis programs cover spatial randomization, structural-variant counting,
replication timing, gene-length matching, transcription, clonality, microhomology
and external-cohort comparisons. Their README files describe the study data,
shared modules and directory structure they require. They are provided as
reference implementations for use with those inputs; installing the listed
dependencies alone does not reconstruct the full analysis environment.

Figure reproduction uses the numerical results included here. Raw sequencing
data, complete external callsets and genomic reference tracks must be obtained
from the originating studies. Available links and access details are listed in
`provenance/input_index.csv`.

## Checks performed

All 17 SVG figures and 16 PDF figures reproduced the reference images exactly
in the recorded environment. Before each run, the verifier removes the 34 CSV
outputs produced by the plotting commands; three of these are copies of saved
tables. A further 18 precomputed tables are checked without being regenerated.
These checks cover data preparation and plotting, rather than the upstream
statistical analyses.

The Porter example reproduced all eight numerical outputs byte for byte. The
other analysis programs were checked for syntax and source provenance. The
integrity check also verifies the seven supplemental tables, 91 mapped values
and 213 indexed table rows. The supplied manifest describes the released files;
reports created by subsequent runs are separate outputs.

## License and citation

The authors' software is distributed under the MIT license; see `LICENSE`.
Data from other studies remain subject to their original terms. Please cite
the accompanying manuscript and the relevant source studies when using these
materials. `CITATION.cff` contains the software citation metadata. The source repository is
[https://github.com/RuleMake/hbv-hpv-integration](https://github.com/RuleMake/hbv-hpv-integration). The permanent archive is available on Zenodo: [10.5281/zenodo.23291379](https://doi.org/10.5281/zenodo.23291379).
