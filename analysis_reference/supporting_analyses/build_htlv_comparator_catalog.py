"""Prepare HTLV-1 integration coordinates from VISDB for spatial comparisons.

The source is Human_T-lymphotropic_virus.csv from VISDB (Tang et al. 2020,
Nucleic Acids Research 48:D633). Coordinates use the hg38_start field and are
converted to hg19 with pyliftover. The build retains numeric, in-range positions
on chromosomes 1-22 and X with unique, chromosome-consistent mappings.

Cell-line records are excluded. ATL, non-malignant and unrecorded clinical states
remain separate, and each filtering step is counted in the output summary."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

VISDB = Path(r"external_data\visdb\Human_T-lymphotropic_virus.csv")
SRC = Path(r"upstream_project")
CHAIN = SRC / "data" / "raw" / "annotation" / "hg38ToHg19.over.chain.gz"
OUT = Path(r"source_project\results")

HG38_LEN = {f"chr{c}": L for c, L in zip(
    list(range(1, 23)) + ["X"],
    [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
     138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
     83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895])}

TUMOUR_DISEASES = {"Adult T-cell leukemia", "Adult T-cell leukemia/lymphoma",
                   "Adult T-cell leukemia-lymphoma"}
NONMALIGNANT_DISEASES = {"Asymptomatic carrier",
                         "HTLV-1 associated myelopathy/tropical spastic paraparesis"}


def build() -> tuple[pd.DataFrame, pd.DataFrame]:
    audit = []

    def note(step, n, detail=""):
        audit.append(dict(step=step, n_records=n, detail=detail))
        print(f"[check]   {step:38s} n = {n:6,}   {detail}")

    d = pd.read_csv(VISDB, low_memory=False)
    note("0. VISDB records", len(d))

    d["sample_type"] = d["Sample type"].astype(str).str.strip()
    is_line = d["sample_type"].str.contains("cell line", case=False, na=False)
    d = d[~is_line].copy()
    note("1. cell lines dropped", len(d), f"{int(is_line.sum()):,} removed")

    d["chrom"] = d["Human chromosome"].astype(str).str.strip()
    d = d[d["chrom"].isin(HG38_LEN)].copy()
    note("2. chr1-22,X only", len(d))

    d["pos38"] = pd.to_numeric(d["hg38_start"], errors="coerce")
    d = d.dropna(subset=["pos38"])
    d["pos38"] = d["pos38"].astype(np.int64)
    inrange = d["pos38"].between(1, d["chrom"].map(HG38_LEN))
    d = d[inrange].copy()
    note("3. numeric + inside hg38 chrom length", len(d))

    # ---------------------------------------------------------------- arm assignment
    # CORRECTED 2026-08-19. The first version split on VISDB's `Sample type` field. That field
    # is annotated PER STUDY, not per sample, and it contradicts the disease field in both
    # directions: PMID 35654946 labels its asymptomatic-carrier samples "Tumor" (6,450 records),
    # while PMID 24735963 and 28811663 label ATL patients' samples "Nontumor" (15,508). A split
    # on it put carriers in the malignant arm and ATL patients in the non-malignant one, and it
    # was additionally collinear with study by construction -- the four "Nontumor" PMIDs and the
    # three "Tumor" PMIDs did not overlap at all, so no contrast between the two arms could be
    # read as biology rather than as batch.
    #
    # The split is now on DISEASE, which is the field that carries the clinical state, with an
    # explicit third group for records whose state is not recorded. Records with no disease
    # given are NOT silently pooled into either arm.
    #
    # Study confounding is reduced but not removed: PMID 35654946 is the ONLY study that
    # contributes to both the malignant and the non-malignant arm (5,045 ATL / 6,450 carrier),
    # so it is the only place a within-study selection contrast can be scored. That contrast is
    # the one to trust; the between-arm contrast is reported with the collinearity stated.
    dis = d["Disease name"].astype(str).str.strip()
    d["arm"] = np.where(
        dis.isin(TUMOUR_DISEASES), "HTLV_ATL",
        np.where(dis.isin(NONMALIGNANT_DISEASES), "HTLV_nonmalignant", "HTLV_unknown"))
    d["arm_sampletype"] = np.where(d["sample_type"].str.lower().eq("tumor"),
                                   "labelled_Tumor", "labelled_Nontumor")
    print()
    print(d.groupby(["arm", "arm_sampletype"]).size().to_string())
    print()
    print("[study]  records per PMID:")
    print(d.groupby(["arm", "PubMed ID"]).size().to_string())
    print()

    from pyliftover import LiftOver
    lo = LiftOver(str(CHAIN))
    rows, fail, multi, xchrom = [], 0, 0, 0
    for c38, p38, arm, armd, pmid, styp, dnam in zip(
            d["chrom"], d["pos38"], d["arm"], d["arm_sampletype"], d["PubMed ID"],
            d["sample_type"], dis):
        res = lo.convert_coordinate(c38, int(p38))
        if not res:
            fail += 1
            continue
        if len(res) > 1:
            multi += 1
            continue
        c19, p19 = res[0][0], int(res[0][1])
        if c19 != c38:
            xchrom += 1
            continue
        rows.append((arm, armd, c19, p19, pmid, styp, dnam))
    out = pd.DataFrame(rows, columns=["arm", "arm_sampletype", "chrom", "pos0", "pmid",
                                      "sample_type", "disease"])
    note("4. lifted hg38 -> hg19", len(out),
         f"unmapped {fail:,}, multi-mapping {multi:,}, chromosome-switch {xchrom:,}")

    # exact-duplicate coordinates within an arm are one observation, not several
    before = len(out)
    out = out.drop_duplicates(subset=["arm", "chrom", "pos0"]).reset_index(drop=True)
    note("5. exact duplicate coordinates collapsed", len(out), f"{before - len(out):,} removed")

    print()
    print(out.groupby("arm").size().to_string())
    return out.sort_values(["arm", "chrom", "pos0"]).reset_index(drop=True), pd.DataFrame(audit)


def main() -> None:
    coords, audit = build()
    OUT.mkdir(parents=True, exist_ok=True)
    coords.to_csv(OUT / "htlv_coordinates.csv", index=False)
    audit.to_csv(OUT / "htlv_build_audit.csv", index=False)
    print(f"\nwrote {OUT / 'htlv_coordinates.csv'}  ({len(coords):,} coordinates)")
    print(f"wrote {OUT / 'htlv_build_audit.csv'}")


if __name__ == "__main__":
    main()
