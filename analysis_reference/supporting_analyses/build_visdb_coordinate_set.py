"""Prepare VISDB HBV and HPV coordinate sets for external comparison.

Coordinates are harmonized to the analysis genome build. Exact matches to the
primary Tian coordinate catalogue are removed, and overlap is summarized by
source study. These external sets remain separate from the primary estimate.
Subsequent analyses assess study composition, leave-one-study-out sensitivity
and genomic-covariate balance."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402

VISDB = Path(r"external_data\visdb")
SRC = Path(r"upstream_project")
CHAIN = SRC / "data" / "raw" / "annotation" / "hg38ToHg19.over.chain.gz"
OUT = Path(r"source_project\results")

SOURCES = {"VISDB_HBV": "Hepatitis-B_virus.csv", "VISDB_HPV": "Human_Papillomavirus.csv"}
TIAN_OF = {"VISDB_HBV": "HBV", "VISDB_HPV": "HPV"}
MATCH_TOL = 10          # bp; an exact re-publication of the same breakpoint, not a coincidence

HG38_LEN = {f"chr{c}": L for c, L in zip(
    list(range(1, 23)) + ["X"],
    [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
     138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
     83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895])}


def tian_index() -> dict[str, dict[str, np.ndarray]]:
    co = F.load_coordinates()
    out = {}
    for v in ("HBV", "HPV"):
        s = co[co["virus_class"] == v]
        out[v] = {str(c): np.sort(g["pos0"].to_numpy(np.int64))
                  for c, g in s.groupby("chrom")}
    return out


def matches_tian(idx: dict[str, np.ndarray], chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    hit = np.zeros(len(pos), dtype=bool)
    for c in np.unique(chrom):
        arr = idx.get(str(c))
        if arr is None or arr.size == 0:
            continue
        m = chrom == c
        p = pos[m]
        j = np.searchsorted(arr, p)
        best = np.full(len(p), np.iinfo(np.int64).max)
        for off in (-1, 0):
            k = np.clip(j + off, 0, len(arr) - 1)
            best = np.minimum(best, np.abs(arr[k] - p))
        hit[m] = best <= MATCH_TOL
    return hit


def build_one(name: str, fn: str, tidx, lo) -> tuple[pd.DataFrame, list]:
    audit = []

    def note(step, n, detail=""):
        audit.append(dict(arm=name, step=step, n_records=n, detail=detail))
        print(f"[{name}] {step:42s} n = {n:6,}   {detail}")

    d = pd.read_csv(VISDB / fn, low_memory=False)
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
    d = d[d["pos38"].between(1, d["chrom"].map(HG38_LEN))].copy()
    note("3. numeric + inside hg38 chrom length", len(d))

    rows, fail, multi, xchrom = [], 0, 0, 0
    for c38, p38, pmid in zip(d["chrom"], d["pos38"], d["PubMed ID"]):
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
        rows.append((c19, p19, pmid))
    out = pd.DataFrame(rows, columns=["chrom", "pos0", "pmid"])
    note("4. lifted hg38 -> hg19", len(out),
         f"unmapped {fail:,}, multi-mapping {multi:,}, chromosome-switch {xchrom:,}")

    before = len(out)
    out = out.drop_duplicates(subset=["chrom", "pos0"]).reset_index(drop=True)
    note("5. exact duplicate coordinates collapsed", len(out), f"{before - len(out):,} removed")

    hit = matches_tian(tidx[TIAN_OF[name]], out["chrom"].to_numpy(str),
                       out["pos0"].to_numpy(np.int64))
    per_study = (pd.DataFrame({"pmid": out["pmid"], "hit": hit})
                 .groupby("pmid")["hit"].agg(["sum", "size"]))
    per_study["rate"] = per_study["sum"] / per_study["size"]
    heavy = per_study[(per_study["rate"] > 0.5) & (per_study["size"] >= 20)]
    if len(heavy):
        print(f"[{name}] studies whose coordinates are >50% identical to the Tian catalogue "
              f"(these ARE Tian, however they are labelled):")
        print(heavy.sort_values("size", ascending=False).to_string())
    else:
        print(f"[{name}] no study is >50% identical to the Tian catalogue")
    out = out[~hit].reset_index(drop=True)
    note("6. coordinates matching Tian dropped (independence)", len(out),
         f"{int(hit.sum()):,} removed at +/-{MATCH_TOL} bp")

    out["arm"] = name
    print(f"[{name}] {out['pmid'].nunique()} contributing studies; "
          f"largest {out['pmid'].value_counts().iloc[0]:,} coordinates "
          f"({100*out['pmid'].value_counts().iloc[0]/len(out):.1f}%)")
    return out, audit


def main() -> None:
    from pyliftover import LiftOver
    lo = LiftOver(str(CHAIN))
    tidx = tian_index()
    frames, audits = [], []
    for name, fn in SOURCES.items():
        f, a = build_one(name, fn, tidx, lo)
        frames.append(f)
        audits.extend(a)
        print()
    coords = pd.concat(frames, ignore_index=True).sort_values(
        ["arm", "chrom", "pos0"]).reset_index(drop=True)
    coords.to_csv(OUT / "visdb_coordinates.csv", index=False)
    pd.DataFrame(audits).to_csv(OUT / "visdb_build_audit.csv", index=False)
    print(coords.groupby("arm").size().to_string())
    print(f"\nwrote {OUT / 'visdb_coordinates.csv'}")


if __name__ == "__main__":
    main()
