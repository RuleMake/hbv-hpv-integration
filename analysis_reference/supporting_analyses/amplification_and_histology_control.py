"""Evaluate duplication enrichment after amplification masking and histology restriction.

Recurrent-amplification peaks come from the PCAWG pan-cancer GISTIC2 catalogue.
Peaks are padded by the analysis window and excluded from both observed coordinates
and null candidate pools. A second comparison restricts the variant background to
Liver-HCC for HBV and Cervix-SCC, Cervix-AdenoCA and Head-SCC for HPV.

The matching strata, event-counting rule, leader clustering and seed are shared
with the primary estimator. EBV has no corresponding histology in the variant
catalogue, so its comparison evaluates peak masking against the full background."""

from __future__ import annotations

import argparse
import gzip
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402

SRC = Path(r"upstream_project")
GISTIC = SRC / "data" / "raw" / "pcawg_cnv" / "all_lesions.conf_95.rmcnv.pt_170207.txt.gz"
HISTOLOGY = SRC / "data" / "raw" / "clinical" / "pcawg_specimen_histology_August2016_v9.xlsx"
OUT = Path(r"source_project\results")
WINDOW = 100_000
BIN = F.BIN_SIZE

HISTO = {
    "HBV": ["Liver-HCC"],
    "HPV": ["Cervix-SCC", "Cervix-AdenoCA", "Head-SCC"],
}


def amp_peaks() -> pd.DataFrame:
    """GISTIC2 pan-cancer amplification peaks, hg19, from Wide Peak Limits."""
    rows = []
    with gzip.open(GISTIC, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        wide = header.index("Wide Peak Limits")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if not p or not p[0].startswith("Amplification Peak"):
                continue
            m = re.match(r"chr([0-9XY]+):(\d+)-(\d+)", p[wide].strip())
            if m:
                rows.append((f"chr{m.group(1)}", int(m.group(2)), int(m.group(3))))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end"]).drop_duplicates()
    print(f"[gistic]  {len(df)} pan-cancer amplification peaks "
          f"(median width {int((df.end - df.start).median()):,} bp)")
    return df


def masked_bins(peaks: pd.DataFrame, pad: int) -> set:
    out = set()
    for c, s, e in peaks.itertuples(index=False):
        lo = max(0, (s - pad) // BIN) * BIN
        hi = ((e + pad) // BIN) * BIN
        for b in range(lo, hi + BIN, BIN):
            out.add((c, b))
    return out


class MaskedNuller(F.Nuller):
    """Published 18-stratum pools with amplification-peak bins removed from the candidate set."""

    def __init__(self, mask: set) -> None:
        super().__init__()
        kept = {}
        removed = total = 0
        for (c, s), arr in self.bins.items():
            keep = np.array([(c, int(b)) not in mask for b in arr], dtype=bool)
            total += len(arr)
            removed += int((~keep).sum())
            if keep.any():
                kept[(c, s)] = arr[keep]
        self.bins = kept
        print(f"[mask]    {removed:,} / {total:,} candidate bins removed "
              f"({100 * removed / total:.1f}%); {len(self.bins)} non-empty pools remain")


MIN_POOL = 5


def viable(frame: pd.DataFrame, nuller: F.Nuller) -> np.ndarray:
    """Coordinates whose matched pool still holds >= MIN_POOL candidates under this nuller.

    Coordinates that fail are DROPPED and counted, never given a fallback pool: falling back to
    the unmasked pool would reintroduce exactly the amplification-peak bins the arm exists to
    remove, and would silently turn the restricted arm back into the baseline.
    """
    ok = np.zeros(len(frame), dtype=bool)
    for i, (c, s, off, p0) in enumerate(zip(frame["chrom"], frame["stratum"],
                                            frame["within_bin_offset"], frame["pos0"])):
        bins = nuller.bins.get((str(c), int(s)))
        if bins is None or bins.size == 0:
            continue
        cand = np.unique(bins + int(off))
        ok[i] = (cand != int(p0)).sum() >= MIN_POOL
    return ok


def restrict_events(histos: list[str] | None):
    """Rebuild the event index, optionally keeping only donors of the given histologies."""
    ev = pd.read_csv(F.BASE / "input_pooled_canonical_events.tsv", sep="\t")
    ev = ev.drop_duplicates("unique_td_id")
    if histos is not None:
        h = pd.read_excel(HISTOLOGY)[["icgc_donor_id", "histology_abbreviation"]]
        h = h.drop_duplicates("icgc_donor_id").set_index("icgc_donor_id")
        keep = ev["donor"].map(h["histology_abbreviation"]).isin(histos)
        ev = ev[keep]
    codes = pd.factorize(ev["unique_td_id"])[0].astype(np.int64)
    ev = ev.assign(_eid=codes)
    idx = {
        str(c): F.ChromEvents(g["left_pos0"].to_numpy(np.int64),
                              g["right_pos0"].to_numpy(np.int64),
                              g["_eid"].to_numpy(np.int64))
        for c, g in ev.groupby("chrom")
    }
    return idx, len(ev), ev["donor"].nunique()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()

    peaks = amp_peaks()
    mask = masked_bins(peaks, WINDOW)
    co = F.load_coordinates()
    nul_full = F.Nuller()
    nul_mask = MaskedNuller(mask)

    idx_full, n_full, d_full = restrict_events(None)
    print(f"[events]  full pool {n_full:,} events / {d_full:,} donors")
    idx_h = {}
    for v, hs in HISTO.items():
        idx_h[v], n, d = restrict_events(hs)
        print(f"[events]  {v} histology-restricted pool {n:,} events / {d:,} donors ({'+'.join(hs)})")

    rows = []
    for vc in ("HBV", "HPV", "EBV"):
        sub_all = co[co["virus_class"] == vc].sort_values(
            ["chrom", "pos0"]).reset_index(drop=True)
        in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                            for c, p in zip(sub_all["chrom"], sub_all["pos0"])])
        sub_keep = sub_all[~in_mask].reset_index(drop=True)
        print(f"\n[{vc}] {len(sub_all)} coords; {int(in_mask.sum())} "
              f"({100 * in_mask.mean():.1f}%) fall in a padded amplification peak")

        arms = [
            ("0_baseline", sub_all, nul_full, idx_full),
            ("1_ampPeakExcluded", sub_keep, nul_mask, idx_full),
            ("2_histologyRestricted", sub_all, nul_full, idx_h.get(vc, idx_full)),
            ("3_both", sub_keep, nul_mask, idx_h.get(vc, idx_full)),
        ]
        for tag, frame, nuller, index in arms:
            if len(frame) == 0:
                continue
            v = viable(frame, nuller)
            if not v.all():
                print(f"        [{tag}] dropping {int((~v).sum())} coordinates whose masked "
                      f"pool holds <{MIN_POOL} candidates", flush=True)
                frame = frame[v].reset_index(drop=True)
            if len(frame) == 0:
                continue
            labels = F.leader_clusters(frame["chrom"].to_numpy(str),
                                       frame["pos0"].to_numpy(np.int64), WINDOW)
            n_clu = int(len(np.unique(labels)))
            t0 = time.time()
            rng = np.random.default_rng(a.seed)
            try:
                o, nl, _, s_obs, s_null = F.est_block(index, nuller, rng, frame, WINDOW,
                                                      a.n_perm, labels, unit="uniqSV")
            except RuntimeError as exc:
                print(f"  [{vc} {tag}] FAILED: {exc}", flush=True)
                continue
            base = {"virus": vc, "arm": tag, "n_coords": len(frame), "n_clusters": n_clu,
                    "histology_pool": vc in idx_h and tag in ("2_histologyRestricted", "3_both")}
            rows.append({**base, "estimator": "bc_uniqSV_blockNull", **F.summarize(o, nl)})
            rows.append({**base, "estimator": "f_clusterSignTest", **F.summarize(s_obs, s_null)})
            print(f"  [{vc} {tag}] n={len(frame)} clusters={n_clu} "
                  f"fold={rows[-2]['fold']} p={rows[-2]['emp_p']:.5g} | "
                  f"sign {rows[-1]['fold']} p={rows[-1]['emp_p']:.5g} "
                  f"({time.time() - t0:.0f}s)", flush=True)

    out = pd.DataFrame(rows)
    dest = OUT / "fig3_amplification_histology_control.csv"
    out.to_csv(dest, index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print("\n" + out[["virus", "arm", "estimator", "n_coords", "n_clusters",
                      "observed", "expected", "fold", "fold_ci_lo", "fold_ci_hi",
                      "emp_p"]].to_string(index=False))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
