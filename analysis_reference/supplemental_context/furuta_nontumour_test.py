"""Prepares Furuta non-tumour coordinates and computes duplication enrichment.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402
import joint_correction as J  # noqa: E402
import amplification_and_histology_control as A  # noqa: E402

DATA = Path(r"source_project\data\furuta2018"
            r"\oncotarget-09-25075-s002.xlsx")
OUT = J.OUT
BIN = F.BIN_SIZE
WINDOW = 100_000
SEED = 7
HOTSPOTS = {"5": (1_000_000, 2_000_000),        # TERT
            "19": (36_000_000, 37_000_000),     # KMT2B
            "2": (216_000_000, 217_000_000)}    # FN1


def load_furuta() -> pd.DataFrame:
    d = pd.read_excel(DATA, sheet_name="ST5", header=None, skiprows=3)
    d = d.iloc[:, :11]
    d.columns = ["sample", "_a", "c1", "p1", "s1", "_b", "c2", "p2", "s2", "_c", "insertion"]
    d = d[d["sample"].notna()].copy()
    d["tissue"] = d["sample"].astype(str).str.extract(r"_(NT|T)$")[0]
    d["donor"] = d["sample"].astype(str).str.replace(r"_(NT|T)$", "", regex=True)

    # the human side is whichever junction is not the HBV reference (AP011098)
    j1_hum = ~d["c1"].astype(str).str.startswith("AP")
    d["h_chr"] = np.where(j1_hum, d["c1"].astype(str), d["c2"].astype(str))
    d["h_pos"] = np.where(j1_hum, d["p1"], d["p2"])
    d["h_strand"] = np.where(j1_hum, d["s1"], d["s2"])
    d["v_pos"] = np.where(j1_hum, d["p2"], d["p1"])
    d["v_strand"] = np.where(j1_hum, d["s2"], d["s1"])

    d = d[~d["h_chr"].str.startswith("AP")].copy()
    d["h_pos"] = pd.to_numeric(d["h_pos"], errors="coerce")
    d = d[d["h_pos"].notna() & d["tissue"].notna()]
    d["chrom"] = "chr" + d["h_chr"].str.replace("chr", "", regex=False)
    d["pos0"] = d["h_pos"].astype(np.int64)
    d = d[d["chrom"] != "chrY"]
    print(f"[furuta] {len(d)} human-side junctions | "
          f"{d.groupby('tissue').size().to_dict()} | {d['donor'].nunique()} donors")
    return d.reset_index(drop=True)


def hotspot_frac(d):
    hit = np.zeros(len(d), dtype=bool)
    for c, (lo, hi) in HOTSPOTS.items():
        hit |= (d["chrom"] == f"chr{c}") & d["pos0"].between(lo, hi)
    return float(hit.mean()), int(hit.sum())


def prep(d, nul, mask):
    in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                        for c, p in zip(d["chrom"], d["pos0"])])
    s = d[~in_mask].copy()
    k = pd.MultiIndex.from_arrays([s["chrom"].astype(str).to_numpy(),
                                   (s["pos0"].to_numpy(np.int64) // BIN) * BIN])
    s["stratum"] = nul.lookup.reindex(k).to_numpy()
    s = s[s["stratum"].notna()].copy()
    s["stratum"] = s["stratum"].astype(int)
    s["within_bin_offset"] = s["pos0"] - (s["pos0"] // BIN) * BIN
    return s.sort_values(["chrom", "pos0"]).reset_index(drop=True)


def run(idx, nul, sub, tag, n_perm):
    labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                               sub["pos0"].to_numpy(np.int64), WINDOW)
    n_clu = int(len(np.unique(labels)))
    nul.n_calls = nul.n_fallback = 0
    rng = np.random.default_rng(SEED)
    o, nl, _, s_obs, s_null = F.est_block(idx, nul, rng, sub, WINDOW, n_perm,
                                          labels, unit="uniqSV")
    fb = round(100 * nul.n_fallback / max(nul.n_calls, 1), 2)
    out = []
    for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                     ("f_clusterSignTest", F.summarize(s_obs, s_null))):
        out.append(dict(subset=tag, estimator=est, n_coords=len(sub), n_clusters=n_clu,
                        pct_pool_fallback=fb, **res))
    print(f"  [{tag:34s}] n={len(sub):5d} clu={n_clu:4d}  bc={out[0]['fold']:.4f} "
          f"p={out[0]['emp_p']:.5g} | sign={out[1]['fold']:.4f} p={out[1]['emp_p']:.5g}",
          flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    a = ap.parse_args()

    d = load_furuta()
    for t, g in d.groupby("tissue"):
        fr, n = hotspot_frac(g)
        print(f"[furuta] {t:3s} hotspot fraction (TERT+KMT2B+FN1 1-Mb bins): "
              f"{100*fr:.2f}% ({n}/{len(g)})")

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    idx, _ = F.load_events()
    nul = J.JointNuller(feat, F.Nuller(), mask)

    rows, t0 = [], time.time()
    print("\n[joint-corrected estimator, +/-100 kb, vs the 31,201-event PCAWG pool]")
    for tag, sel in (("NON-TUMOUR (unselected)", d.tissue == "NT"),
                     ("TUMOUR (same donors)", d.tissue == "T"),
                     ("ALL Furuta junctions", d.tissue.notna())):
        sub = prep(d[sel], nul, mask)
        if len(sub) < 30:
            continue
        rows += run(idx, nul, sub, tag, a.n_perm)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "furuta_nontumour_test.csv", index=False)
    pd.set_option("display.width", 240)
    print("\n" + df[["subset", "estimator", "n_coords", "n_clusters", "observed",
                     "expected", "fold", "fold_ci_lo", "fold_ci_hi", "emp_p"]]
          .to_string(index=False))
    print("\nreference under the SAME joint correction:")
    print("  Tian HBV (tumour)        bc 1.0186 [0.989, 1.049] p=0.110 | sign 1.0071 p=0.386")
    print("  Tian HBV non-hotspot 93% bc 1.0126 [0.983, 1.043] p=0.205 | sign 1.0002 p=0.501")
    print("  Tian HPV (survives)      bc 1.0836 p<5e-5          | sign 1.0791 p=3.5e-4")
    print(f"\nwrote {OUT / 'furuta_nontumour_test.csv'}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
