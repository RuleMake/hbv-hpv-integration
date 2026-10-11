"""Compare HTLV-1 duplication-boundary geometry with the HBV and HPV references.

All coordinate sets use the same interval geometry, estimator and junction-exclusion
radii of 0, 1, 10 and 50 kb. HBV and HPV are recomputed and checked against saved
reference results. Minimum detectable displacements accompany the HTLV-1 estimates."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F           # noqa: E402
import joint_correction as J          # noqa: E402
import amplification_and_histology_control as A   # noqa: E402
import boundary_geometry_junction_exclusion as BOUNDARY       # noqa: E402

OUT = J.OUT
BIN = F.BIN_SIZE
COORDS = OUT / "htlv_coordinates.csv"
REFERENCE = ("HBV", "HPV")
HTLV_ARMS = ("HTLV_all", "HTLV_ATL", "HTLV_nonmalignant", "HTLV_unknown",
             "HTLV_ATL_s35654946", "HTLV_nonmalig_s35654946")


def asym_mde(o_num, n_num, o_den, n_den) -> dict:
    """Calculate the 80%-power minimum detectable displacement from the null mean.
    
    The signed asymmetry can cross zero, so the bound is an absolute displacement
    of 2.80 null standard deviations. It is calculated from the full-precision draws
    rather than a rounded z score."""
    tot_n = np.maximum(n_num + n_den, 1)
    a_n = (n_num - n_den) / tot_n
    a_o = (o_num - o_den) / max(o_num + o_den, 1)
    sd = float(a_n.std(ddof=1))
    return {"displacement": round(float(a_o - a_n.mean()), 6),
            "null_sd": round(sd, 6),
            "mde_displacement": round(2.80 * sd, 6)}


def load_all_coordinates() -> pd.DataFrame:
    co = F.load_coordinates()[["virus_class", "chrom", "pos0"]].copy()
    h = pd.read_csv(COORDS)
    # PMID 35654946 is the only study contributing to both clinical groups, so it carries the
    # only malignant-vs-non-malignant contrast in this dataset that is not confounded with
    # study. Both of its groups are emitted as their own arms.
    # PRIMARY arm: every HTLV-1 coordinate, pooled. The premise being tested is "does a virus
    # that cuts its own site carry the host break grammar", and that question does not need the
    # clinical split at all -- every record here is integrase-chosen whatever the donor's disease
    # state. Pooling is therefore the primary analysis and the maximum-power one; the clinical
    # arms below are a sensitivity analysis on whether clonal selection matters, and they are
    # confounded with study (see htlv_build.py).
    allarm = h.copy()
    allarm["arm"] = "HTLV_all"

    s = h[h["pmid"] == 35654946]
    extra = []
    for grp, name in (("HTLV_ATL", "HTLV_ATL_s35654946"),
                      ("HTLV_nonmalignant", "HTLV_nonmalig_s35654946")):
        e = s[s["arm"] == grp].copy()
        e["arm"] = name
        extra.append(e)
    h = pd.concat([h, allarm] + extra, ignore_index=True)
    htlv = h.rename(columns={"arm": "virus_class"})[["virus_class", "chrom", "pos0"]]
    print(f"[htlv]    {len(htlv):,} HTLV-1 coordinates "
          f"{htlv['virus_class'].value_counts().to_dict()}")
    return pd.concat([co, htlv], ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--viruses", default=",".join(REFERENCE + HTLV_ARMS))
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), BOUNDARY.WINDOW)
    base = F.Nuller()
    nul = J.JointNuller(feat, base, mask)
    geo = BOUNDARY.build_geometry(base)
    co = load_all_coordinates()

    rows = []
    for virus in a.viruses.split(","):
        sub = co[co["virus_class"] == virus].sort_values(
            ["chrom", "pos0"]).reset_index(drop=True)
        in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                            for c, p in zip(sub["chrom"], sub["pos0"])])
        sub = sub[~in_mask].reset_index(drop=True)
        k = pd.MultiIndex.from_arrays([sub["chrom"].astype(str).to_numpy(),
                                       (sub["pos0"].to_numpy(np.int64) // BIN) * BIN])
        sub = sub.assign(stratum=nul.lookup.reindex(k).to_numpy())
        sub = sub[sub["stratum"].notna()].reset_index(drop=True)
        sub["stratum"] = sub["stratum"].astype(int)
        sub["within_bin_offset"] = sub["pos0"] - (sub["pos0"] // BIN) * BIN

        labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                   sub["pos0"].to_numpy(np.int64), BOUNDARY.WINDOW)
        n_clu = int(len(np.unique(labels)))
        nul.n_calls = nul.n_fallback = 0
        rng = np.random.default_rng(BOUNDARY.SEED)
        obs, null = BOUNDARY.run_block(geo, nul, rng, sub, a.n_perm, labels)
        fb = round(100 * nul.n_fallback / max(nul.n_calls, 1), 2)
        meta = dict(virus=virus, window_kb=BOUNDARY.WINDOW // 1000, n_coords=len(sub),
                    n_clusters=n_clu, pct_pool_fallback=fb, n_perm=a.n_perm)

        o_int, n_int = obs[f"depth_{BOUNDARY.EDGE_BP}"], null[f"depth_{BOUNDARY.EDGE_BP}"]
        for jx in BOUNDARY.EXCL:
            o_edge = obs[f"depth_{jx}"] - o_int
            n_edge = null[f"depth_{jx}"] - n_int
            r = {**meta, "excl_kb": jx // 1000, "statistic": "edge_vs_interior",
                 "edge_fold": round(o_edge / max(float(n_edge.mean()), 1e-9), 4),
                 "interior_fold": round(o_int / max(float(n_int.mean()), 1e-9), 4),
                 **BOUNDARY.asym(o_edge, n_edge, o_int, n_int),
                 **asym_mde(o_edge, n_edge, o_int, n_int)}
            rows.append(r)
            print(f"  [{virus:18s} excl {jx//1000:>2} kb] edge {o_edge:>8} "
                  f"(fold {r['edge_fold']:.4f}) vs interior {o_int:>8} "
                  f"(fold {r['interior_fold']:.4f}) | A {r['observed']:+.4f} "
                  f"null {r['expected']:+.4f} z {r['z']:+.2f} p {r['emp_p']:.5g}", flush=True)
        for m in BOUNDARY.MARGINS:
            rows.append({**meta, "excl_kb": m // 1000,
                         "statistic": f"depth_margin_{m//1000}kb",
                         **F.summarize(obs[f"depth_{m}"], null[f"depth_{m}"])})

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "integrase_positive_geometry.csv", index=False)

    pub = Path(OUT / "boundary_junction_excluded.csv")
    if pub.exists():
        p = pd.read_csv(pub)
        key = ["virus", "excl_kb", "statistic"]
        m = d.merge(p[key + ["observed", "z"]], on=key, suffixes=("", "_pub"))
        m = m[m["virus"].isin(REFERENCE) & (m["statistic"] == "edge_vs_interior")]
        if len(m):
            dev = (m["observed"] - m["observed_pub"]).abs().max()
            print(f"\n[check]    {len(m)} reference rows; max |A - published| = {dev:.6f}")
            if dev > 1e-4:
                raise SystemExit("reproduction gate FAILED on the geometry table")

    pd.set_option("display.width", 260)
    print("\n=== edge_vs_interior by junction-exclusion radius (0 kb = published statistic) ===")
    print(d[d["statistic"] == "edge_vs_interior"][
        ["virus", "excl_kb", "n_clusters", "edge_fold", "interior_fold",
         "observed", "expected", "z", "emp_p"]].to_string(index=False))
    print(f"\nwrote {OUT/'integrase_positive_geometry.csv'}  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
