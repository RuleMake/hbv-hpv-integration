"""Estimate class-specific enrichment at HTLV-1 integration sites.

Clinical-state groups are analyzed separately, with unknown states retained as a
separate group. Within-study subsets of PMID 35654946 support comparisons that
are not confounded by between-study differences. Leave-one-study-out estimates
assess the contribution of individual source studies.

HBV, HPV and EBV reference results are recomputed through the same matching and
counting procedures and checked before the HTLV-1 comparisons are reported."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F          # noqa: E402
import joint_correction as J         # noqa: E402
import amplification_and_histology_control as A   # noqa: E402
import discriminating_tests as ND     # noqa: E402

PROJ = Path(r"source_project")
OUT = PROJ / "results"
COORDS = OUT / "htlv_coordinates.csv"
WINDOW = 100_000
CLASSES = ["DUP", "DEL", "INV", "TRA"]
REFERENCE = ("HBV", "HPV", "EBV")
HTLV_ARMS = ("HTLV_all", "HTLV_ATL", "HTLV_nonmalignant", "HTLV_unknown",
             "HTLV_ATL_s35654946", "HTLV_nonmalig_s35654946")
SEED = 7
Z_FLOOR = 0.05


def mde(fold: float, z: float) -> float:
    """1 + 2.80 |(fold-1)/z| -- the project's standard 80%-power minimum detectable fold
    (mde_audit.py). Infinite when z is at the floor: such a null carries no scale
    information and must be written as "not resolvable", never as "no effect"."""
    if not np.isfinite(fold) or not np.isfinite(z) or abs(z) < Z_FLOOR:
        return np.inf
    return 1.0 + 2.80 * abs((fold - 1.0) / z)


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
    htlv = h.rename(columns={"arm": "virus_class"})[["virus_class", "chrom", "pos0", "pmid"]]
    print(f"[htlv]    {len(htlv):,} HTLV-1 coordinates "
          f"{htlv['virus_class'].value_counts().to_dict()}")
    co["pmid"] = np.nan
    return pd.concat([co, htlv], ignore_index=True)


def atlas(allco, nul, mask, idx_by_class, n_ev_by_class, viruses, n_perm, seed=SEED):
    rows = []
    for cls in CLASSES:
        idx, n_ev = idx_by_class[cls], n_ev_by_class[cls]
        for virus in viruses:
            sub = ND.eligible(allco, nul, mask, virus)
            labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), WINDOW)
            n_clu = int(len(np.unique(labels)))
            rng = np.random.default_rng(seed)
            o, nl, _, so, sn = F.est_block(idx, nul, rng, sub, WINDOW, n_perm,
                                           labels, unit="uniqSV")
            for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                             ("f_clusterSignTest", F.summarize(so, sn))):
                res["mde_fold"] = round(mde(res["fold"], res["z"]), 4)
                rows.append(dict(svclass=cls, n_events=n_ev, virus=virus, estimator=est,
                                 n_coords=len(sub), n_clusters=n_clu, **res))
            b, s = rows[-2], rows[-1]
            print(f"  [{cls:3s} {virus:18s}] n={b['n_coords']:6,}/{n_clu:6,}cl  "
                  f"bc={b['fold']:.4f} [{b['fold_ci_lo']:.4f},{b['fold_ci_hi']:.4f}] "
                  f"p={b['emp_p']:.4g} MDE={b['mde_fold']:.3f} | "
                  f"sign={s['fold']:.4f} p={s['emp_p']:.4g}", flush=True)
    return pd.DataFrame(rows)


def check_reproduction(df: pd.DataFrame) -> None:
    """The reference viruses must reproduce the published atlas or the HTLV cells mean nothing."""
    pub = pd.read_csv(OUT / "sv_class_locus_level.csv")
    key = ["svclass", "virus", "estimator"]
    m = df.merge(pub[key + ["fold", "observed"]], on=key, suffixes=("", "_pub"))
    m = m[m["virus"].isin(REFERENCE)]
    if m.empty:
        print("[check]    No reference cells in this run; reference check unavailable")
        return
    bad_obs = m[m["observed"] != m["observed_pub"]]
    dev = (m["fold"] - m["fold_pub"]).abs()
    print(f"[check]    {len(m)} reference cells; max |fold - published| = {dev.max():.5f}; "
          f"observed-count mismatches = {len(bad_obs)}")
    if len(bad_obs):
        print(bad_obs[key + ["observed", "observed_pub"]].to_string(index=False))
        raise SystemExit("reproduction gate FAILED: observed counts differ from published atlas")
    if dev.max() > 0.01:
        raise SystemExit(f"reproduction gate FAILED: fold deviates by {dev.max():.4f}")


def loso(allco, nul, mask, idx_by_class, n_perm) -> pd.DataFrame:
    """Omit each contributing study in turn to assess study-level sensitivity in HTLV-1."""
    rows = []
    h = allco[allco["virus_class"].isin(HTLV_ARMS)]
    for arm in HTLV_ARMS:
        pmids = sorted(h.loc[h["virus_class"] == arm, "pmid"].dropna().unique())
        if len(pmids) < 2:
            # A single-study arm has no leave-one-study-out: dropping its only study leaves an
            # empty coordinate set, which is what crashed the first corrected run inside
            # ND.eligible. The within-study arms (s35654946) are single-study by construction.
            print(f"  [LOSO {arm:22s}] single study ({int(pmids[0]) if pmids else 'none'}) "
                  f"-- no leave-one-study-out is defined, skipped")
            continue
        for drop in pmids:
            sel = allco[(allco["virus_class"] == arm) & (allco["pmid"] != drop)].copy()
            for cls in CLASSES:
                sub = ND.eligible(sel, nul, mask, arm)
                labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                           sub["pos0"].to_numpy(np.int64), WINDOW)
                rng = np.random.default_rng(SEED)
                o, nl, _, so, sn = F.est_block(idx_by_class[cls], nul, rng, sub, WINDOW,
                                               n_perm, labels, unit="uniqSV")
                r = F.summarize(o, nl)
                r["mde_fold"] = round(mde(r["fold"], r["z"]), 4)
                rows.append(dict(virus=arm, dropped_pmid=int(drop), svclass=cls,
                                 n_coords=len(sub),
                                 n_clusters=int(len(np.unique(labels))), **r))
                print(f"  [LOSO {arm:18s} -{int(drop)} {cls:3s}] n={len(sub):6,} "
                      f"fold={r['fold']:.4f} p={r['emp_p']:.4g}", flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--loso-perm", type=int, default=5000)
    ap.add_argument("--skip-loso", action="store_true")
    ap.add_argument("--viruses", default=",".join(REFERENCE + HTLV_ARMS))
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    nul = J.JointNuller(feat, F.Nuller(), mask)

    canon = pd.read_csv(F.BASE / "input_pooled_canonical_events.tsv", sep="\t")
    led = ND.build_ledger(set(canon["aliquot_metadata"].astype(str)))
    idx_by_class, n_ev_by_class = {}, {}
    for cls in CLASSES:
        idx_by_class[cls], n_ev_by_class[cls] = ND.to_index(led[led.svclass == cls])

    allco = load_all_coordinates()
    viruses = tuple(a.viruses.split(","))

    print("\n=== Class enrichment: HBV, HPV and EBV compared with "
          "HTLV-1 ===")
    df = atlas(allco, nul, mask, idx_by_class, n_ev_by_class, viruses, a.n_perm)
    df.to_csv(OUT / "integrase_positive_atlas.csv", index=False)
    check_reproduction(df)

    if not a.skip_loso:
        print("\n=== HTLV-1 leave-one-study-out comparisons ===")
        dl = loso(allco, nul, mask, idx_by_class, a.loso_perm)
        dl.to_csv(OUT / "integrase_positive_loso.csv", index=False)

    cols = ["svclass", "virus", "estimator", "n_coords", "n_clusters", "fold",
            "fold_ci_lo", "fold_ci_hi", "z", "emp_p", "mde_fold"]
    print("\n" + df[df.estimator == "bc_uniqSV_blockNull"][cols].to_string(index=False))
    print(f"\nwrote {OUT / 'integrase_positive_atlas.csv'}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
