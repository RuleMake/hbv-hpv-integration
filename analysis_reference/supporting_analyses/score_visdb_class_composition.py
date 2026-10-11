"""Estimate structural-variant class composition at VISDB HBV and HPV coordinates.

The external coordinate sets are scored on the shared integration-negative host
reference and remain separate from the primary estimate. Class-atlas, composition
and leave-one-study-out results are written separately. Source-study heterogeneity
and the absence of patient identifiers limit patient-level sensitivity analyses."""

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
import locus_class_propensity as COMPOSITION  # noqa: E402

OUT = Path(r"source_project\results")
COORDS = OUT / "visdb_coordinates.csv"
WINDOW = 100_000
CLASSES = ["DUP", "DEL", "INV", "TRA"]
REFERENCE = ("HBV", "HPV")
VISDB_ARMS = ("VISDB_HBV", "VISDB_HPV")
SEED = 7
Z_FLOOR = 0.05


def mde(fold: float, z: float) -> float:
    if not np.isfinite(fold) or not np.isfinite(z) or abs(z) < Z_FLOOR:
        return np.inf
    return 1.0 + 2.80 * abs((fold - 1.0) / z)


def load_all_coordinates() -> pd.DataFrame:
    co = F.load_coordinates()[["virus_class", "chrom", "pos0"]].copy()
    co["pmid"] = np.nan
    v = pd.read_csv(COORDS).rename(columns={"arm": "virus_class"})
    print(f"[visdb]   {len(v):,} coordinates "
          f"{v['virus_class'].value_counts().to_dict()}")
    return pd.concat([co, v[["virus_class", "chrom", "pos0", "pmid"]]], ignore_index=True)


# ------------------------------------------------------------------ atlas
def run_atlas(allco, nul, mask, idx_by_class, viruses, n_perm):
    rows = []
    for cls in CLASSES:
        idx, n_ev = idx_by_class[cls]
        for virus in viruses:
            sub = ND.eligible(allco, nul, mask, virus)
            # Exclude and count coordinates whose matched pool contains no alternative position.
            _calls, _fb = nul.n_calls, nul.n_fallback
            sub, n_unscorable = F.drop_unscorable(nul, sub, label=f"{cls} {virus}")
            nul.n_calls, nul.n_fallback = _calls, _fb
            labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), WINDOW)
            n_clu = int(len(np.unique(labels)))
            rng = np.random.default_rng(SEED)
            o, nl, _, so, sn = F.est_block(idx, nul, rng, sub, WINDOW, n_perm,
                                           labels, unit="uniqSV")
            for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                             ("f_clusterSignTest", F.summarize(so, sn))):
                res["mde_fold"] = round(mde(res["fold"], res["z"]), 4)
                rows.append(dict(svclass=cls, n_events=n_ev, virus=virus, estimator=est,
                                 n_coords=len(sub), n_clusters=n_clu,
                                 n_unscorable_dropped=n_unscorable, **res))
            b = rows[-2]
            print(f"  [{cls:3s} {virus:11s}] n={b['n_coords']:6,}/{n_clu:6,}cl "
                  f"fold={b['fold']:.4f} [{b['fold_ci_lo']:.4f},{b['fold_ci_hi']:.4f}] "
                  f"p={b['emp_p']:.4g} MDE={b['mde_fold']:.4f}", flush=True)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ axis
def run_axis(allco, joint, mask, cum, nb, p_global, viruses, n_perm, window, seed=SEED):
    half = window // COMPOSITION.BIN
    dens = COMPOSITION.DensityNuller(joint, cum, nb, half)
    own_i = [COMPOSITION.CLASSES.index(c) for c in COMPOSITION.AXIS[0]]
    oth_i = [COMPOSITION.CLASSES.index(c) for c in COMPOSITION.AXIS[1]]
    rows, draws = [], {}
    for vi, vc in enumerate(viruses):
        sub0, sub, tier, labels = COMPOSITION.prepare(allco, vc, joint, dens, mask)
        # Exclude and count coordinates whose matched pool contains no alternative position.
        _calls, _fb = joint.n_calls, joint.n_fallback
        keep = np.ones(len(sub), dtype=bool)
        for i, (c, p, st) in enumerate(zip(sub["chrom"], sub["pos0"], sub["stratum"])):
            try:
                joint.pool(str(c), int(st), int(p) % COMPOSITION.BIN, int(p))
            except F.EmptyMatchedPool:
                keep[i] = False
        n_unscorable = int((~keep).sum())
        joint.n_calls, joint.n_fallback = _calls, _fb
        if n_unscorable:
            print(f"[unscorable] {vc}: dropped {n_unscorable:,} of {len(sub):,} coordinates with "
                  f"no matched position to draw against", flush=True)
            sub = sub.loc[keep].reset_index(drop=True)
            labels = COMPOSITION.F.leader_clusters(sub["chrom"].to_numpy(str),
                                           sub["pos0"].to_numpy(np.int64),
                                           COMPOSITION.CLUSTER_RADIUS) if hasattr(COMPOSITION, "CLUSTER_RADIUS")                 else F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), 100_000)
        n_clu = int(len(np.unique(labels)))
        print(f"\n[{window//1000} kb | {vc}] {len(sub0):,} -> {len(sub):,} coordinates, "
              f"{n_clu:,} clusters", flush=True)
        for arm, nuller, dt in (("stratum", joint, None), ("density", dens, tier)):
            rng = np.random.default_rng([seed, vi])
            nuller.n_calls = nuller.n_fallback = 0
            obs, obs_d, null, null_d, _ = COMPOSITION.score(sub, labels, nuller, cum, nb, half,
                                                    p_global, rng, n_perm, density_tier=dt)
            base = dict(virus=vc, arm=arm, window_bp=window, n_coords=len(sub),
                        n_clusters=n_clu, n_perm=n_perm,
                        null_fallback_pct=round(100 * nuller.n_fallback
                                                / max(nuller.n_calls, 1), 2))
            o_ax = obs[own_i].sum() - obs[oth_i].sum()
            n_ax = null[:, own_i].sum(axis=1) - null[:, oth_i].sum(axis=1)
            rows.append({**base, **COMPOSITION.summarize(o_ax, n_ax,
                                                 "class_axis_TRAINV_minus_DELDUP", signed=True)})
            rows.append({**base, **COMPOSITION.summarize(obs_d, null_d, "breakpoint_density")})
            draws[(arm, vc)] = (o_ax - n_ax.mean(), n_ax - n_ax.mean())
            print(f"  [{arm:7s}] axis {o_ax - n_ax.mean():+.5f}  z={rows[-2]['z']:+.2f} "
                  f"MDE={rows[-2]['mde_excess']:.5f} clears={rows[-2]['clears_mde']} "
                  f"fallback {base['null_fallback_pct']}%", flush=True)
    for arm in ("stratum", "density"):
        for i in range(len(viruses)):
            for k in range(i + 1, len(viruses)):
                va, vb = viruses[i], viruses[k]
                (ea, na), (eb, nbv) = draws[(arm, va)], draws[(arm, vb)]
                r = float(np.corrcoef(na, nbv)[0, 1])
                rows.append({"virus": f"{va}-vs-{vb}", "arm": arm, "window_bp": window,
                             "n_coords": np.nan, "n_clusters": np.nan, "n_perm": n_perm,
                             "null_fallback_pct": np.nan, "null_stream_corr": round(r, 4),
                             **COMPOSITION.summarize(ea - eb, na - nbv,
                                             f"interaction_{va}_minus_{vb}", signed=True)})
                print(f"  [{arm:7s}] {va} - {vb}: {ea - eb:+.5f}  "
                      f"z={rows[-1]['z']:+.2f} clears={rows[-1]['clears_mde']}", flush=True)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ leave-one-study-out
def run_loso(allco, joint, mask, cum, nb, p_global, n_perm, window, top_n) -> pd.DataFrame:
    half = window // COMPOSITION.BIN
    dens = COMPOSITION.DensityNuller(joint, cum, nb, half)
    own_i = [COMPOSITION.CLASSES.index(c) for c in COMPOSITION.AXIS[0]]
    oth_i = [COMPOSITION.CLASSES.index(c) for c in COMPOSITION.AXIS[1]]
    rows = []
    for arm_name in VISDB_ARMS:
        sel_all = allco[allco["virus_class"] == arm_name]
        top = sel_all["pmid"].value_counts().head(top_n).index
        for drop in top:
            sel = pd.concat([allco[allco["virus_class"] != arm_name],
                             sel_all[sel_all["pmid"] != drop]], ignore_index=True)
            _, sub, tier, labels = COMPOSITION.prepare(sel, arm_name, joint, dens, mask)
            rng = np.random.default_rng(SEED)
            obs, _, null, _, _ = COMPOSITION.score(sub, labels, joint, cum, nb, half,
                                           p_global, rng, n_perm, density_tier=None)
            o_ax = obs[own_i].sum() - obs[oth_i].sum()
            n_ax = null[:, own_i].sum(axis=1) - null[:, oth_i].sum(axis=1)
            rows.append(dict(virus=arm_name, dropped_pmid=int(drop), window_bp=window,
                             n_coords=len(sub), n_clusters=int(len(np.unique(labels))),
                             n_perm=n_perm,
                             **COMPOSITION.summarize(o_ax, n_ax, "class_axis_LOSO", signed=True)))
            print(f"  [LOSO {arm_name} -{int(drop)}] n={len(sub):,} "
                  f"excess={rows[-1]['excess']:+.5f} z={rows[-1]['z']:+.2f}", flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas-perm", type=int, default=20000)
    ap.add_argument("--axis-perm", type=int, default=10000)
    ap.add_argument("--loso-perm", type=int, default=4000)
    ap.add_argument("--loso-top", type=int, default=5)
    ap.add_argument("--window", type=int, default=1_000_000)
    ap.add_argument("--skip-atlas", action="store_true")
    ap.add_argument("--skip-loso", action="store_true")
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    joint = J.JointNuller(feat, F.Nuller(), mask)
    allco = load_all_coordinates()
    viruses = REFERENCE + VISDB_ARMS

    if not a.skip_atlas:
        canon = pd.read_csv(F.BASE / "input_pooled_canonical_events.tsv", sep="\t")
        led = ND.build_ledger(set(canon["aliquot_metadata"].astype(str)))
        idx_by_class = {c: ND.to_index(led[led.svclass == c]) for c in CLASSES}
        print("\n=== Class enrichment in Tian and VISDB coordinate sets ===")
        da = run_atlas(allco, joint, mask, idx_by_class, viruses, a.atlas_perm)
        da.to_csv(OUT / "visdb_atlas.csv", index=False)
        pub = pd.read_csv(OUT / "sv_class_locus_level.csv")
        key = ["svclass", "virus", "estimator"]
        m = da.merge(pub[key + ["fold"]], on=key, suffixes=("", "_pub"))
        m = m[m["virus"].isin(REFERENCE)]
        print(f"[check]    {len(m)} reference cells; "
              f"max |fold - published| = {(m['fold'] - m['fold_pub']).abs().max():.5f}")

    cum, nb, p_global = COMPOSITION.propensity_map()
    print("\n=== VISDB class-composition comparison ===")
    dx = run_axis(allco, joint, mask, cum, nb, p_global, viruses, a.axis_perm, a.window)
    dx.to_csv(OUT / "visdb_axis.csv", index=False)

    if not a.skip_loso:
        print("\n=== VISDB leave-one-study-out results (largest studies per group) ===")
        dl = run_loso(allco, joint, mask, cum, nb, p_global, a.loso_perm, a.window, a.loso_top)
        dl.to_csv(OUT / "visdb_loso.csv", index=False)

    pd.set_option("display.width", 240)
    cols = ["virus", "arm", "statistic", "observed", "excess", "z", "emp_p_two_sided",
            "mde_excess", "clears_mde"]
    print("\n" + dx[dx.statistic.str.startswith(("class_axis", "interaction"))][cols]
          .to_string(index=False))
    print(f"\n[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
