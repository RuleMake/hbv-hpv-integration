"""Computes structural-variant class composition using a germline-variant reference.
"""

from __future__ import annotations

import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F                    # noqa: E402
import joint_correction as J                   # noqa: E402
import amplification_and_histology_control as A    # noqa: E402
import discriminating_tests as N                # noqa: E402

PROJ = Path(r"source_project")
DGV = PROJ / "data" / "external_sv" / "dgvMerged_hg19.txt.gz"
OUT = PROJ / "results" / "probe_germline_grammar.csv"

CLUSTER_WIN = 100_000                       # leader clustering, as everywhere else
WINDOWS = (200_000, 1_000_000, 2_000_000)
VIRUSES = ("HBV", "HPV")
DEL_TYPES = {"deletion", "loss"}
DUP_TYPES = {"duplication", "gain", "tandem duplication"}
P3 = "1000 Genomes Consortium Phase 3"


def load_dgv():
    """Return {track: {chrom: (sorted breakpoint positions, is_dup flag)}}."""
    rec = {"p3": [], "all": []}
    with gzip.open(DGV, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vt = f[10]
            if vt in DEL_TYPES:
                dup = 0
            elif vt in DUP_TYPES:
                dup = 1
            else:
                continue
            chrom, start, end, study = f[1], int(f[2]), int(f[3]), f[11]
            rec["all"].append((chrom, start, end, dup))
            if study == P3:
                rec["p3"].append((chrom, start, end, dup))

    tracks = {}
    for k, v in rec.items():
        d = pd.DataFrame(v, columns=["chrom", "start", "end", "dup"])
        # each record contributes two breakpoints, mirroring the somatic estimator
        bp = pd.concat([d[["chrom", "dup"]].assign(pos=d["start"]),
                        d[["chrom", "dup"]].assign(pos=d["end"])], ignore_index=True)
        bp = bp.sort_values(["chrom", "pos"]).reset_index(drop=True)
        tracks[k] = {c: (g["pos"].to_numpy(np.int64), g["dup"].to_numpy(np.int8))
                     for c, g in bp.groupby("chrom", sort=False)}
        print(f"[dgv]    {k:4s} {len(d):>7,} records -> {len(bp):>7,} breakpoints, "
              f"{d['dup'].mean():.3f} duplication share, {d['chrom'].nunique()} chromosomes")
    return tracks


def window_counts(track, chrom, pos, half):
    """(n_del, n_dup) germline breakpoints within +/- half of each position."""
    n_del = np.zeros(len(pos), np.int64)
    n_dup = np.zeros(len(pos), np.int64)
    for c in pd.unique(chrom):
        m = chrom == c
        if c not in track:
            continue
        p, dup = track[c]
        cum = np.concatenate([[0], np.cumsum(dup)])          # cumulative dup count
        lo = np.searchsorted(p, pos[m] - half, "left")
        hi = np.searchsorted(p, pos[m] + half, "right")
        tot = hi - lo
        nd = cum[hi] - cum[lo]
        n_dup[m] = nd
        n_del[m] = tot - nd
    return n_del, n_dup


def main() -> None:
    t0 = time.time()
    tracks = load_dgv()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), CLUSTER_WIN)
    nul = J.JointNuller(feat, F.Nuller(), mask)
    co = F.load_coordinates()

    rows = []
    for track_name in ("p3", "all"):
        track = tracks[track_name]
        for win in WINDOWS:
            half = win // 2
            per_virus = {}
            for vc in VIRUSES:
                sub = N.eligible(co, nul, mask, vc)
                labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                           sub["pos0"].to_numpy(np.int64), CLUSTER_WIN)
                # one row per cluster: leader coordinate, its stratum, its chromosome
                first = pd.DataFrame({"cluster": labels,
                                      "chrom": sub["chrom"].astype(str).to_numpy(),
                                      "pos0": sub["pos0"].to_numpy(np.int64),
                                      "stratum": sub["stratum"].to_numpy()}
                                     ).groupby("cluster", sort=True).first().reset_index()
                nd, nu = window_counts(track, first["chrom"].to_numpy(str),
                                       first["pos0"].to_numpy(np.int64), half)
                tot = nd + nu
                first["n_del"], first["n_dup"], first["n_tot"] = nd, nu, tot
                first["density"] = tot / (win / 1e6)                 # breakpoints per Mb
                ok = tot > 0
                first["prop_dup"] = np.where(ok, nu / np.maximum(tot, 1), np.nan)
                per_virus[vc] = first

                # ---- Q1 dynamic range: observed spread vs binomial null at the same depth
                v = first.loc[ok]
                pbar = v["n_dup"].sum() / v["n_tot"].sum()
                obs_sd = float(v["prop_dup"].std())
                exp_sd = float(np.sqrt((pbar * (1 - pbar) / v["n_tot"]).mean()))
                # ---- Q2 absorption by chromosome x joint stratum
                g = v.groupby(["chrom", "stratum"])["prop_dup"]
                resid = v["prop_dup"] - g.transform("mean")
                absorbed = 1 - float(resid.var()) / float(v["prop_dup"].var())
                dg = v.groupby(["chrom", "stratum"])["density"]
                d_absorbed = 1 - float((v["density"] - dg.transform("mean")).var()) \
                    / float(v["density"].var())

                rows.append(dict(
                    track=track_name, window_bp=win, virus=vc, gate="Q1_Q2",
                    n_clusters=len(first), n_with_data=int(ok.sum()),
                    pct_with_data=round(100 * ok.mean(), 1),
                    median_n_tot=float(v["n_tot"].median()),
                    mean_prop_dup=round(float(v["prop_dup"].mean()), 5),
                    mean_density=round(float(v["density"].mean()), 3),
                    sd_obs=round(obs_sd, 5), sd_binomial=round(exp_sd, 5),
                    heterogeneity_ratio=round(obs_sd / exp_sd, 3),
                    absorbed_prop_dup=round(absorbed, 4),
                    absorbed_density=round(d_absorbed, 4)))
                print(f"[{track_name} {win//1000:>4}kb {vc}] {int(ok.sum()):>5}/{len(first)} clusters "
                      f"with data, median {v['n_tot'].median():.0f} bp; "
                      f"prop_dup {v['prop_dup'].mean():.4f} (sd {obs_sd:.4f} vs binom {exp_sd:.4f} "
                      f"= {obs_sd/exp_sd:.2f}x); stratum absorbs {100*absorbed:.0f}% of prop_dup, "
                      f"{100*d_absorbed:.0f}% of density", flush=True)

            # ---- Q3 raw HPV - HBV difference, no null
            h, p = per_virus["HBV"], per_virus["HPV"]
            hv, pv = h[h["n_tot"] > 0], p[p["n_tot"] > 0]
            d_prop = float(pv["prop_dup"].mean() - hv["prop_dup"].mean())
            d_dens = float(pv["density"].mean() - hv["density"].mean())
            # crude two-sample scale for reference only -- the real test needs the block null
            se_prop = float(np.hypot(pv["prop_dup"].sem(), hv["prop_dup"].sem()))
            se_dens = float(np.hypot(pv["density"].sem(), hv["density"].sem()))
            rows.append(dict(
                track=track_name, window_bp=win, virus="HPV-minus-HBV", gate="Q3",
                n_clusters=len(pv) + len(hv), n_with_data=len(pv) + len(hv),
                mean_prop_dup=round(d_prop, 5), mean_density=round(d_dens, 3),
                naive_z_prop_dup=round(d_prop / se_prop, 2) if se_prop else np.nan,
                naive_z_density=round(d_dens / se_dens, 2) if se_dens else np.nan,
                predicted_sign="positive"))
            print(f"[{track_name} {win//1000:>4}kb  Q3] HPV - HBV: prop_dup {d_prop:+.5f} "
                  f"(naive z {d_prop/se_prop:+.2f}), density {d_dens:+.3f}/Mb "
                  f"(naive z {d_dens/se_dens:+.2f})   [predicted: positive]\n", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    pd.set_option("display.width", 260)
    print(d.to_string(index=False))
    print(f"\nwrote {OUT}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
