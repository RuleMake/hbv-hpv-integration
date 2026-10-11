"""Evaluate structural-variant classes and reference replication-initiation zones.

The class analysis rebuilds class-specific variant tables from a common aliquot
set. Intrachromosomal classes use intervals; translocations use endpoints.
The initiation-zone analysis uses fork-directionality gradients from the HeLa
OK-seq reference. It describes the location of integration coordinates relative
to that reference, rather than origin activity in individual tumours.

Both analyses use joint genomic matching, amplification-peak exclusion and rigid
cluster placement, with 20,000 permutations, seed 7 and three virus coordinate sets."""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402
import joint_correction as J  # noqa: E402
import amplification_and_histology_control as A  # noqa: E402
from minimal_bigwig_reader import BigWig  # noqa: E402

SRC = Path(r"upstream_project")
OUT = J.OUT
CACHE = OUT / "cache" / "sv_class_ledger.parquet"
SV_DIRS = [SRC / "data/raw/pcawg_sv/icgc_open", SRC / "data/raw/pcawg_sv/tcga_open"]
RFD_BW = SRC / "data/raw/okseq/Petryk2016_HeLa_OkaSeq_RFD_hg19.bw"
BIN = F.BIN_SIZE
RES = 50_000               # RFD native resolution
SLOPE_HALF = 4             # +/- 4 bins = +/- 200 kb for the local RFD slope
SEED = 7
CLASSES = {"DUP": ("DUP",), "DEL": ("DEL",), "INV": ("H2HINV", "T2TINV"), "TRA": ("TRA",)}
SPAN = (100_000, 7_500_000)


# Structural-variant class catalogue
def build_ledger(aliquots: set[str]) -> pd.DataFrame:
    if CACHE.exists():
        d = pd.read_parquet(CACHE)
        print(f"[ledger]  cache hit: {len(d):,} events, classes {d.svclass.value_counts().to_dict()}")
        return d
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    paths = {f.name.split(".")[0]: f for dd in SV_DIRS for f in dd.rglob("*.bedpe.gz")}
    use = {a: paths[a] for a in aliquots if a in paths}
    print(f"[ledger]  parsing {len(use)} canonical aliquots for DUP/DEL/INV/TRA")
    keep = []
    for i, (al, p) in enumerate(sorted(use.items()), 1):
        with gzip.open(p, "rt") as fh:
            d = pd.read_csv(fh, sep="\t", low_memory=False)
        if d.empty:
            continue
        cls = d["svclass"].astype(str).str.upper()
        same = d["chrom1"].astype(str) == d["chrom2"].astype(str)
        c1 = "chr" + d["chrom1"].astype(str).str.replace("chr", "", regex=False)
        c2 = "chr" + d["chrom2"].astype(str).str.replace("chr", "", regex=False)
        s1, s2 = d["start1"].to_numpy(np.int64), d["start2"].to_numpy(np.int64)
        for name, labs in CLASSES.items():
            m = cls.isin(labs).to_numpy()
            if name == "TRA":
                # a translocation has no interval: score each endpoint as its own zero-length event
                mm = m & (~same).to_numpy()
                if mm.any():
                    for cc, ss in ((c1[mm], s1[mm]), (c2[mm], s2[mm])):
                        keep.append(pd.DataFrame({"aliquot": al, "svclass": name,
                                                  "chrom": cc.to_numpy(),
                                                  "left_pos0": ss, "right_pos0": ss}))
                continue
            mm = m & same.to_numpy()
            if not mm.any():
                continue
            lo, hi = np.minimum(s1[mm], s2[mm]), np.maximum(s1[mm], s2[mm])
            sp = hi - lo
            k = (sp > SPAN[0]) & (sp <= SPAN[1])
            if k.any():
                keep.append(pd.DataFrame({"aliquot": al, "svclass": name,
                                          "chrom": c1[mm].to_numpy()[k],
                                          "left_pos0": lo[k], "right_pos0": hi[k]}))
        if i % 400 == 0:
            print(f"          {i}/{len(use)}", flush=True)
    out = pd.concat(keep, ignore_index=True)
    out = out[out["chrom"].isin([f"chr{c}" for c in list(range(1, 23)) + ["X"]])]
    out["unique_td_id"] = (out["aliquot"] + ":" + out["svclass"] + ":" + out["chrom"] + ":"
                           + out["left_pos0"].astype(str) + "-" + out["right_pos0"].astype(str))
    out = out.drop_duplicates("unique_td_id").reset_index(drop=True)
    out.to_parquet(CACHE, index=False)
    print(f"[ledger]  built {len(out):,} events {out.svclass.value_counts().to_dict()}")
    return out


def to_index(ev: pd.DataFrame):
    ev = ev.sort_values(["chrom", "left_pos0"])
    codes = pd.factorize(ev["unique_td_id"])[0].astype(np.int64)
    ev = ev.assign(_eid=codes)
    return {str(c): F.ChromEvents(g["left_pos0"].to_numpy(np.int64),
                                  g["right_pos0"].to_numpy(np.int64),
                                  g["_eid"].to_numpy(np.int64))
            for c, g in ev.groupby("chrom")}, len(ev)


# Replication-initiation arrays
def initiation_arrays(chrom_sizes: dict[str, int]):
    """Per-chromosome 50-kb arrays: RFD, its local slope, and an initiation-like indicator."""
    bw = BigWig(RFD_BW)
    v = bw.validate()
    print(f"[rfd]     bigWig validate rel_diff {v['rel_diff']:.3g} passed={v['passed']}")
    binned, _, _ = bw.binned_means(RES)
    bw.close()
    out = {}
    tot_init = tot_bins = 0
    for c, size in chrom_sizes.items():
        if c not in binned:
            continue
        _, m = binned[c]
        r = np.nan_to_num(m, nan=0.0).astype(np.float32)
        n = len(r)
        k = SLOPE_HALF
        pad = np.pad(r, k, mode="edge")
        # local slope over +/- SLOPE_HALF bins: (downstream mean) - (upstream mean)
        down = np.array([pad[i + k + 1:i + 2 * k + 1].mean() for i in range(n)], dtype=np.float32)
        up = np.array([pad[i:i + k].mean() for i in range(n)], dtype=np.float32)
        slope = down - up
        # initiation-like: RFD rises through the middle AND is not already saturated either side
        init = ((slope > 0.15) & (up < 0.1) & (down > -0.1)).astype(np.int8)
        term = ((slope < -0.15) & (up > -0.1) & (down < 0.1)).astype(np.int8)
        out[c] = dict(rfd=r, slope=slope, init=init, term=term, n=n)
        tot_init += int(init.sum())
        tot_bins += n
    print(f"[rfd]     initiation-like bins: {tot_init}/{tot_bins} "
          f"({100*tot_init/max(tot_bins,1):.1f}% of the genome at {RES//1000} kb)")
    return out


IK = ("init", "term", "slope_pos", "abs_slope")


def init_stats(arr: dict, pos: np.ndarray) -> dict[str, np.ndarray]:
    q = pos // RES
    ok = (q >= 0) & (q < arr["n"])
    qq = np.clip(q, 0, arr["n"] - 1)
    z = ok.astype(np.float64)
    return {"init": arr["init"][qq] * z, "term": arr["term"][qq] * z,
            "slope_pos": (arr["slope"][qq] > 0).astype(np.float64) * z,
            "abs_slope": np.abs(arr["slope"][qq]) * z}


def run_init(arrs, nul, rng, sub, n_perm, labels):
    df = sub.copy()
    df["_cluster"] = labels
    obs = {k: 0.0 for k in IK}
    null = {k: np.zeros(n_perm) for k in IK}
    for _, g in df.groupby("_cluster", sort=True):
        c = str(g["chrom"].iloc[0])
        arr = arrs.get(c)
        if arr is None:
            continue
        p = g["pos0"].to_numpy(np.int64)
        st, off = g["stratum"].to_numpy(), g["within_bin_offset"].to_numpy()
        j = rng.integers(0, len(p), size=n_perm)
        delta = np.empty(n_perm, dtype=np.int64)
        for k in range(len(p)):
            m = j == k
            nk = int(m.sum())
            if nk:
                delta[m] = nul.draw(rng, c, int(st[k]), int(off[k]), int(p[k]), nk) - int(p[k])
        for pi in p:
            o = init_stats(arr, np.array([pi], dtype=np.int64))
            nl = init_stats(arr, pi + delta)
            for k in IK:
                obs[k] += float(o[k][0])
                null[k] += nl[k]
    return obs, null


def eligible(co, nul, mask, virus, win=100_000):
    sub = co[co["virus_class"] == virus].sort_values(["chrom", "pos0"]).reset_index(drop=True)
    in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                        for c, p in zip(sub["chrom"], sub["pos0"])])
    sub = sub[~in_mask].reset_index(drop=True)
    k = pd.MultiIndex.from_arrays([sub["chrom"].astype(str).to_numpy(),
                                   (sub["pos0"].to_numpy(np.int64) // BIN) * BIN])
    sub = sub.assign(stratum=nul.lookup.reindex(k).to_numpy())
    sub = sub[sub["stratum"].notna()].reset_index(drop=True)
    sub["stratum"] = sub["stratum"].astype(int)
    sub["within_bin_offset"] = sub["pos0"] - (sub["pos0"] // BIN) * BIN
    return sub


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--viruses", default="HBV,HPV,EBV")
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), 100_000)
    base = F.Nuller()
    nul = J.JointNuller(feat, base, mask)
    co = F.load_coordinates()
    viruses = a.viruses.split(",")

    # Class-specific spatial enrichment
    canon = pd.read_csv(F.BASE / "input_pooled_canonical_events.tsv", sep="\t")
    led = build_ledger(set(canon["aliquot_metadata"].astype(str)))
    rows = []
    print("\n=== Joint-matched class-specific spatial enrichment ===")
    for cls in CLASSES:
        idx, n_ev = to_index(led[led.svclass == cls])
        for virus in viruses:
            sub = eligible(co, nul, mask, virus)
            labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), 100_000)
            rng = np.random.default_rng(SEED)
            o, nl, _, so, sn = F.est_block(idx, nul, rng, sub, 100_000, a.n_perm,
                                           labels, unit="uniqSV")
            for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                             ("f_clusterSignTest", F.summarize(so, sn))):
                rows.append(dict(svclass=cls, n_events=n_ev, virus=virus, estimator=est,
                                 n_coords=len(sub), **res))
            print(f"  [{cls} {virus}] n_ev={n_ev:6d} bc={rows[-2]['fold']:.4f} "
                  f"[{rows[-2]['fold_ci_lo']:.4f},{rows[-2]['fold_ci_hi']:.4f}] "
                  f"p={rows[-2]['emp_p']:.5g} | sign={rows[-1]['fold']:.4f} "
                  f"p={rows[-1]['emp_p']:.5g}", flush=True)
    d37 = pd.DataFrame(rows)
    d37.to_csv(OUT / "sv_class_locus_level.csv", index=False)

    # Reference initiation-zone comparison
    arrs = initiation_arrays(base.chrom_size)
    rows = []
    print("\n=== Integration coordinates and reference replication-initiation zones ===")
    for virus in viruses:
        sub = eligible(co, nul, mask, virus)
        labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                   sub["pos0"].to_numpy(np.int64), 100_000)
        rng = np.random.default_rng(SEED)
        obs, null = run_init(arrs, nul, rng, sub, a.n_perm, labels)
        for k in IK:
            mu = float(null[k].mean())
            sd = float(null[k].std(ddof=1))
            p = (int(np.sum(np.abs(null[k] - mu) >= abs(obs[k] - mu))) + 1) / (a.n_perm + 1)
            fold = obs[k] / mu if mu else np.nan
            rows.append(dict(virus=virus, statistic=k, n_coords=len(sub),
                             observed=round(obs[k], 2), expected=round(mu, 2),
                             fold=round(fold, 4) if np.isfinite(fold) else np.nan,
                             z=round((obs[k] - mu) / sd, 2) if sd else np.nan,
                             emp_p=round(p, 5),
                             mde_fold_80pct=round(1 + 2.8 * sd / mu, 3) if mu else np.inf))
        print(f"  [{virus}] init {rows[-4]['fold']:.4f} (p={rows[-4]['emp_p']:.4g}, "
              f"MDE {rows[-4]['mde_fold_80pct']}) | term {rows[-3]['fold']:.4f} "
              f"(p={rows[-3]['emp_p']:.4g})", flush=True)
    d38 = pd.DataFrame(rows)
    d38.to_csv(OUT / "initiation_zone_test.csv", index=False)

    pd.set_option("display.width", 250)
    print("\n=== Class-specific spatial enrichment ===")
    print(d37[["svclass", "n_events", "virus", "estimator", "fold", "fold_ci_lo", "fold_ci_hi",
               "emp_p"]].to_string(index=False))
    print("\n=== Reference initiation zones ===")
    print(d38.to_string(index=False))
    print(f"\nwrote both ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
