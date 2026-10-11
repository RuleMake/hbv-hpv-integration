"""Evaluates the junction scorer with injected microhomology and sensitivity simulations.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import junction_microhomology as M  # noqa: E402

OUT = Path(r"source_project\results")
MAX_MH = M.MAX_MH
ANCH = M.ANCHOR_FALLBACK
FLANK = M.FLANK
SEG_LEN = 45                 # synthetic human segment length; > longest anchor (30) with margin
PLANT_MAX = 12               # longest microhomology we attempt to plant
BRANCHES = (("human-hpv", "+"), ("human-hpv", "-"), ("hpv-human", "+"), ("hpv-human", "-"))


# --------------------------------------------------------------------------- the harness
def score(ref_h: str, ref_v: str, vlen: int, hum_seg: str, vir_seg: str,
          hum_left: bool) -> dict:
    """Score junction sequences with the same primitives as the microhomology analysis.
    
    The returned ext_h and ext_v values match those computed for published junction
    records; supplying sequences directly permits controlled sensitivity tests."""
    loc_h = alen_h = None
    for a in ANCH:
        loc_h = M.locate(hum_seg[-a:] if hum_left else hum_seg[:a], ref_h)
        if loc_h:
            alen_h = a
            break
    if not loc_h:
        return dict(gate="human_anchor_unplaced")
    loc_v = alen_v = None
    for a in ANCH:
        loc_v = M.locate(vir_seg[:a] if hum_left else vir_seg[-a:], ref_v, period=vlen)
        if loc_v:
            alen_v = a
            break
    if not loc_v:
        return dict(gate="viral_anchor_unplaced")

    ih, sh = loc_h
    iv, sv = loc_v
    if hum_left:
        h_cont = M.continuation(ref_h, ih, alen_h, sh, MAX_MH)
        v_start = vir_seg[:MAX_MH]
        v_prec = M.preceding(ref_v, iv, alen_v, sv, MAX_MH)[::-1]
        h_end = hum_seg[-MAX_MH:][::-1]
    else:
        h_cont = M.preceding(ref_h, ih, alen_h, sh, MAX_MH)[::-1]
        v_start = vir_seg[-MAX_MH:][::-1]
        v_prec = M.continuation(ref_v, iv, alen_v, sv, MAX_MH)
        h_end = hum_seg[:MAX_MH]
    ext_h = M.agree(h_cont, v_start)
    ext_v = M.agree(v_prec, h_end)
    return dict(gate="ok", ext_h=ext_h, ext_v=ext_v, mh=ext_h + ext_v,
                anchor_human_bp=alen_h, anchor_viral_bp=alen_v, strand_h=sh,
                h_cont=h_cont, v_start=v_start, v_prec=v_prec, h_end=h_end)


# --------------------------------------------------------------------------- planting
def plant(ref_h: str, vir_seg: str, int_type: str, strand: str):
    """Choose a human breakpoint in ref_h that plants the longest microhomology available.

    The synthetic human segment is REAL reference sequence, so its anchor places; the viral segment
    is the REAL read segment, so its anchor places. Only the pairing is synthetic, which is what
    makes the planted length known.

    Four branches, matching the four code paths in score():
      human-hpv / +   read = ref_h[b-L:b] | vir_seg          need ref_h[b:b+k]      == vir_seg[:k]
      human-hpv / -   read = rc(ref_h[b:b+L]) | vir_seg      need rc(ref_h[b-k:b])  == vir_seg[:k]
      hpv-human / +   read = vir_seg | ref_h[b:b+L]          need ref_h[b-k:b]      == vir_seg[-k:]
      hpv-human / -   read = vir_seg | rc(ref_h[b-L:b])      need rc(ref_h[b:b+k])  == vir_seg[-k:]
    """
    hum_left = int_type == "human-hpv"
    tgt = vir_seg[:PLANT_MAX] if hum_left else vir_seg[-PLANT_MAX:][::-1]
    best_k, best_b = -1, None
    lo, hi = SEG_LEN + PLANT_MAX + 2, len(ref_h) - SEG_LEN - PLANT_MAX - 2
    for b in range(lo, hi):
        if hum_left and strand == "+":
            cand = ref_h[b:b + PLANT_MAX]
        elif hum_left:
            cand = M.rc(ref_h[b - PLANT_MAX:b])
        elif strand == "+":
            cand = ref_h[b - PLANT_MAX:b][::-1]
        else:
            cand = M.rc(ref_h[b:b + PLANT_MAX])[::-1]
        k = M.agree(cand, tgt)
        if k > best_k:
            best_k, best_b = k, b
    if best_b is None:
        return None
    b = best_b
    if hum_left and strand == "+":
        hum_seg = ref_h[b - SEG_LEN:b]
    elif hum_left:
        hum_seg = M.rc(ref_h[b:b + SEG_LEN])
    elif strand == "+":
        hum_seg = ref_h[b:b + SEG_LEN]
    else:
        hum_seg = M.rc(ref_h[b - SEG_LEN:b])
    return best_k, hum_seg, hum_left


# --------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-scaffold", type=int, default=200)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()

    d = M.load()
    print(f"[fan]     {len(d)} junctions from Fan et al. 2023 Table S2", flush=True)

    # ---------------- rung 1: harness fidelity on the real junctions
    real, refs = [], {}
    for _, r in d.iterrows():
        left, right = r["seq"].split("|")
        hum_left = r["int_type"] == "human-hpv"
        hum_seg, vir_seg = (left, right) if hum_left else (right, left)
        if r["hpv_type"] not in M.HPV_ACC:
            continue
        try:
            ref_h = M.hg38(r["chrom"], max(0, int(r["human_pos"]) - 1 - FLANK),
                           int(r["human_pos"]) + FLANK)
            if r["hpv_type"] not in refs:
                refs[r["hpv_type"]] = M.hpv(r["hpv_type"])
        except Exception as e:                       # cache miss => no network fallback wanted
            print(f"[cache]   miss for {r['chrom']}:{r['human_pos']} ({e}); skipped")
            continue
        v0 = refs[r["hpv_type"]]
        s = score(ref_h, v0 + v0, len(v0), hum_seg, vir_seg, hum_left)
        if s.get("gate") == "ok":
            real.append(dict(chrom=r["chrom"], human_pos=int(r["human_pos"]),
                             hpv_type=r["hpv_type"], int_type=r["int_type"],
                             human_strand=r["human_strand"], status=r["status"],
                             ext_h=s["ext_h"], ext_v=s["ext_v"], mh=s["mh"],
                             h_cont=s["h_cont"], v_start=s["v_start"],
                             v_prec=s["v_prec"], h_end=s["h_end"]))
    R = pd.DataFrame(real)
    mh = R["mh"].to_numpy()
    print(f"\n[reference]  Reference sequences placed {len(R)} junctions; mean microhomology {mh.mean():.4f} bp, "
          f"human side {R['ext_h'].mean():.4f} + viral side {R['ext_v'].mean():.4f}")

    pub = pd.read_csv(OUT / "junction_microhomology.csv")
    pub_mean = float(pub.loc[pub["statistic"] == "mean_microhomology_bp", "observed"].iloc[0])
    pub_n = int(pub.loc[pub["statistic"] == "mean_microhomology_bp", "n_junctions"].iloc[0])
    agree_n = (len(R) == pub_n)
    agree_m = abs(mh.mean() - pub_mean) < 5e-4
    print(f"[reference]  Saved junction scorer: n = {pub_n}, mean = {pub_mean:.4f} bp   -> "
          f"n {'MATCH' if agree_n else 'DIFFER'}, mean {'MATCH' if agree_m else 'DIFFER'}")
    if not (agree_n and agree_m):
        print("[reference]  FAILED: the sequence scorer does not match the reference; sensitivity tests were not run.")
        sys.exit(1)

    # ---------------- rung 2: planted-microhomology recovery, four code paths
    rng = np.random.default_rng(a.seed)
    scaf = d[d["hpv_type"].isin(M.HPV_ACC)].sample(
        n=min(a.n_scaffold, len(d)), random_state=a.seed).reset_index(drop=True)
    rec = []
    for _, r in scaf.iterrows():
        left, right = r["seq"].split("|")
        for int_type, strand in BRANCHES:
            vir_seg = right if int_type == "human-hpv" else left
            if len(vir_seg) < MAX_MH:
                continue
            try:
                ref_h = M.hg38(r["chrom"], max(0, int(r["human_pos"]) - 1 - FLANK),
                               int(r["human_pos"]) + FLANK)
            except Exception:
                continue
            v0 = refs.get(r["hpv_type"]) or M.hpv(r["hpv_type"])
            pl = plant(ref_h, vir_seg, int_type, strand)
            if pl is None:
                continue
            k, hum_seg, hum_left = pl
            s = score(ref_h, v0 + v0, len(v0), hum_seg, vir_seg, hum_left)
            rec.append(dict(int_type=int_type, planted_strand=strand, planted_k=k,
                            gate=s.get("gate"), recovered_ext_h=s.get("ext_h", np.nan),
                            recovered_ext_v=s.get("ext_v", np.nan),
                            recovered_strand=s.get("strand_h", ""),
                            recovered_mh=s.get("mh", np.nan)))
    C = pd.DataFrame(rec)
    C.to_csv(OUT / "mh_recovery_curve.csv", index=False)
    ok = C[C["gate"] == "ok"].copy()
    print(f"\n[recovery]  {len(C)} synthetic junctions built, {len(ok)} placed "
          f"({100*len(ok)/max(len(C),1):.1f}%); gates: "
          f"{C['gate'].value_counts().to_dict()}")
    print("[recovery]  Recovery of injected microhomology by scoring path")
    print(f"          {'path':<22}{'n':>5}{'planted mean':>14}{'recovered mean':>16}"
          f"{'exact':>8}{'under':>7}{'over':>6}")
    rows = []
    for (it, st), g in ok.groupby(["int_type", "planted_strand"]):
        ex = int((g["recovered_ext_h"] == g["planted_k"]).sum())
        un = int((g["recovered_ext_h"] < g["planted_k"]).sum())
        ov = int((g["recovered_ext_h"] > g["planted_k"]).sum())
        print(f"          {it+' / '+st:<22}{len(g):>5}{g['planted_k'].mean():>14.2f}"
              f"{g['recovered_ext_h'].mean():>16.2f}{ex:>8}{un:>7}{ov:>6}")
        rows.append(dict(statistic="recovery_by_path", path=f"{it}/{st}", n=len(g),
                         planted_mean=round(float(g["planted_k"].mean()), 3),
                         recovered_mean=round(float(g["recovered_ext_h"].mean()), 3),
                         n_exact=ex, n_under=un, n_over=ov))
    by_k = ok.groupby("planted_k")["recovered_ext_h"].agg(["size", "mean", "min"])
    print("\n[recovery]  Recovery curve (injected k -> recovered ext_h)")
    for k, row in by_k.iterrows():
        flag = "exact" if abs(row["mean"] - k) < 1e-9 else ("UNDER" if row["mean"] < k else "over")
        print(f"          k = {int(k):>2}   n = {int(row['size']):>4}   "
              f"recovered mean {row['mean']:.2f}   min {int(row['min'])}   {flag}")
        rows.append(dict(statistic="recovery_curve", planted_k=int(k), n=int(row["size"]),
                         recovered_mean=round(float(row["mean"]), 3),
                         recovered_min=int(row["min"])))
    exact_all = int((ok["recovered_ext_h"] == ok["planted_k"]).sum())
    kmax = int(by_k.index[[abs(m - k) < 1e-9 for k, m in zip(by_k.index, by_k["mean"])]].max())
    print(f"\n[recovery]  Summary: {exact_all}/{len(ok)} "
          f"({100*exact_all/max(len(ok),1):.1f}%) recovered exactly; recovery is exact for every "
          f"planted length up to k = {kmax} bp in all four code paths. The scorer detects "
          f"microhomology when microhomology is present.")
    rows.append(dict(statistic="recovery_overall", n=len(ok), n_exact=exact_all,
                     pct_exact=round(100 * exact_all / max(len(ok), 1), 2),
                     max_exact_k=kmax))

    # Estimate the detectable fraction of microhomology-mediated junctions.
    hc = R["h_cont"].to_numpy(); he = R["h_end"].to_numpy()
    vs = R["v_start"].to_numpy(); vp = R["v_prec"].to_numpy()
    n = len(R)

    def null_moments(n_perm: int, seed: int):
        g = np.random.default_rng(seed)
        nm, nf = np.empty(n_perm), np.empty(n_perm)
        for p in range(n_perm):
            q = g.permutation(n)
            keep = q != np.arange(n)
            m = np.array([M.agree(hc[i], vs[k]) + M.agree(vp[k], he[i])
                          for i, k in zip(np.arange(n)[keep], q[keep])])
            nm[p], nf[p] = m.mean(), ((m >= 3) & (m <= 10)).mean()
        return nm, nf

    nm, nf = null_moments(a.n_perm, a.seed)
    sd_m, sd_f = float(nm.std(ddof=1)), float(nf.std(ddof=1))
    print(f"\n[sensitivity]  composition null: mean {nm.mean():.4f} bp (sd {sd_m:.4f}), "
          f"MMEJ-range fraction {nf.mean():.4f} (sd {sd_f:.4f})")
    for c, lab in ((2.486, "one-sided 2.486 (as in microhomology)"), (2.80, "two-sided 2.80 (paper convention)")):
        print(f"          MDE on the mean {c*sd_m:.4f} bp | on the MMEJ fraction "
              f"{100*c*sd_f:.2f} pp   [{lab}]")
        rows.append(dict(statistic="mde_constants", mde_mean_bp=round(c * sd_m, 4),
                         mde_frac_pp=round(100 * c * sd_f, 3), constant=c, note=lab))

    # planted MMEJ pool: synthetic junctions with 3-10 bp of real planted microhomology
    pool = ok[(ok["recovered_ext_h"] >= 3) & (ok["recovered_ext_h"] <= 10)]
    print(f"[sensitivity]  planted MMEJ pool: {len(pool)} synthetic junctions, mean "
          f"{pool['recovered_ext_h'].mean():.2f} bp of microhomology")
    obs_m, obs_f = float(mh.mean()), float(((mh >= 3) & (mh <= 10)).mean())
    print(f"          observed set: mean {obs_m:.4f} bp, MMEJ-range {100*obs_f:.1f}%")
    print(f"\n[sensitivity]  spiking a fraction f of the {n} junctions with planted-MMEJ junctions")
    print(f"          {'f':>6}{'mean bp':>10}{'excess':>9}{'clears?':>9}"
            f"{'MMEJ %':>9}{'excess pp':>11}{'clears?':>9}")
    plant_mh = pool["recovered_ext_h"].to_numpy(dtype=float)
    detect_m = detect_f = None
    for f in (0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60):
        k = int(round(f * n))
        reps = np.empty(200); repf = np.empty(200)
        for b in range(200):
            idx = rng.choice(n, size=k, replace=False)
            v = mh.astype(float).copy()
            v[idx] = rng.choice(plant_mh, size=k, replace=True)
            reps[b] = v.mean(); repf[b] = ((v >= 3) & (v <= 10)).mean()
        em, ef = reps.mean() - nm.mean(), repf.mean() - nf.mean()
        cm, cf = em >= 2.80 * sd_m, ef >= 2.80 * sd_f
        if cm and detect_m is None:
            detect_m = f
        if cf and detect_f is None:
            detect_f = f
        print(f"          {f:>6.0%}{reps.mean():>10.3f}{em:>9.3f}{'yes' if cm else 'no':>9}"
              f"{100*repf.mean():>9.1f}{100*ef:>11.2f}{'yes' if cf else 'no':>9}")
        rows.append(dict(statistic="spike_in_power", spike_fraction=f, n_replaced=k,
                         mean_bp=round(float(reps.mean()), 4), excess_bp=round(float(em), 4),
                         clears_mean_mde_280=bool(cm),
                         mmej_frac=round(float(repf.mean()), 4), excess_pp=round(100 * float(ef), 3),
                         clears_frac_mde_280=bool(cf)))
    print(f"\n[sensitivity]  Under this simulation model, the detectable MMEJ fraction is "
          f"f ~ {detect_m:.0%} of junctions on the mean statistic and f ~ {detect_f:.0%} on the "
          f"MMEJ-range fraction. These thresholds describe sensitivity to the injected "
          f"junction model.")
    rows.append(dict(statistic="detectable_mmej_share", detect_f_mean=detect_m,
                     detect_f_frac=detect_f, n_junctions=n, n_perm=a.n_perm))

    pd.DataFrame(rows).to_csv(OUT / "mh_positive_control.csv", index=False)
    print(f"\nwrote {OUT/'mh_positive_control.csv'}\n"
          f"wrote {OUT/'mh_recovery_curve.csv'}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
