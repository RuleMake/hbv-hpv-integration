"""Scores junction microhomology from published sequence reads against a composition-matched null.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import locus_class_propensity as E              # noqa: E402

SRC = Path(r"source_project")
XLS = SRC / "data" / "external" / "fan2023" / "Supplemental_Tables.xlsx"
OUT = E.PROJ / "results" / "junction_microhomology.csv"
OUT_J = E.PROJ / "results" / "junction_microhomology_per_junction.csv"
CACHE = E.PROJ / "results" / "cache" / "junction_sequence_cache"

FLANK = 80                       # hg38 window half-width around each human breakpoint
ANCHOR_FALLBACK = (30, 25, 20)   # read bases used to place a breakpoint by exact match
MAX_MH = 25                      # longest microhomology scored

HPV_ACC = {"HPV16": "K02718.1", "HPV18": "X05015.1", "HPV31": "J04353.1",
           "HPV33": "M12732.1", "HPV34": "X74476.1", "HPV45": "X74479.1",
           "HPV58": "D90400.1", "HPV66": "U31794.1"}

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def _get(url: str, key: str, parser) -> str:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{key}.txt"
    if f.exists():
        return f.read_text()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                out = parser(r.read().decode())
            f.write_text(out)
            return out
        except Exception:
            if attempt == 3:
                raise
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(url)


def hg38(chrom: str, start0: int, end0: int) -> str:
    url = (f"https://api.genome.ucsc.edu/getData/sequence?"
           f"genome=hg38;chrom={chrom};start={start0};end={end0}")
    return _get(url, f"hg38_{chrom}_{start0}_{end0}", lambda t: json.loads(t)["dna"].upper())


def hpv(vtype: str) -> str:
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
           f"db=nuccore&id={HPV_ACC[vtype]}&rettype=fasta&retmode=text")
    return _get(url, f"hpv_{vtype}", lambda t: "".join(t.split("\n")[1:]).upper())


def locate(anchor: str, ref: str, period: int | None = None):
    """Unique exact placement of `anchor` or its reverse complement. Circular refs pass `period`."""
    hits = []
    for a, strand in ((anchor, "+"), (rc(anchor), "-")):
        i, seen = ref.find(a), set()
        while i >= 0:
            seen.add(i % period if period else i)
            i = ref.find(a, i + 1)
        if len(seen) == 1:
            hits.append((seen.pop(), strand))
        elif len(seen) > 1:
            return None
    return hits[0] if len(hits) == 1 else None


def continuation(ref: str, idx: int, alen: int, strand: str, n: int) -> str:
    """Reference bases AFTER the anchor in the read's direction, read-oriented."""
    return ref[idx + alen: idx + alen + n] if strand == "+" else rc(ref[max(0, idx - n): idx])


def preceding(ref: str, idx: int, alen: int, strand: str, n: int) -> str:
    """Reference bases BEFORE the anchor in the read's direction, farthest-first.

    The minus-strand branch MUST offset by the anchor length: an anchor whose reverse complement
    occupies ref[idx : idx+alen] is read 3'->5' along the reference, so the bases preceding it in
    read direction are at ref[idx+alen ...], not at ref[idx ...]. Omitting that offset silently
    corrupts ext_v for every minus-strand anchor and biases microhomology downward -- it is the
    bug `selftest_bp` below exists to catch.
    """
    return (ref[max(0, idx - n): idx] if strand == "+"
            else rc(ref[idx + alen: idx + alen + n]))


def agree(a: str, b: str) -> int:
    k = 0
    for x, y in zip(a, b):
        if x != y or x == "N" or y == "N":
            break
        k += 1
    return k


def contexts(r) -> dict:
    """Place both breakpoints by alignment; return the four junction contexts, read-oriented."""
    left, right = r["seq"].split("|")
    hum_seg, vir_seg, hum_left = ((left, right, True) if r["int_type"] == "human-hpv"
                                  else (right, left, False))
    out = dict(sample=r["sample"], hpv_type=r["hpv_type"], chrom=r["chrom"],
               human_pos=int(r["human_pos"]), hpv_pos=int(r["hpv_pos"]),
               human_strand=r["human_strand"], int_type=r["int_type"], status=r["status"],
               gate="", h_cont="", h_end="", v_prec="", v_start="")

    ref_h = hg38(r["chrom"], max(0, int(r["human_pos"]) - 1 - FLANK), int(r["human_pos"]) + FLANK)
    loc_h = alen_h = None
    for a in ANCHOR_FALLBACK:
        loc_h = locate(hum_seg[-a:] if hum_left else hum_seg[:a], ref_h)
        if loc_h:
            alen_h = a
            break
    if not loc_h:
        out["gate"] = "human_anchor_unplaced"
        return out

    if r["hpv_type"] not in HPV_ACC:
        out["gate"] = "no_reference_for_type"
        return out
    ref_v0 = hpv(r["hpv_type"])
    ref_v, vlen = ref_v0 + ref_v0, len(ref_v0)
    loc_v = alen_v = None
    for a in ANCHOR_FALLBACK:
        loc_v = locate(vir_seg[:a] if hum_left else vir_seg[-a:], ref_v, period=vlen)
        if loc_v:
            alen_v = a
            break
    if not loc_v:
        out["gate"] = "viral_anchor_unplaced"
        return out

    ih, sh = loc_h
    iv, sv = loc_v
    if hum_left:
        out["h_cont"] = continuation(ref_h, ih, alen_h, sh, MAX_MH)
        out["v_start"] = vir_seg[:MAX_MH]
        out["v_prec"] = preceding(ref_v, iv, alen_v, sv, MAX_MH)[::-1]
        out["h_end"] = hum_seg[-MAX_MH:][::-1]
    else:
        out["h_cont"] = preceding(ref_h, ih, alen_h, sh, MAX_MH)[::-1]
        out["v_start"] = vir_seg[-MAX_MH:][::-1]
        out["v_prec"] = continuation(ref_v, iv, alen_v, sv, MAX_MH)
        out["h_end"] = hum_seg[:MAX_MH]
    # SELF-TEST. Walking backwards from the human anchor, away from the junction, must reproduce
    # the read's own human bases -- same parent, so agreement should run to the tested length.
    # Any orientation or offset error in continuation()/preceding() collapses this to ~0.
    if hum_left:
        back_ref = preceding(ref_h, ih, alen_h, sh, 40)[::-1]
        back_read = hum_seg[:-alen_h][::-1][:40]
    else:
        back_ref = continuation(ref_h, ih, alen_h, sh, 40)
        back_read = hum_seg[alen_h:alen_h + 40]
    out["selftest_bp"] = agree(back_ref, back_read)
    out["selftest_len"] = min(40, len(back_read))
    # same test on the VIRAL side, which is where ext_v comes from and where the circular
    # doubled reference is used; strain divergence from the type reference shows up here
    if hum_left:
        vb_ref = continuation(ref_v, iv, alen_v, sv, 40)
        vb_read = vir_seg[alen_v:alen_v + 40]
    else:
        vb_ref = preceding(ref_v, iv, alen_v, sv, 40)[::-1]
        vb_read = vir_seg[:-alen_v][::-1][:40]
    out["selftest_v_bp"] = agree(vb_ref, vb_read)
    out["selftest_v_len"] = min(40, len(vb_read))
    out.update(gate="ok", anchor_human_bp=alen_h, anchor_viral_bp=alen_v)
    return out


def load() -> pd.DataFrame:
    d = pd.read_excel(XLS, sheet_name="Table S2", header=2)
    d.columns = ["sample", "hpv_type", "hpv_pos", "chrom", "human_pos",
                 "human_strand", "int_type", "reads", "seq", "status"]
    d = d[d["sample"].notna()].reset_index(drop=True)
    d["seq"] = d["seq"].astype(str).str.replace(" ", "", regex=False)
    return d[d["seq"].str.count(r"\|") == 1].reset_index(drop=True)


def bucket(mh: np.ndarray) -> dict:
    return {"NHEJ-like 0-2 bp": float(((mh >= 0) & (mh <= 2)).mean()),
            "MMEJ-like 3-10 bp": float(((mh >= 3) & (mh <= 10)).mean()),
            "long >10 bp": float((mh > 10).mean())}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()

    d = load()
    if a.limit:
        d = d.head(a.limit)
    print(f"[fan]    {len(d)} junctions, {d.groupby(['chrom','human_pos']).ngroups} unique human "
          f"breakpoints, {d['hpv_type'].nunique()} papillomavirus types", flush=True)

    rows = []
    for i, (_, r) in enumerate(d.iterrows(), 1):
        rows.append(contexts(r))
        if i % 100 == 0:
            print(f"          placed {sum(1 for x in rows if x['gate']=='ok')}/{i}", flush=True)
    j = pd.DataFrame(rows)
    print("\n[check]   " + j["gate"].value_counts().to_string().replace("\n", "\n         "))

    ok = j[j["gate"] == "ok"].reset_index(drop=True)
    if ok.empty:
        print("nothing placed")
        return
    st = ok["selftest_bp"].to_numpy()
    full = (st >= ok["selftest_len"].to_numpy()).mean()
    print(f"[selftest] walking back into the human segment agrees for "
          f"median {np.median(st):.0f}/40 bp; {100*full:.1f}% reach the full tested length")
    if full < 0.90:
        print("           WARNING: human-side geometry self-test below 90% -- scores untrustworthy")
    stv = ok["selftest_v_bp"].to_numpy()
    fullv = (stv >= ok["selftest_v_len"].to_numpy()).mean()
    print(f"[selftest] walking back into the viral segment agrees for "
          f"median {np.median(stv):.0f}/40 bp; {100*fullv:.1f}% reach the full tested length "
          f"(shortfall here is HPV strain divergence from the type reference)")
    ext_h = np.array([agree(r.h_cont, r.v_start) for r in ok.itertuples()])
    ext_v = np.array([agree(r.v_prec, r.h_end) for r in ok.itertuples()])
    mh = ext_h + ext_v
    ok["ext_h"], ok["ext_v"], ok["microhomology_bp"] = ext_h, ext_v, mh
    ok.drop(columns=["h_cont", "h_end", "v_prec", "v_start"]).to_csv(OUT_J, index=False)

    print(f"\n[obs]    n = {len(mh)}   mean {mh.mean():.2f} bp   median {np.median(mh):.0f}   "
          f"(human side {ext_h.mean():.2f} + viral side {ext_v.mean():.2f})")
    for k, v in bucket(mh).items():
        print(f"         {k:20s} {100*v:5.1f}%")

    # ---- null: re-pair each human breakpoint with a viral breakpoint from the observed set
    rng = np.random.default_rng(a.seed)
    hc, he = ok["h_cont"].to_numpy(), ok["h_end"].to_numpy()
    vs, vp = ok["v_start"].to_numpy(), ok["v_prec"].to_numpy()
    n = len(ok)
    nm, nf = np.empty(a.n_perm), np.empty(a.n_perm)
    for p in range(a.n_perm):
        q = rng.permutation(n)
        keep = q != np.arange(n)
        m = np.array([agree(hc[i], vs[k]) + agree(vp[k], he[i])
                      for i, k in zip(np.arange(n)[keep], q[keep])])
        nm[p], nf[p] = m.mean(), ((m >= 3) & (m <= 10)).mean()

    res = []
    for name, o, nu in (("mean_microhomology_bp", float(mh.mean()), nm),
                        ("frac_MMEJ_3_10bp", float(((mh >= 3) & (mh <= 10)).mean()), nf)):
        sd = float(nu.std(ddof=1))
        e = o - float(nu.mean())
        res.append(dict(statistic=name, n_junctions=n, n_perm=a.n_perm,
                        observed=round(o, 4), null_mean=round(float(nu.mean()), 4),
                        excess=round(e, 4), z=round(e / sd, 2) if sd else np.nan,
                        mde_excess=round(2.80 * sd, 4),
                        clears_mde=bool(abs(e) >= 2.80 * sd)))
        print(f"[null]   {name:22s} obs {o:.4f}  null {nu.mean():.4f}  z {e/sd:+.2f}  "
              f"MDE {2.80*sd:.4f}  {'clears' if abs(e) >= 2.80*sd else 'FAILS'}")

    hist = pd.Series(mh).value_counts().sort_index()
    for k, v in hist.items():
        res.append(dict(statistic=f"count_mh_{int(k)}bp", observed=int(v),
                        n_junctions=n, n_perm=a.n_perm))
    for st in ok["status"].dropna().unique():
        m = mh[ok["status"].to_numpy() == st]
        res.append(dict(statistic=f"mean_mh_{st}", observed=round(float(m.mean()), 3),
                        n_junctions=len(m), n_perm=a.n_perm))
        print(f"[split]  {st:12s} n={len(m):>4}  mean {m.mean():.2f} bp")
    pd.DataFrame(res).to_csv(OUT, index=False)

    print(f"\nmicrohomology length histogram:\n{hist.to_string()}")
    print(f"\nwrote {OUT}\nwrote {OUT_J}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
