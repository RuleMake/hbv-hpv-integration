"""Builds junction-support and copy-number comparisons used by the amplicon analyses.
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

SRC_ROOT = Path(r"source_project")
sys.path.insert(0, str(SRC_ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from scipy.stats import spearmanr, mannwhitneyu  # noqa: E402

from src.common.config import (  # noqa: E402
    chrom_norm, LARGE_TD_MIN, LARGE_TD_MAX,
)

OUT = Path(r"source_project\results")
CLCA = SRC_ROOT / "data" / "interim" / "clca"
WIN = 1_000_000


def parse_intervals(s):
    out = []
    if not isinstance(s, str):
        return out
    for tok in s.split(","):
        m = re.match(r"^(chr[\dXYM]+):(\d+)-(\d+)$", tok.strip())
        if m:
            out.append((chrom_norm(m.group(1)), int(m.group(2)), int(m.group(3))))
    return out


def main() -> None:
    ints = pd.read_csv(OUT / "ascertainment_integration_annot.csv")
    print("=" * 78)
    print("EVENT-ORDER INFERENCE FROM READ SUPPORT")
    print("=" * 78)
    print(f"integrations: {len(ints)} in {ints['sample'].nunique()} samples")

    # within-sample normalisation removes depth / purity / sample clonality
    med = ints.groupby("sample")["support"].transform("median")
    ints["support_norm"] = ints["support"] / med.replace(0, np.nan)
    ints["log_support"] = np.log2(ints["support"].clip(lower=1))

    # ---------------------------------------------------------------- ANALYSIS A
    print("\n" + "-" * 78)
    print("A. Does junction support scale with local copy number?")
    print("   integration-first -> slope > 0 ; integration-last -> slope ~ 0")
    print("-" * 78)
    sub = ints.dropna(subset=["amp_avg_copy"]).copy()
    print(f"  integrations inside an amplicon with known avg_copy: n={len(sub)} "
          f"in {sub['sample'].nunique()} samples")
    if len(sub) >= 20:
        rho, p = spearmanr(sub["amp_avg_copy"], sub["support"])
        print(f"  raw           Spearman rho={rho:.3f}  p={p:.4g}")
        rho2, p2 = spearmanr(sub["amp_avg_copy"], sub["support_norm"])
        print(f"  within-sample Spearman rho={rho2:.3f}  p={p2:.4g}")

        sub["log_cn"] = np.log2(sub["amp_avg_copy"].clip(lower=0.1))
        try:
            m = smf.ols("log_support ~ log_cn", data=sub).fit(
                cov_type="cluster", cov_kwds={"groups": sub["sample"]})
            b = m.params["log_cn"]
            ci = m.conf_int().loc["log_cn"]
            print(f"\n  log2(support) ~ log2(copy number), sample-clustered SE:")
            print(f"    slope = {b:.3f}   95% CI [{ci[0]:.3f}, {ci[1]:.3f}]   p = {m.pvalues['log_cn']:.4g}")
            print(f"    slope 0 = integration AFTER amplification (junction on one allele)")
            print(f"    slope 1 = integration BEFORE amplification (fully co-amplified)")
            if m.pvalues["log_cn"] < 0.05 and b > 0:
                print("    -> SUPPORTS INTEGRATION-FIRST for amplified loci")
            else:
                print("    -> no significant scaling; cannot infer order")
        except Exception as e:
            print(f"  [model failed: {e}]")
    else:
        print("  too few amplicon-resident integrations to model")

    # ---------------------------------------------------------------- ANALYSIS B
    print("\n" + "-" * 78)
    print("B. Are TD-enclosed integrations elevated in support?")
    print("-" * 78)
    sv = pd.read_csv(CLCA / "sv_api.csv")
    sv["c1"] = sv["chra"].map(chrom_norm)
    sv["c2"] = sv["chrb"].map(chrom_norm)
    sv = sv[sv["c1"] == sv["c2"]].copy()
    sv["span"] = (sv["posb"] - sv["posa"]).abs()
    sv["sv_class"] = ["DUP" if (a == "-" and b == "+") else
                      "DEL" if (a == "+" and b == "-") else
                      "INV" if a == b else "OTHER"
                      for a, b in zip(sv["strand_a"], sv["strand_b"])]
    tds = sv[(sv.sv_class == "DUP") & (sv.span >= LARGE_TD_MIN) & (sv.span <= LARGE_TD_MAX)]
    td_by = {}
    for smp, g in tds.groupby("caseid"):
        d = {}
        for _, t in g.iterrows():
            lo, hi = sorted([int(t["posa"]), int(t["posb"])])
            d.setdefault(t["c1"], []).append((lo, hi))
        td_by[smp] = d

    enc = []
    for _, r in ints.iterrows():
        hit = False
        for lo, hi in td_by.get(r["sample"], {}).get(r["chrom"], []):
            if lo < r["pos"] < hi:
                hit = True
                break
        enc.append(hit)
    ints["td_enclosed"] = enc
    a = ints.loc[ints.td_enclosed, "support_norm"].dropna()
    b = ints.loc[~ints.td_enclosed, "support_norm"].dropna()
    print(f"  TD-enclosed    : n={len(a):4d}  median norm support={a.median():.3f}")
    print(f"  not enclosed   : n={len(b):4d}  median norm support={b.median():.3f}")
    if len(a) >= 10:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        print(f"  Mann-Whitney p={p:.4g}  -> "
              f"{'enclosed integrations carry MORE support (consistent with co-duplication)' if a.median() > b.median() and p < 0.05 else 'no significant difference'}")
    else:
        print("  too few enclosed integrations to test")

    # ---------------------------------------------------------------- ANALYSIS C
    print("\n" + "-" * 78)
    print("C. Does CN-residualised clonality predict proximal rearrangement?")
    print("-" * 78)
    # residualise support on the CN handles we have, so the remainder indexes clonality
    ints["in_amp"] = ints["in_amplicon"].astype(int)
    base = ints.dropna(subset=["support_norm"]).copy()
    try:
        r_m = smf.ols("np.log2(support_norm.clip(lower=0.01)) ~ in_amp + td_enclosed",
                      data=base).fit()
        base["clonality_resid"] = r_m.resid
        print(f"  residualised support on amplicon membership + TD enclosure "
              f"(R2={r_m.rsquared:.3f}, n={len(base)})")
    except Exception as e:
        print(f"  [residualisation failed: {e}]")
        base["clonality_resid"] = np.nan

    svL = sv[(sv.span >= LARGE_TD_MIN) & (sv.span <= LARGE_TD_MAX)]
    bymap = {}
    for smp, g in svL.groupby("caseid"):
        m = {}
        for _, t in g.iterrows():
            m.setdefault(t["c1"], []).extend([int(t["posa"]), int(t["posb"])])
        bymap[smp] = {c: np.array(sorted(v)) for c, v in m.items()}
    prox = []
    for _, r in base.iterrows():
        arr = bymap.get(r["sample"], {}).get(r["chrom"])
        if arr is None or arr.size == 0:
            prox.append(0)
        else:
            prox.append(int(np.searchsorted(arr, r["pos"] + WIN, "right")
                            - np.searchsorted(arr, r["pos"] - WIN, "left")))
    base["n_proximal_sv"] = prox
    base["any_proximal"] = (base["n_proximal_sv"] > 0).astype(int)

    ok = base.dropna(subset=["clonality_resid"])
    print(f"  integrations with >=1 proximal large SV: "
          f"{int(ok['any_proximal'].sum())}/{len(ok)}")
    if len(ok) > 50:
        rho, p = spearmanr(ok["clonality_resid"], ok["n_proximal_sv"])
        print(f"  clonality residual vs proximal SV count: Spearman rho={rho:.3f}, p={p:.4g}")
        try:
            lg = smf.logit("any_proximal ~ clonality_resid", data=ok).fit(disp=0)
            b_ = lg.params["clonality_resid"]
            print(f"  logistic any_proximal ~ clonality: beta={b_:.3f}, "
                  f"OR={np.exp(b_):.3f}, p={lg.pvalues['clonality_resid']:.4g}")
            if lg.pvalues["clonality_resid"] < 0.05:
                direction = ("clonal/early integrations carry MORE proximal SV "
                             "(consistent with integration-first)") if b_ > 0 else \
                            ("subclonal/late integrations carry more proximal SV "
                             "(consistent with rearrangement-first)")
                print(f"    -> {direction}")
            else:
                print("    -> NULL: clonality does not predict proximal rearrangement")
        except Exception as e:
            print(f"  [logit failed: {e}]")

    base.to_csv(OUT / "clonality_timing.csv", index=False)
    print(f"\nWrote {OUT / 'clonality_timing.csv'}")


if __name__ == "__main__":
    main()
