"""Compare amplicon and non-amplicon slopes with a common sample-level estimator.

Both comparisons use sample means. A patient sign-flip test compares the observed
score with its null on the same scale; leave-one-sample-out estimates describe
sensitivity to individual observations. Reference calculations are checked before
results are written.

This source retains normal-based intervals. Figure 6 preparation uses the saved
estimates and standard errors to calculate the t-based intervals displayed in
the manuscript."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

RES = Path(r"source_project\results")
SEED, N_BOOT = 7, 20_000
K = 2.80


def ols(m: pd.DataFrame) -> tuple[float, float, float, float]:
    r = sm.OLS(m.log_support, sm.add_constant(m.log_cn)).fit()
    b, se = float(r.params["log_cn"]), float(r.bse["log_cn"])
    return b, b - 1.96 * se, b + 1.96 * se, float(r.pvalues["log_cn"])


def signflip(m: pd.DataFrame, rng) -> tuple[float, float, float, float, float]:
    """Return (p on the score scale, p as the shipped code computes it, max |t_null|,
    the observed value on the score scale, the largest single sample's share of sum(g^2), and the
    observed value on the standard-error scale that the shipped code wrongly compares against it)."""
    y, x = m.log_support.to_numpy(float), m.log_cn.to_numpy(float)
    xc, yc = x - x.mean(), y - y.mean()
    denom = float((xc * xc).sum())
    b = float((xc * y).sum() / denom)
    se = float(np.sqrt(((yc - b * xc) ** 2).sum() / (len(y) - 2) / denom))
    g = xc * yc / denom
    sd = float(np.sqrt((g * g).sum()))
    t = (rng.choice([-1.0, 1.0], size=(N_BOOT, len(y))) @ g) / sd
    return (float((np.abs(t) >= abs(b / sd)).mean()),
            float((np.abs(t) >= abs(b / se)).mean()),
            float(np.abs(t).max()), b / sd, float((g * g).max() / (g * g).sum()), b / se)


def _share(m: pd.DataFrame) -> np.ndarray:
    """Each sample's share of sum(g^2) -- its leverage over the sign-flip null's support."""
    y, x = m.log_support.to_numpy(float), m.log_cn.to_numpy(float)
    xc, yc = x - x.mean(), y - y.mean()
    g = xc * yc / float((xc * xc).sum())
    return (g * g) / float((g * g).sum())


def main() -> None:
    d = pd.read_csv(RES / "clonality_timing.csv")
    amp = d[d.in_amplicon == True].dropna(subset=["amp_avg_copy"]).copy()  # noqa: E712
    amp["log_cn"] = np.log2(amp.amp_avg_copy.clip(lower=0.1))
    non = d[d.in_amplicon == False].copy()                                 # noqa: E712
    non["cn"] = non["sample"].map(amp.groupby("sample").amp_avg_copy.median())
    non = non.dropna(subset=["cn"])
    non["log_cn"] = np.log2(non.cn.clip(lower=0.1))

    def means(t):
        m = t.groupby("sample")[["log_support", "log_cn"]].mean().reset_index()
        m["n"] = t.groupby("sample").size().to_numpy()
        return m
    ma, mp = means(amp), means(non)

    # ---- reproduction rung, four values, before anything new is computed
    b, lo, hi, p = ols(ma)
    rng = np.random.default_rng(SEED)
    pf = np.random.default_rng(SEED)
    xc = ma.log_cn.to_numpy(float) - ma.log_cn.mean()
    yc = ma.log_support.to_numpy(float) - ma.log_support.mean()
    _ = pf.choice([-1.0, 1.0], size=(N_BOOT, len(ma)))          # the shipped code draws these first
    p_perm = float((np.abs([float((xc * pf.permutation(yc)).sum() / (xc * xc).sum())
                            for _ in range(N_BOOT)]) >= abs(b)).mean())
    r2 = sm.OLS(non.log_support, sm.add_constant(non.log_cn)).fit(
        cov_type="cluster", cov_kwds={"groups": non["sample"]})
    b2, se2 = float(r2.params["log_cn"]), float(r2.bse["log_cn"])
    want = {"slope": (round(b, 4), 1.0275), "ci_lo": (round(lo, 4), 0.3764),
            "ci_hi": (round(hi, 4), 1.6787), "p": (round(p, 5), 0.00417),
            "perm_p": (round(p_perm, 5), 0.0042), "placebo_slope": (round(b2, 4), 0.0309),
            "placebo_lo": (round(b2 - 1.96 * se2, 4), -0.2204),
            "placebo_hi": (round(b2 + 1.96 * se2, 4), 0.2823),
            "placebo_p": (round(float(r2.pvalues["log_cn"]), 5), 0.80938)}
    for k, (got, exp) in want.items():
        if abs(got - exp) > 1e-4:
            raise SystemExit(f"RUNG FAILED  {k}: got {got}, timing_between_sample_check says {exp}")
    print(f"Reference check: {len(want)} stored values reproduced exactly")

    p_score, p_ascoded, tmax, t_score, share, t_se = signflip(ma, np.random.default_rng(SEED))
    if p_ascoded != 0.0:
        raise SystemExit("the shipped bootstrap no longer returns 0.0; re-examine before claiming a defect")
    print(f"Reference implementation: original sign-flip p {p_ascoded}, null support tops out at {tmax:.3f}")

    # ---- leave-one-sample-out, the robustness the bootstrap cannot speak to
    fits = [ols(ma[ma["sample"] != s]) for s in ma["sample"]]
    loo = np.array([f[0] for f in fits])
    loo_p = np.array([f[3] for f in fits])
    pd.DataFrame(dict(sample=ma["sample"], n_integrations=ma.n, slope_without=loo.round(4),
                      p_without=loo_p.round(5),
                      share_of_score=_share(ma))).to_csv(
        RES / "coamplification_loo.csv", index=False)
    bp, lop, hip, pp = ols(mp)
    p2_score = signflip(mp, np.random.default_rng(SEED))[0]

    rows = [
        dict(arm="amplicon", estimator="OLS on sample means", unit="sample mean", n_units=len(ma),
             n_rows=len(amp), slope=round(b, 4), ci_lo=round(lo, 4), ci_hi=round(hi, 4),
             p=round(p, 5), note="the panel's estimate"),
        dict(arm="placebo", estimator="OLS on sample means", unit="sample mean", n_units=len(mp),
             n_rows=len(non), slope=round(bp, 4), ci_lo=round(lop, 4), ci_hi=round(hip, 4),
             p=round(pp, 5), note="like-for-like control; same estimator as the arm above"),
        dict(arm="placebo", estimator="OLS per row, cluster SE", unit="integration", n_units=30,
             n_rows=len(non), slope=round(b2, 4), ci_lo=round(b2 - 1.96 * se2, 4),
             ci_hi=round(b2 + 1.96 * se2, 4), p=round(float(r2.pvalues["log_cn"]), 5),
             note="as published; NOT like-for-like, interval borrows within-sample replication"),
        dict(arm="amplicon", estimator="exact pairing permutation", unit="sample mean",
             n_units=len(ma), n_rows=len(amp), slope=round(b, 4), ci_lo=np.nan, ci_hi=np.nan,
             p=round(p_perm, 5), note=f"{N_BOOT} draws"),
        dict(arm="amplicon", estimator="sign-flip bootstrap, score scale", unit="sample mean",
             n_units=len(ma), n_rows=len(amp), slope=round(b, 4), ci_lo=np.nan, ci_hi=np.nan,
             p=round(p_score, 5),
             note=f"observed on the score scale {t_score:.3f}; the correct comparison"),
        dict(arm="amplicon", estimator="sign-flip bootstrap, as shipped", unit="sample mean",
             n_units=len(ma), n_rows=len(amp), slope=round(b, 4), ci_lo=np.nan, ci_hi=np.nan,
             p=p_ascoded,
             note=f"DEFECT: observed studentised at {t_se:.3f} on the standard-error scale and "
                  f"compared against a score-scale null whose support ends at {tmax:.3f}; "
                  f"p is forced by that mismatch, not measured"),
        dict(arm="placebo", estimator="sign-flip bootstrap, score scale", unit="sample mean",
             n_units=len(mp), n_rows=len(non), slope=round(bp, 4), ci_lo=np.nan, ci_hi=np.nan,
             p=round(p2_score, 5), note="same estimator as the arm above"),
        dict(arm="amplicon", estimator="leave-one-sample-out", unit="sample mean", n_units=len(ma),
             n_rows=len(amp), slope=round(float(loo.mean()), 4), ci_lo=round(float(loo.min()), 4),
             ci_hi=round(float(loo.max()), 4), p=round(float(loo_p.max()), 5),
             note="ci columns are the LOO slope range; p is the worst LOO p"),
        dict(arm="amplicon", estimator="leverage", unit="sample mean", n_units=len(ma),
             n_rows=len(amp), slope=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=np.nan,
             note=f"largest single sample carries {share:.3f} of sum(g^2); "
                  f"sign-flip null support ends at |t| = {tmax:.3f}"),
    ]
    out = pd.DataFrame(rows)
    out.to_csv(RES / "coamplification_like_for_like.csv", index=False)

    pts = pd.concat([ma.assign(arm="amplicon"), mp.assign(arm="placebo")], ignore_index=True)
    pts = pts[["arm", "sample", "log_cn", "log_support", "n"]].rename(columns={"n": "n_integrations"})
    pts.to_csv(RES / "coamplification_points.csv", index=False)

    pd.set_option("display.width", 220)
    print()
    print(out.to_string(index=False))
    print(f"\nwrote coamplification_like_for_like.csv ({len(out)} rows) and "
          f"coamplification_points.csv ({len(pts)} points)")


if __name__ == "__main__":
    main()
