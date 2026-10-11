"""Parses CLCA amplicon intervals and summarizes the fields available for comparison.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

CLCA = Path(r"source_project\data\interim\clca")
IV = re.compile(r"(chr)?([0-9XY]+)\s*:\s*(\d+)\s*-\s*(\d+)")


def parse_intervals(s: str) -> list[tuple[str, int, int]]:
    """Amplicon `intervals` is a free-text list like 'chr1:179385001-203717000,chr3:129...'."""
    out = []
    for _p, c, a, b in IV.findall(str(s)):
        lo, hi = int(a), int(b)
        out.append((f"chr{c}", min(lo, hi), max(lo, hi)))
    return out


def main() -> None:
    amp = pd.read_csv(CLCA / "amplicons.csv", low_memory=False)
    ints = pd.read_csv(CLCA / "integrations.csv", low_memory=False)
    ints["chrom"] = "chr" + ints["chrom"].astype(str).str.replace("chr", "", regex=False)
    ints["pos"] = ints["pos"].astype(np.int64)

    print(f"[amp]   {len(amp):,} amplicon records, {amp['sample'].nunique()} samples")
    print(f"        amp_class: {amp['amp_class'].value_counts().head(6).to_dict()}")
    print(f"        avg_copy: median {amp['avg_copy'].median():.1f}, "
          f"range {amp['avg_copy'].min():.1f}-{amp['avg_copy'].max():.1f}, "
          f"missing {int(amp['avg_copy'].isna().sum())}")
    print(f"[ints]  {len(ints):,} integration records, {ints['sample'].nunique()} samples; "
          f"read support present: split {int(ints['split_reads'].notna().sum())}, "
          f"pair {int(ints['pair_reads'].notna().sum())}")

    shared = set(amp["sample"]) & set(ints["sample"])
    print(f"[join]  samples with BOTH amplicons and integrations: {len(shared)}")

    # expand amplicon intervals
    rows = []
    for _, r in amp.iterrows():
        for c, lo, hi in parse_intervals(r["intervals"]):
            rows.append((r["sample"], c, lo, hi, r["avg_copy"], r["amp_class"]))
    seg = pd.DataFrame(rows, columns=["sample", "chrom", "start", "end",
                                      "avg_copy", "amp_class"])
    print(f"[seg]   {len(seg):,} amplicon intervals parsed "
          f"({seg['sample'].nunique()} samples); median width "
          f"{(seg['end'] - seg['start']).median()/1e6:.2f} Mb")

    # integration inside an amplicon interval of the SAME sample
    hit_cn, hit_cls = [], []
    by = {(s, c): g[["start", "end", "avg_copy", "amp_class"]].to_numpy(object)
          for (s, c), g in seg.groupby(["sample", "chrom"])}
    for s, c, p in zip(ints["sample"], ints["chrom"], ints["pos"]):
        arr = by.get((s, c))
        cn, cl = np.nan, None
        if arr is not None:
            for st, en, v, k in arr:
                if st <= p <= en:
                    if not np.isfinite(cn) or (np.isfinite(v) and v > cn):
                        cn, cl = float(v), k
        hit_cn.append(cn); hit_cls.append(cl)
    ints["amp_cn"] = hit_cn
    ints["amp_class"] = hit_cls

    inside = ints["amp_cn"].notna()
    print(f"\n[hit]   {int(inside.sum()):,}/{len(ints):,} integrations fall inside an amplicon "
          f"of their own sample ({100*inside.mean():.1f}%)")
    if inside.sum():
        d = ints[inside]
        print(f"        across {d['sample'].nunique()} samples; "
              f"amp_cn median {d['amp_cn'].median():.1f}, "
              f"IQR {d['amp_cn'].quantile(.25):.1f}-{d['amp_cn'].quantile(.75):.1f}, "
              f"max {d['amp_cn'].max():.1f}")
        d = d.assign(support=d["split_reads"].fillna(0) + d["pair_reads"].fillna(0))
        print(f"        read support median {d['support'].median():.0f}, "
              f"IQR {d['support'].quantile(.25):.0f}-{d['support'].quantile(.75):.0f}")
        per = d.groupby("sample").size()
        print(f"        integrations per sample: median {per.median():.0f}, "
              f"max {per.max()}, samples with >=2 distinct amp_cn: "
              f"{int((d.groupby('sample')['amp_cn'].nunique() >= 2).sum())}")
        print("\n        Within-sample dose-response estimation requires samples carrying integrations")
        print("        in amplicons of DIFFERENT copy number. That count is the number above.")


if __name__ == "__main__":
    main()
