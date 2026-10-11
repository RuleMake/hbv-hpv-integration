# GISTIC amplification-peak masking for spatial-null candidate pools.
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import pandas as pd
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WS = ROOT.parents[1]
DATA01 = WS / 'analysis_inputs/structural_variant_catalogue/donor_duplication_statistics/data'
sys.path.insert(0, str(WS / 'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
import gzip
import re
import spatial_null_core as F
GISTIC = DATA01 / 'gistic_peaks.txt.gz'
BIN = F.BIN_SIZE
MIN_POOL = 5
WINDOW = 100_000
def amp_peaks() -> pd.DataFrame:
    """GISTIC2 pan-cancer amplification peaks, hg19, from Wide Peak Limits."""
    rows = []
    with gzip.open(GISTIC, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        wide = header.index("Wide Peak Limits")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if not p or not p[0].startswith("Amplification Peak"):
                continue
            m = re.match(r"chr([0-9XY]+):(\d+)-(\d+)", p[wide].strip())
            if m:
                rows.append((f"chr{m.group(1)}", int(m.group(2)), int(m.group(3))))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end"]).drop_duplicates()
    print(f"[gistic]  {len(df)} pan-cancer amplification peaks "
          f"(median width {int((df.end - df.start).median()):,} bp)")
    return df

def masked_bins(peaks: pd.DataFrame, pad: int) -> set:
    out = set()
    for c, s, e in peaks.itertuples(index=False):
        lo = max(0, (s - pad) // BIN) * BIN
        hi = ((e + pad) // BIN) * BIN
        for b in range(lo, hi + BIN, BIN):
            out.add((c, b))
    return out

class MaskedNuller(F.Nuller):
    """Published 18-stratum pools with amplification-peak bins removed from the candidate set."""

    def __init__(self, mask: set) -> None:
        super().__init__()
        kept = {}
        removed = total = 0
        for (c, s), arr in self.bins.items():
            keep = np.array([(c, int(b)) not in mask for b in arr], dtype=bool)
            total += len(arr)
            removed += int((~keep).sum())
            if keep.any():
                kept[(c, s)] = arr[keep]
        self.bins = kept
        print(f"[mask]    {removed:,} / {total:,} candidate bins removed "
              f"({100 * removed / total:.1f}%); {len(self.bins)} non-empty pools remain")

def viable(frame: pd.DataFrame, nuller: F.Nuller) -> np.ndarray:
    """Coordinates whose matched pool still holds >= MIN_POOL candidates under this nuller.

    Coordinates that fail are DROPPED and counted, never given a fallback pool: falling back to
    the unmasked pool would reintroduce exactly the amplification-peak bins the arm exists to
    remove, and would silently turn the restricted arm back into the baseline.
    """
    ok = np.zeros(len(frame), dtype=bool)
    for i, (c, s, off, p0) in enumerate(zip(frame["chrom"], frame["stratum"],
                                            frame["within_bin_offset"], frame["pos0"])):
        bins = nuller.bins.get((str(c), int(s)))
        if bins is None or bins.size == 0:
            continue
        cand = np.unique(bins + int(off))
        ok[i] = (cand != int(p0)).sum() >= MIN_POOL
    return ok