"""Assess genomic-covariate balance for HTLV-1 spatial comparisons.

The same balance estimator is applied to the HTLV-1, HBV, HPV and EBV coordinate
sets. The coordinate loader adds the HTLV-1 groups and redirects output to this
comparison's result file."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F           # noqa: E402
import covariate_balance_check as BALANCE        # noqa: E402

PROJ = Path(r"source_project")
OUT = PROJ / "results"
COORDS = OUT / "htlv_coordinates.csv"

_original_loader = F.load_coordinates


def load_with_htlv() -> pd.DataFrame:
    """Load the same six HTLV-1 groups used for class-enrichment comparisons.
    
    The three recorded clinical-state groups are supplemented by the pooled group
    and two within-study subsets from PMID 35654946. Groups are derived locally to
    avoid running the class-atlas module during import."""
    co = _original_loader()[["virus_class", "chrom", "pos0"]].copy()
    h = pd.read_csv(COORDS)

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

    htlv = h.rename(columns={"arm": "virus_class"})[["virus_class", "chrom", "pos0"]]
    print(f"[htlv]    {len(htlv):,} HTLV-1 coordinates "
          f"{htlv['virus_class'].value_counts().to_dict()}")
    return pd.concat([co, htlv], ignore_index=True)


class _Redirect:
    """Redirect balance-summary output to the current comparator result file."""

    def __truediv__(self, name: str) -> Path:
        return OUT / "integrase_positive_balance.csv"


def main() -> None:
    F.load_coordinates = load_with_htlv
    BALANCE.F.load_coordinates = load_with_htlv
    BALANCE.OUT = _Redirect()
    if not any(arg.startswith("--viruses") for arg in sys.argv):
        sys.argv += ["--viruses", "HBV,HPV,EBV,HTLV_all,HTLV_ATL,HTLV_nonmalignant,"
                     "HTLV_unknown,HTLV_ATL_s35654946,HTLV_nonmalig_s35654946"]
    BALANCE.main()


if __name__ == "__main__":
    main()
