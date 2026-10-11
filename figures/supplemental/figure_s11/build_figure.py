"""Plot HPV clonality results in Supplemental Figure S11.

Read the two bundled count and test-summary tables, then write the figure
and displayed counts to this directory.
"""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

OUT = Path(__file__).resolve().parent
SRC = OUT / "clonality_rebinning.csv"
DON = OUT / "direction_band_decomposition_donors.csv"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def rows(p):
    with p.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))

source_files = [SRC, DON]
source_hashes = {p.name: sha(p) for p in source_files}
saved = rows(SRC)
donors = rows(DON)
presence = [r for r in saved if r["rule"] == "rule_invariant_presence_contrast"]
permutation = next(r for r in saved if r["rule"] == "presence_contrast_within_donor_permutation")
subset_test = next(r for r in saved if r["rule"] == "presence_contrast_donors_in_both_arms")
near_ids = {r["donor"] for r in donors if r["band"] != "NONE"}
far_ids = {r["donor"] for r in donors if r["band"] == "NONE"}
shared = near_ids & far_ids
assert len(near_ids | far_ids) == 10 and len(shared) == 4

displayed = []
for scope in ("All 10 tumours", "4 tumours contributing to both groups"):
    for nearby in (True, False):
        selected = [r for r in donors if (r["band"] != "NONE") == nearby
                    and (scope == "All 10 tumours" or r["donor"] in shared)]
        n = sum(int(r["n"]) for r in selected)
        sub = sum(int(r["n_subclonal"]) for r in selected)
        displayed.append(dict(scope=scope, nearby_sv=nearby, n_subclonal=sub,
                              n_junctions=n, pct_subclonal=100 * sub / n,
                              n_tumours=len({r["donor"] for r in selected})))
assert [(r["n_subclonal"], r["n_junctions"]) for r in displayed] == [(25,99),(12,19),(24,61),(10,17)]
assert [(int(r["n_subclonal"]),int(r["n"])) for r in presence] == [(25,99),(12,19)]
assert int(permutation["n"]) == 118 and int(subset_test["n"]) == 78
assert float(presence[0]["p_omnibus"]) == 0.00225021
assert float(permutation["p_omnibus"]) == 0.137543
assert float(subset_test["p_omnibus"]) == 0.175779

with (OUT / "displayed_presence_counts.csv").open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(displayed[0])); w.writeheader();w.writerows(displayed)
tests = [
    {"test":"pooled Fisher", "scope":"all 118 junction records from 10 tumours", "p":float(presence[0]["p_omnibus"])},
    {"test":"within-donor permutation", "scope":"all 118 junction records; proximity labels shuffled within each of 10 donors; only 4 donors have both labels", "n_permutations":20000,"seed":7,"p":float(permutation["p_omnibus"])},
    {"test":"pooled Fisher", "scope":"78 junction records from the 4 tumours contributing to both groups", "p":float(subset_test["p_omnibus"])}]
(OUT / "saved_tests.json").write_text(json.dumps(tests,indent=2)+"\n",encoding="utf-8")

style = {"font.family":"Arial", "font.size":8, "axes.labelsize":8,
         "axes.titlesize":9,"xtick.labelsize":8,"ytick.labelsize":8,
         "svg.fonttype":"none", "axes.linewidth":0.7,
         "figure.facecolor":"white","axes.facecolor":"white",
         "savefig.facecolor":"white", "text.color":"#20262D",
         "axes.labelcolor":"#20262D","xtick.color":"#20262D","ytick.color":"#20262D"}
with mpl.rc_context(style):
    fig, axes=plt.subplots(1,2,figsize=(178/25.4,105/25.4),sharey=True)
    fig.subplots_adjust(left=.09,right=.985,bottom=.28,top=.80,wspace=.22)
    colours=["#0072B2","#D55E00"]
    for ix,ax in enumerate(axes):
        data=displayed[ix*2:ix*2+2]
        bars=ax.bar([0,1],[r["pct_subclonal"] for r in data],width=.56,
                    color=colours,edgecolor="#20262D",linewidth=.65,zorder=3)
        bars[1].set_hatch("///")
        for x,r in enumerate(data):
            ax.text(x,r["pct_subclonal"]+3,f'{r["pct_subclonal"]:.1f}%\n({r["n_subclonal"]}/{r["n_junctions"]})',
                    ha="center",va="bottom",fontsize=8,linespacing=1.3)
        ax.set(ylim=(0,100),xlim=(-.60,1.60),yticks=[0,25,50,75,100])
        ax.set_xticks([0,1],["SV within 100 kb","No SV within 100 kb"])
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100,decimals=0))
        ax.grid(axis="y",color="#DFE3E8",linewidth=.55,zorder=0)
        ax.tick_params(axis="x",length=0,pad=7)
        ax.tick_params(axis="y",length=3)
        ax.spines[["top","right"]].set_visible(False)
        ax.text(-.20,1.28,chr(65+ix),transform=ax.transAxes,fontsize=12,fontweight="bold",va="top")
        header="All 10 tumours" if ix==0 else "4 tumours contributing to both groups"
        ax.text(0,1.26,header,transform=ax.transAxes,fontsize=9,fontweight="bold",va="top")
        ax.text(0,1.12,"118 junction records" if ix==0 else "78 junction records",transform=ax.transAxes,fontsize=8,va="top",color="#4C5661")
    axes[0].set_ylabel("Subclonal junctions (%)",labelpad=7)
    axes[0].text(0,-.31,"Pooled Fisher  P = 0.0023\nWithin-donor permutation  P = 0.138",transform=axes[0].transAxes,fontsize=7.8,va="top",linespacing=1.55)
    axes[1].text(0,-.31,"Descriptive subset of panel A",transform=axes[1].transAxes,fontsize=7.8,va="top",color="#4C5661")
    fig.savefig(OUT/"Supplemental_Fig_S11.svg")
    fig.savefig(OUT/"Supplemental_Fig_S11.png",dpi=600)
    fig.savefig(OUT/"Supplemental_Fig_S11_preview.png",dpi=180)
    plt.close(fig)

assert all(sha(OUT/name)==digest for name,digest in source_hashes.items())
print(json.dumps({"counts":displayed,"source_hashes_unchanged":True},indent=2))
