"""Figure 1. Study design and integration-coordinate clustering."""
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, PathPatch, Rectangle
from matplotlib.path import Path as MPath

import nstyle as ns
from nstyle import C, FS

HERE = Path(__file__).resolve().parent
DATA = HERE / "source_data" / "main"
PLOT = HERE / "plotted_data"
OUT = HERE / "rendered"

cv = ns.Canvas(180, 229)
fig = cv.fig
cv.text("Figure 1", 178, 1.5, fontsize=9, fontweight="bold", ha="right")

SEG = dict(fc=ns.tint(C["mid"], 0.80), ec=C["mid"], lw=0.55)
H = 4.6


def host_block(ax, x, y, w, lab):
    ax.add_patch(Rectangle((x, y), w, H, **SEG, zorder=2))
    ax.text(x + w / 2, y + H / 2, lab, ha="center", va="center", fontsize=FS["small"],
            style="italic", color=C["ink"], zorder=3)


def viral_block(ax, x, y, w, seams=True):
    ax.add_patch(Rectangle((x, y), w, H, fc=ns.tint(C["viral"], 0.55), ec=C["viral"],
                           lw=0.55, hatch="////", zorder=2))
    if seams:
        for jx in (x, x + w):
            ax.plot([jx, jx], [y - 1.2, y + H + 1.2], color=C["seam"], lw=1.3,
                    solid_capstyle="butt", zorder=4)


# ============================================================== A  two routes, one outcome
cv.letter("A", 1.0, 1.5, "Integration and host rearrangement")
ax = cv.canvas_ax(1.0, 7.0, 178, 60)
ax.add_patch(Rectangle((0.6, 0.4), 176.8, 59.2, fc=C["wash"], ec="none", zorder=0))

# Two schematic histories, with separate label and connector lanes.
ax.text(4.0, 5.0, "Integration followed by rearrangement", fontsize=FS["label"],
        fontweight="bold", va="center")
for x, lab in [(4, "a"), (14, "b"), (33, "c")]:
    host_block(ax, x, 13.2, 10, lab)
viral_block(ax, 24, 13.2, 9)

ax.text(4.0, 31.0, "Viral capture during break repair", fontsize=FS["label"],
        fontweight="bold", va="center")
viral_block(ax, 21.5, 37.5, 9, seams=False)
ax.add_patch(FancyArrowPatch((26.0, 42.5), (26.0, 45.7), arrowstyle="-|>", mutation_scale=5,
                             color=C["viral"], lw=0.8, zorder=4))
yb = 46.5
for x, lab in [(4, "a"), (14, "b"), (28, "b"), (38, "c")]:
    host_block(ax, x, yb, 10, lab)
for jx in (24.4, 27.6):
    ax.add_patch(PathPatch(MPath([(jx, yb), (jx + 0.9, yb + 1.5), (jx - 0.9, yb + 3.1),
                                  (jx, yb + H)], [1, 2, 2, 2]), fc="none", ec=C["seam"],
                           lw=1.1, zorder=4))

for y0, rad, step, label_y in ((15.5, 0.14, "local restructuring", 12.0),
                              (48.8, -0.14, "repair", 49.0)):
    ax.add_patch(FancyArrowPatch((55, y0), (91, 28.2), arrowstyle="-|>", mutation_scale=6,
                                 color=C["mid"], lw=0.9, zorder=3,
                                 connectionstyle=f"arc3,rad={rad}"))
    ax.text(73, label_y, step, fontsize=FS["note"], color=C["mid"],
            ha="center", va="center")

ax.add_patch(Rectangle((96, 13.0), 78, 28.0, fc="white", ec=C["rule"], lw=0.7, zorder=1))
ax.text(135, 18.3, "Observed structure", ha="center", va="center", fontsize=FS["label"],
        fontweight="bold")
for x, lab in [(106, "a"), (118, "b"), (141, "b"), (153, "c")]:
    host_block(ax, x, 25.3, 12, lab)
viral_block(ax, 130, 25.3, 11)
ax.text(135, 36.8, "Viral DNA between host junctions", ha="center",
        va="center", fontsize=FS["note"], color=C["mid"])

h_viral = Rectangle((0, 0), 1, 1, fc=ns.tint(C["viral"], 0.55), ec=C["viral"], lw=0.55,
                    hatch="////")
h_seam = mpl.lines.Line2D([], [], color=C["seam"], lw=1.3)
h_brk = mpl.lines.Line2D([], [], color=C["seam"], ls="none", mew=1.1,
                         marker=MPath([(0, -1), (.4, -.35), (-.4, .35), (0, 1)]),
                         ms=5, mfc="none")
ax.legend([h_viral, h_seam, h_brk], ["viral DNA", "host junction", "unrepaired break"],
          loc="lower right", bbox_to_anchor=(0.985, 0.035), ncol=3, fontsize=FS["note"],
          handlelength=1.3, columnspacing=1.2, borderaxespad=0.0)

# ======================================================================= B  three analyses
cv.letter("B", 1.0, 71.5, "Analysis framework")
bx = cv.canvas_ax(1.0, 78.0, 178, 25)
for x0, title in ((4, "Counting unit"), (63, "Position in a duplication"),
                  (122, "Class composition")):
    bx.text(x0, 2.0, title, fontsize=FS["label"], fontweight="bold", va="center")
for sx in (57, 116):
    bx.plot([sx, sx], [0.5, 23.5], color=C["rule"], lw=0.6)

# B1 clustered records collapse to loci
y = 12.0
bx.plot([4, 52], [y, y], color=C["mid"], lw=0.8)
groups = [(9, 2), (22, 1), (33, 2), (42, 1)]
for cx, k in groups:
    for t in ([cx] if k == 1 else np.linspace(cx - 1.4, cx + 1.4, k)):
        bx.plot([t, t], [y - 2.4, y + 2.4], color=C["mid"], lw=0.9)
bx.text(52, 6.4, "Coordinates", fontsize=FS["note"], color=C["mid"],
        va="center", ha="right")
y2 = 19.4
bx.plot([4, 52], [y2, y2], color=C["mid"], lw=0.8)
for cx, k in groups:
    bx.plot([cx], [y2], marker="o", ms=3.4, mfc=C["mid"], mec=C["mid"], ls="none")
    bx.add_patch(FancyArrowPatch((cx, y + 3.0), (cx, y2 - 2.0), arrowstyle="-|>",
                                 mutation_scale=4, color=C["soft"], lw=0.6))
bx.text(52, 22.6, "Genomic clusters", fontsize=FS["note"], color=C["mid"],
        va="center", ha="right")

# B2 boundary versus interior of a duplicated interval
dy = 13.0
bx.add_patch(Rectangle((63, dy - 3.0), 48, 6.0, fc=ns.tint(C["mid"], 0.88), ec="none"))
for sx in (63, 103):
    bx.add_patch(Rectangle((sx, dy - 3.0), 8, 6.0, fc=ns.tint(C["mid"], 0.62), ec="none"))
bx.plot([63, 111], [dy, dy], color=C["mid"], lw=0.8)
for ix in (67, 87):
    bx.plot([ix, ix], [dy - 4.8, dy + 4.8], color=C["seam"], lw=1.3)
bx.text(77, 6.2, "Integration sites", fontsize=FS["note"], color=C["seam"],
        ha="center", fontweight="bold")
bx.text(67, 20.6, "boundary", fontsize=FS["note"], color=C["mid"], ha="center")
bx.text(87, 20.6, "interior", fontsize=FS["note"], color=C["mid"], ha="center")
bx.text(107, 20.6, "boundary", fontsize=FS["note"], color=C["mid"], ha="center")
for xa, xb in ((63, 71), (103, 111)):
    bx.annotate("", (xa, 17.6), (xb, 17.6),
                arrowprops=dict(arrowstyle="<->", color=C["soft"], lw=0.6, shrinkA=0, shrinkB=0))

# B3 four rearrangement classes on a signed composition axis
GL = C["mid"]
# DEL : a segment removed
bx.plot([122, 126], [10.4, 10.4], color=GL, lw=1.1)
bx.plot([126, 130], [10.4, 10.4], color=C["rule"], lw=1.1, ls=(0, (1.6, 1.2)))
bx.plot([130, 134], [10.4, 10.4], color=GL, lw=1.1)
# DUP : a segment present twice
bx.plot([136, 148], [10.4, 10.4], color=GL, lw=1.1)
bx.plot([140, 144], [8.6, 8.6], color=GL, lw=1.1)
# INV : opposing orientations
bx.annotate("", (150, 11.0), (156, 11.0),
            arrowprops=dict(arrowstyle="->", color=GL, lw=0.9, shrinkA=0, shrinkB=0))
bx.annotate("", (162, 9.6), (156, 9.6),
            arrowprops=dict(arrowstyle="->", color=GL, lw=0.9, shrinkA=0, shrinkB=0))
# TRA : two chromosomes exchanged
bx.plot([164, 174], [9.0, 11.6], color=GL, lw=1.0)
bx.plot([164, 174], [11.6, 9.0], color=GL, lw=1.0)
for lab, cx in (("DEL", 128), ("DUP", 142), ("INV", 156), ("TRA", 169)):
    bx.text(cx, 5.8, lab, fontsize=FS["note"], fontweight="bold", ha="center", va="center")

bx.annotate("", (122, 17.0), (176, 17.0),
            arrowprops=dict(arrowstyle="<->", color=C["mid"], lw=0.7, shrinkA=0, shrinkB=0))
bx.text(122, 21.4, "DEL + DUP", fontsize=FS["note"], color=C["mid"], ha="left")
bx.text(176, 21.4, "TRA + INV", fontsize=FS["note"], color=C["mid"], ha="right")

# =========================================================== C  spacing within sample-locus
cv.letter("C", 1.0, 111.5, "Coordinate spacing within source samples and loci")
sp = pd.read_csv(DATA / "Figure_1C.csv")
summ = pd.read_csv(PLOT / "Figure_1C_plotted_values.csv")

axc = cv.ax(15, 127, 68, 40)
ns.clean(axc)
for virus in ("HBV", "HPV"):
    for kind, ls in (("observed", "-"), ("independent_within_locus", (0, (3.0, 1.6)))):
        s = np.sort(np.maximum(sp[(sp.virus == virus) & (sp.kind == kind)].spacing_bp.values, 1))
        axc.step(s, np.arange(1, len(s) + 1) / len(s), where="post", color=C[virus], lw=1.1,
                 ls=ls, zorder=3)
axc.set_xscale("log")
axc.set_xlim(1, 1.6e5)
axc.set_ylim(-0.02, 1.02)
axc.set_xticks([1, 10, 100, 1e3, 1e4, 1e5])
axc.set_xticklabels(["1", "10", "100", "1 kb", "10 kb", "100 kb"])
axc.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axc.set_xlabel("Adjacent spacing (bp)")
axc.set_ylabel("Cumulative fraction")
ns.href(axc, 0.5, ls=(0, (2.2, 1.8)))
for virus, fx, fy, ha in (("HBV", 0.40, 0.085, "right"), ("HPV", 9.0, -0.215, "left")):
    m = np.median(sp[(sp.virus == virus) & (sp.kind == "observed")].spacing_bp)
    axc.plot([m], [0.5], marker=ns.MARK[virus], ms=3.8, mfc=C[virus], mec=C[virus], ls="none",
             zorder=5)
    axc.annotate(f"{m:,.1f} bp" if m % 1 else f"{m:,.0f} bp", (m, 0.5), (m * fx, 0.5 + fy), fontsize=FS["note"],
                 color=C[virus], ha=ha, va="center",
                 arrowprops=dict(arrowstyle="-", color=C[virus], lw=0.5, shrinkA=0, shrinkB=1.5))
legend_ax = cv.ax(15, 117, 68, 8)
legend_ax.set_axis_off()
leg = legend_ax.legend(handles=[
    mpl.lines.Line2D([], [], color=C["HBV"], lw=1.1, label="HBV"),
    mpl.lines.Line2D([], [], color=C["HPV"], lw=1.1, label="HPV"),
    mpl.lines.Line2D([], [], color=C["mid"], lw=1.1, label="observed"),
    mpl.lines.Line2D([], [], color=C["mid"], lw=1.1, ls=(0, (3.0, 1.6)), label="independent"),
], loc="center left", ncol=2, fontsize=FS["note"],
    handlelength=1.6, columnspacing=1.0, borderaxespad=0)

# C right : median and spread, observed against independent placement
cv.text("Median and 10th–90th percentiles", 108, 122, fontsize=8, color=C['mid'])
axd = cv.ax(122, 127, 55, 40)
ns.clean(axd, left=False)
rows, ticks = [], []
for virus in ("HBV", "HPV"):
    for kind, filled in (("observed", True), ("independent_within_locus", False)):
        r = summ[(summ.virus == virus) & (summ.placement == ('Observed' if filled else 'Independent'))].iloc[0]
        rows.append((virus, filled, r.displayed_median_bp, r.displayed_p10_bp, r.displayed_p90_bp))
        ticks.append("observed" if filled else "independent")
for i, (virus, filled, med, lo, hi) in enumerate(rows):
    yy = len(rows) - 1 - i
    axd.plot([lo, hi], [yy, yy], color=ns.tint(C[virus], 0.0 if filled else 0.45), lw=1.1,
             solid_capstyle="butt", zorder=3)
    axd.vlines([lo, hi], yy - 0.09, yy + 0.09,
               color=ns.tint(C[virus], 0.0 if filled else 0.45), lw=0.75, zorder=3)
    axd.plot([med], [yy], marker=ns.MARK[virus], ms=4.2, mfc=C[virus] if filled else "white",
             mec=C[virus], mew=0.9, ls="none", zorder=4)
    axd.text(med, yy + 0.22, f"{med:,.1f} bp" if med % 1 else f"{med:,.0f} bp", fontsize=FS["note"],
             color=C[virus] if filled else C["mid"], va="bottom", ha="center")
axd.set_xscale("log")
axd.set_xlim(1, 1.6e5)
axd.set_ylim(-0.7, 3.7)
axd.set_yticks(range(len(rows))[::-1])
axd.set_yticklabels(ticks, fontsize=FS["note"])
axd.set_xticks([1, 100, 1e4])
axd.set_xticklabels(["1", "100", "10 kb"])
axd.set_xlabel("Spacing (bp)")
for i, virus in enumerate(("HBV", "HPV")):
    axd.text(-0.36, 3 - i * 2 - 0.5, virus, transform=axd.get_yaxis_transform(),
             fontsize=FS["note"], fontweight="bold", color=C[virus], ha="center",
             va="center", rotation=90)

# ============================================================ D  genome-wide distributions
cv.letter("D", 1.0, 179.0, "Genome-wide distribution in 5-Mb bins")
dd = pd.read_csv(PLOT / "Figure_1D_plotted_values.csv", dtype={'chrom':str})
peaks = pd.read_csv(PLOT / "Figure_1D_query_exclusions.csv", dtype={'chrom':str})
order = [str(i) for i in range(1, 23)] + ["X"]
lens = dd.groupby("chrom").chrom_length.first()
GAP = 12e6
offset, cuts, run = {}, [], 0.0
for ch in order:
    offset[ch] = run
    run += lens[ch] + GAP
    cuts.append(run - GAP / 2)
TOT = run - GAP
dd["gx"] = dd.chrom.map(offset) + (dd.bin_start + dd.bin_end) / 2
peaks["gx0"] = peaks.chrom.astype(str).map(offset) + peaks.start
peaks["gx1"] = peaks.chrom.astype(str).map(offset) + peaks.end


def round_tick(v):
    step = 10 ** int(np.floor(np.log10(v)))
    for m in (1, 2, 2.5, 5, 10):
        if m * step <= v:
            best = m * step
    return int(best)


top = 187.0
for name, col, hgt in (("SV breakends", C["mid"], 10.0), ("HBV", C["HBV"], 8.0),
                       ("HPV", C["HPV"], 8.0)):
    a = cv.ax(32, top, 146, hgt)
    ns.clean(a, left=False, bottom=False)
    g = dd[dd.track == name].sort_values("gx")
    a.bar(g.gx, g["count"], width=(g.bin_end-g.bin_start), color=col, lw=0,
          align="center", zorder=2)
    a.set_xlim(-GAP / 2, TOT + GAP / 2)
    a.set_ylim(0, g["count"].max() * 1.14)
    a.set_yticks([])
    a.tick_params(axis="x", length=0, labelbottom=False)
    for c in cuts[:-1]:
        a.axvline(c, color=C["rule"], lw=0.4, zorder=0)
    a.text(-0.012, 0.76, name, transform=a.transAxes, ha="right", va="center",
           fontsize=FS["note"], color=col, fontweight="bold")
    a.text(-0.012, 0.24, f"0–{int(g['count'].max()):,}", transform=a.transAxes,
           ha="right", va="center", fontsize=8, color=C["soft"])
    top += hgt + 1.3

mask_ax = cv.ax(32, top + 1.2, 146, 2.2)
mask_ax.set_axis_off()
mask_ax.set_xlim(-GAP / 2, TOT + GAP / 2)
mask_ax.set_ylim(0, 1)
for r in peaks.itertuples(index=False):
    mask_ax.add_patch(Rectangle((r.gx0, .28), r.gx1-r.gx0, .44,
                               fc=C['mid'],ec='none'))
mask_ax.text(-.012,.5,'Viral query exclusions',transform=mask_ax.transAxes,
             fontsize=8,color=C['mid'],ha='right',va='center')

lab_ax = cv.ax(32, top + 5.2, 146, 5.0)
lab_ax.set_axis_off()
lab_ax.set_xlim(-GAP / 2, TOT + GAP / 2)
lab_ax.set_ylim(0, 1)
chrom_labels = [lab_ax.text(offset[ch] + lens[ch] / 2, .80, ch,
                ha="center", va="center", fontsize=8, color=C["mid"]) for ch in order]

# One baseline; minimally adjust crowded labels horizontally, within their chromosomes.
# Bounded isotonic projection separates labels while keeping their centres within their chromosomes.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
centres = np.array([lab_ax.transData.transform(t.get_position())[0] for t in chrom_labels])
widths = np.array([t.get_window_extent(renderer).width for t in chrom_labels])
gap_px = .55 / 25.4 * fig.dpi
required_steps = (widths[:-1]+widths[1:])/2 + gap_px
cumulative = np.r_[0., np.cumsum(required_steps)]
target = centres - cumulative
margin_px = .15 / 25.4 * fig.dpi
lower=np.array([lab_ax.transData.transform((offset[ch],.8))[0]+margin_px for ch in order])-cumulative
upper=np.array([lab_ax.transData.transform((offset[ch]+lens[ch],.8))[0]-margin_px for ch in order])-cumulative
pools = []
for i, value in enumerate(target):
    pools.append([i, i+1, float(value), 1, lower[i], upper[i]])
    while len(pools)>1 and np.clip(pools[-2][2]/pools[-2][3],pools[-2][4],pools[-2][5]) > np.clip(pools[-1][2]/pools[-1][3],pools[-1][4],pools[-1][5]):
        right=pools.pop()
        left=pools.pop()
        low=max(left[4],right[4])
        high=min(left[5],right[5])
        assert low<=high, 'Label-spacing constraints cannot fit within chromosome spans'
        pools.append([left[0],right[1],left[2]+right[2],left[3]+right[3],low,high])
fitted = np.empty(len(order))
for start,end,total,weight,low,high in pools:
    fitted[start:end] = np.clip(total/weight,low,high)
adjusted = fitted+cumulative
chrom_label_positions = []
for ch,t,x0,x1 in zip(order,chrom_labels,centres,adjusted):
    px,py = lab_ax.transData.transform(t.get_position())
    new_x = lab_ax.transData.inverted().transform((x1,py))[0]
    assert offset[ch] <= new_x <= offset[ch]+lens[ch]
    t.set_x(new_x)
    chrom_label_positions.append(dict(chromosome=ch, displacement_mm=float((x1-x0)/fig.dpi*25.4),
                                      within_chromosome=True, baseline=.8))
fig.canvas.draw()
chrom_boxes=[t.get_window_extent(fig.canvas.get_renderer()) for t in chrom_labels]
chrom_label_min_gap_mm=min((b.x0-a.x1)/fig.dpi*25.4 for a,b in zip(chrom_boxes,chrom_boxes[1:]))
assert chrom_label_min_gap_mm >= .549

for aligned_axis in (axc, axd):
    aligned_axis.xaxis.set_label_coords(.5, -.18)
    aligned_axis.xaxis.label.set_verticalalignment("baseline")
fig.canvas.draw()
legend_clearance_mm = (leg.get_window_extent(fig.canvas.get_renderer()).y0 - axc.bbox.y1) / fig.dpi * 25.4
assert legend_clearance_mm >= 1.0, legend_clearance_mm
c_axis_baseline_difference_mm = abs(axc.bbox.y0-axd.bbox.y0) / fig.dpi * 25.4
c_xlabel_baseline_difference_mm = abs(axc.xaxis.label.get_transform().transform(axc.xaxis.label.get_position())[1]-axd.xaxis.label.get_transform().transform(axd.xaxis.label.get_position())[1]) / fig.dpi * 25.4
c_plot_gutter_mm = (axd.bbox.x0-axc.bbox.x1) / fig.dpi * 25.4
assert c_axis_baseline_difference_mm < .01
assert c_xlabel_baseline_difference_mm < .01
assert c_plot_gutter_mm >= 38.9

qa = ns.export(cv, OUT, "Yan_Fig1", extra={
    "panel_sources": ["source_data/main/Figure_1C.csv",
                      "plotted_data/Figure_1C_plotted_values.csv", "plotted_data/Figure_1D_plotted_values.csv",
                      "plotted_data/Figure_1D_query_exclusions.csv"],
    "cohort_caption_source": "plotted_data/Figure_1B_legend_values.csv",
    "schematic_panels": "Two possible histories and three analysis components; cohort values retained in caption only",
    "spacing_rows_plotted": int(len(sp)),
    "zero_spacings_clipped_to_1bp": int((sp.spacing_bp == 0).sum()),
    "bins_per_track": int((dd.track == "HBV").sum()),
    "query_exclusion_display": "Effective 50-kb filter bins with 100-kb peak padding, in a separate aligned annotation track explicitly labelled Viral query exclusions",
    "legend_to_plot_clearance_mm": float(legend_clearance_mm),
    "c_axis_baseline_difference_mm": c_axis_baseline_difference_mm,
    "c_xlabel_baseline_difference_mm": c_xlabel_baseline_difference_mm,
    "c_plot_gutter_mm": c_plot_gutter_mm,
    "a_height_mm": 60,
    "chromosome_label_positions": chrom_label_positions,
    "chromosome_label_min_gap_mm": chrom_label_min_gap_mm,
    "terminal_bar_widths": "Actual bin_end minus bin_start",
    "new_statistical_tests": False,
})
print(qa["min_font_pt"], "| overlaps:", qa["text_overlaps"], "| leaks:",
      qa["label_warnings"], "| outside:", qa["text_outside_canvas"])
