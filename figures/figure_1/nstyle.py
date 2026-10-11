"""Typography, colour palette, millimetre layout and export checks.

Panel positions are specified in millimetres from the top-left corner
of the canvas.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb

# --------------------------------------------------------------------------- palette
# Use blue and orange to distinguish the viruses, with neutral colours for
# supporting elements.
C = {
    "HBV": "#17557F",
    "HPV": "#D1762F",
    "EBV": "#9BA19E",
    "HTLV-1": "#4A4A52",
    "ink": "#1A1F24",
    "mid": "#5E666D",
    "soft": "#8C9398",
    "rule": "#D6D4D0",
    "wash": "#F5F3EF",
    "wash2": "#EAE7E1",
    "viral": "#5C8F72",
    "seam": "#BE7D24",
    "mask": "#4A4A52",
    "white": "#FFFFFF",
}
MARK = {"HBV": "o", "HPV": "D", "EBV": "s", "HTLV-1": "^"}

MM = 25.4
FS = {"letter": 10.0, "title": 9.0, "label": 8.5, "tick": 8.0, "small": 8.0, "note": 8.0}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": FS["label"],
    "axes.labelsize": FS["label"],
    "axes.titlesize": FS["title"],
    "xtick.labelsize": FS["tick"],
    "ytick.labelsize": FS["tick"],
    "legend.fontsize": FS["small"],
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.0,
    "lines.solid_capstyle": "round",
    "patch.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.2,
    "ytick.major.size": 2.2,
    "xtick.minor.width": 0.45,
    "ytick.minor.width": 0.45,
    "xtick.minor.size": 1.3,
    "ytick.minor.size": 1.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "legend.handletextpad": 0.5,
    "legend.labelspacing": 0.45,
    "legend.borderpad": 0.0,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "text.color": C["ink"],
    "axes.labelcolor": C["ink"],
    "axes.edgecolor": C["ink"],
    "xtick.color": C["ink"],
    "ytick.color": C["ink"],
    "hatch.linewidth": 0.4,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.default": "regular",
})


def tint(color: str, amount: float = 0.5) -> tuple:
    """Blend towards white. amount=0 keeps the colour, amount=1 is white."""
    return tuple(np.asarray(to_rgb(color)) * (1 - amount) + amount)


def shade(color: str, amount: float = 0.25) -> tuple:
    """Blend towards the ink colour."""
    return tuple(np.asarray(to_rgb(color)) * (1 - amount) + np.asarray(to_rgb(C["ink"])) * amount)


# ---------------------------------------------------------------------------- canvas
class Canvas:
    """A figure whose panels are positioned in millimetres from the top-left."""

    def __init__(self, width_mm: float, height_mm: float):
        self.W = float(width_mm)
        self.H = float(height_mm)
        self.fig = plt.figure(figsize=(self.W / MM, self.H / MM))

    # -- geometry ----------------------------------------------------------------
    def _rect(self, x, y, w, h):
        return [x / self.W, 1.0 - (y + h) / self.H, w / self.W, h / self.H]

    def ax(self, x, y, w, h, **kw):
        """Axes with its top-left corner at (x, y) mm and size (w, h) mm."""
        a = self.fig.add_axes(self._rect(x, y, w, h), **kw)
        a.tick_params(pad=1.6)
        return a

    def canvas_ax(self, x, y, w, h):
        """A blank drawing surface in millimetre coordinates, for schematics."""
        a = self.fig.add_axes(self._rect(x, y, w, h))
        a.set_xlim(0, w)
        a.set_ylim(h, 0)  # y grows downwards, matching the layout convention
        a.set_axis_off()
        a.set_aspect("equal", adjustable="box")
        return a

    # -- labels ------------------------------------------------------------------
    def letter(self, ch, x, y, title=None, title_dx=4.6):
        """Panel letter at (x, y) mm, optionally followed by a short panel title."""
        self.fig.text(x / self.W, 1 - y / self.H, ch, fontsize=FS["letter"],
                      fontweight="bold", ha="left", va="top", color=C["ink"])
        if title:
            self.fig.text((x + title_dx) / self.W, 1 - (y + 0.15) / self.H, title,
                          fontsize=FS["title"], fontweight="bold", ha="left", va="top",
                          color=C["ink"])

    def text(self, s, x, y, **kw):
        kw.setdefault("fontsize", FS["small"])
        kw.setdefault("ha", "left")
        kw.setdefault("va", "top")
        kw.setdefault("color", C["ink"])
        return self.fig.text(x / self.W, 1 - y / self.H, s, **kw)

    def rule(self, x0, x1, y, color=None, lw=0.5):
        color = color or C["rule"]
        ln = mpl.lines.Line2D([x0 / self.W, x1 / self.W], [1 - y / self.H] * 2,
                              transform=self.fig.transFigure, color=color, lw=lw,
                              solid_capstyle="butt", zorder=0)
        self.fig.add_artist(ln)
        return ln


# ------------------------------------------------------------------------- utilities
def clean(ax, left=True, bottom=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    if not left:
        ax.tick_params(axis="y", length=0)
    return ax


def vref(ax, x=1.0, color=None, lw=0.65, ls="-"):
    ax.axvline(x, color=color or C["rule"], lw=lw, ls=ls, zorder=0.5)


def href(ax, y=0.0, color=None, lw=0.65, ls="-"):
    ax.axhline(y, color=color or C["rule"], lw=lw, ls=ls, zorder=0.5)


def thousands(v) -> str:
    return f"{int(round(v)):,}"


def fmt(v, nd=3) -> str:
    return f"{v:.{nd}f}"


def pstar(p) -> str:
    """Significance mark used only where the legend defines the family."""
    return "*" if (p is not None and p < 0.05) else ""


IT = {k: r"$\mathit{%s}$" % k for k in "Pznt"}


def ptext(p) -> str:
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return ""
    if p < 1e-4:
        return IT["P"] + " < 0.0001"
    if p < 0.001:
        return IT["P"] + f" = {p:.5f}".rstrip("0")
    return IT["P"] + f" = {p:.4f}".rstrip("0")


def dot_ci(ax, y, value, lo, hi, color, marker="o", size=4.0, filled=True, lw=1.0, zorder=3):
    """One horizontal estimate with its reference range."""
    face = color if filled else "white"
    ax.plot([lo, hi], [y, y], color=color, lw=lw, solid_capstyle="round", zorder=zorder)
    ax.plot([value], [y], marker=marker, ms=size, mfc=face, mec=color,
            mew=0.9, ls="none", zorder=zorder + 1)


def vdot_ci(ax, x, value, lo, hi, color, marker="o", size=4.0, filled=True, lw=1.0, zorder=3):
    face = color if filled else "white"
    ax.plot([x, x], [lo, hi], color=color, lw=lw, solid_capstyle="round", zorder=zorder)
    ax.plot([x], [value], marker=marker, ms=size, mfc=face, mec=color,
            mew=0.9, ls="none", zorder=zorder + 1)


def null_band(ax, lo, hi, y, half=0.26, colour=None):
    """A spatial-null reference range: light fill with end rules so it reads as an interval."""
    colour = colour or tint(C["rule"], 0.30)
    ax.add_patch(mpl.patches.Rectangle((lo, y - half), hi - lo, 2 * half, fc=colour,
                                       ec="none", zorder=1))
    for x in (lo, hi):
        ax.plot([x, x], [y - half, y + half], color=tint(C["mid"], 0.45), lw=0.6, zorder=2,
                solid_capstyle="butt")


def legend_items(entries, size=4.0):
    """entries: list of (label, color, marker, filled). Returns handles for ax.legend."""
    hs = []
    for label, color, marker, filled in entries:
        hs.append(mpl.lines.Line2D([], [], marker=marker, ls="none", ms=size,
                                   mfc=color if filled else "white", mec=color, mew=0.9,
                                   label=label))
    return hs


# ----------------------------------------------------------------------------- export
_LEAK = re.compile(
    r"\b(?:E\d{2,3}|G[12][′']?|rung|bc_uniqSV\w*|f_cluster\w*|paired_sample|donor_signature"
    r"|qc_pass|TODO|TBD|FIXME|draft_instruction|editorial_note|canonical|snapshot|stratum_eligible)\b"
    r"|[A-Za-z]:[/\\]|\.csv\b|\.parquet\b|working[/\\]|reference_path|_v\d{2}\b",
    re.I,
)


def _bbox_pairs(fig):
    """Overlapping pairs among visible, non-blank text artists."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    items = []
    for t in fig.findobj(mpl.text.Text):
        if not t.get_visible():
            continue
        s = t.get_text().strip()
        if not s:
            continue
        try:
            b = t.get_window_extent(rend)
        except Exception:
            continue
        if b.width <= 0 or b.height <= 0:
            continue
        items.append((s, b))
    hits = []
    for i in range(len(items)):
        si, bi = items[i]
        for j in range(i + 1, len(items)):
            sj, bj = items[j]
            ox = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
            oy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            if ox > 0.8 and oy > 0.8:  # ignore sub-pixel touching
                area = ox * oy
                small = min(bi.width * bi.height, bj.width * bj.height)
                if area / max(small, 1.0) > 0.12:
                    hits.append({"a": si, "b": sj, "overlap_frac": round(area / small, 3)})
    return hits, items


def sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export(canvas: Canvas, folder, stem, extra=None, allow_overlap=()):
    """Write svg/pdf/png plus a QA record. Canvas size is preserved; no tight crop."""
    fig = canvas.fig
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)

    # Measure the rendered typography before export.
    hits, items = _bbox_pairs(fig)
    allow = {tuple(sorted(p)) for p in allow_overlap}
    hits = [h for h in hits if tuple(sorted((h["a"], h["b"]))) not in allow]

    rend = fig.canvas.get_renderer()
    outside, sizes, leaks = [], [], []
    for s, b in items:
        if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > fig.bbox.width + 0.5 or b.y1 > fig.bbox.height + 0.5:
            outside.append(s)
        if _LEAK.search(s):
            leaks.append(s)
    for t in fig.findobj(mpl.text.Text):
        if t.get_visible() and t.get_text().strip():
            sizes.append(float(t.get_fontsize()))

    fig.savefig(folder / f"{stem}.svg")
    fig.savefig(folder / f"{stem}.pdf")
    fig.savefig(folder / f"{stem}.png", dpi=600)
    fig.savefig(folder / f"{stem}_preview.png", dpi=200)

    qa = {
        "figure": stem,
        "size_mm": [round(canvas.W, 2), round(canvas.H, 2)],
        "min_font_pt": round(min(sizes), 2) if sizes else None,
        "max_font_pt": round(max(sizes), 2) if sizes else None,
        "n_text_artists": len(items),
        "text_outside_canvas": outside,
        "text_overlaps": hits,
        "label_warnings": leaks,
        "svg_text_editable": True,
        "pdf_fonttype": 42,
        "png_dpi": 600,
        "files": {e: sha(folder / f"{stem}{e}") for e in (".svg", ".pdf", ".png", "_preview.png")},
    }
    if extra:
        qa.update(extra)
    (folder / f"{stem}_QA.json").write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    plt.close(fig)
    return qa
