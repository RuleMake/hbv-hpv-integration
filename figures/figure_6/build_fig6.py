"""Figure 6: plot preserved observations and prepared OLS summaries."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import nstyle as ns
from nstyle import C

HERE = Path(__file__).resolve().parent
DATA = HERE/'plotted_data'
arm = pd.read_csv(DATA/'Figure_6A_plotted_values.csv', float_precision='round_trip')
amp = pd.read_csv(DATA/'Figure_6B_plotted_values.csv', float_precision='round_trip')
points = pd.read_csv(DATA/'Figure_6C_plotted_values.csv', float_precision='round_trip')
fits = pd.read_csv(DATA/'Figure_6C_fits_plotted_values.csv', float_precision='round_trip')

cv = ns.Canvas(180, 183)
cv.text('Figure 6', 178, 2, ha='right', fontsize=9, fontweight='bold')
geometry = []

cv.letter('A', 2, 2, 'HPV duplication enrichment')
cv.text('Fusion transcript', 53.5, 14, ha='center', va='center', fontsize=8.5)
a = cv.ax(21, 25, 65, 46)
ns.clean(a); ns.href(a, 1, color=C['soft'], lw=.65)
for i, (key, label) in enumerate([('PRODUCTIVE', 'Detected'), ('SILENT', 'Not detected')]):
    r = arm[arm.subset == key].iloc[0]
    err = a.errorbar(i, r.fold, yerr=[[r.fold-r.fold_ci_lo], [r.fold_ci_hi-r.fold]],
                    fmt='none', ecolor=C['HPV'], elinewidth=.95, capsize=3, capthick=.75, zorder=2)
    err.lines[2][0].set_gid(f'A_range_{i}')
    a.plot(i, r.fold, marker='D', ms=5.1, mfc=C['HPV'], mec=C['HPV'], mew=.8,
           ls='none', zorder=3, gid=f'A_point_{i}')
    geometry.append(dict(panel='A', subset=key, x=i, estimate=float(r.fold), lower=float(r.fold_ci_lo),
                         upper=float(r.fold_ci_hi), axis_min=.9, axis_max=1.85))
a.set(xlim=(-.65, 1.65), ylim=(.9, 1.85), xticks=[0, 1], yticks=[1, 1.2, 1.4, 1.6, 1.8])
a.set_yticklabels(['1.0', '1.2', '1.4', '1.6', '1.8'])
a.set_xticklabels(['Detected\n73 clusters', 'Not detected\n274 clusters'])
a.set_ylabel('Observed / matched null', fontsize=8.5)

cv.letter('B', 94, 2, 'HBV amplicon overlap')
cv.text('625 clusters · 193 samples', 143.5, 11, fontsize=8, color=C['mid'], ha='center', va='center')
key = cv.ax(108, 16, 69, 6); key.set_axis_off()
key.legend(handles=[mpl.patches.Patch(fc=ns.tint(C['HBV'], .22), ec=C['HBV'], lw=.5, label='Observed'),
                    mpl.patches.Patch(fc=ns.tint(C['mid'], .73), ec=ns.tint(C['mid'], .3), lw=.5, label='Null mean')],
           loc='center', ncol=2, fontsize=8, handlelength=1, columnspacing=1.5)
b = cv.ax(111, 25, 65, 46); ns.clean(b)
for i, r in amp.iterrows():
    for status, dx, value in [('observed', -.18, r.observed), ('expected', .18, r.expected)]:
        detected = status == 'observed'
        color = ns.tint(C['HBV'], .22) if detected else ns.tint(C['mid'], .73)
        edge = C['HBV'] if detected else ns.tint(C['mid'], .3)
        bars = b.bar(i+dx, value, width=.29, color=color, ec=edge, lw=.45, zorder=2)
        bars[0].set_gid(f'B_bar_{i}_{status}')
        b.annotate(f'{value:g}', (i+dx, value), xytext=(0, 4), textcoords='offset points',
                   ha='center', va='bottom', fontsize=8, color=C['HBV'] if detected else C['mid'])
        geometry.append(dict(panel='B', group=int(i), status=status, estimate=float(value), axis_min=0, axis_max=50))
b.set(xlim=(-.6, 1.6), ylim=(0, 50), xticks=[0, 1], xticklabels=['All amplicons', 'Circular class'],
      yticks=[0, 10, 20, 30, 40, 50])
b.set_ylabel('Integration clusters', fontsize=8.5)

cv.letter('C', 2, 91, 'HBV junction support and amplicon copy number')
for i, (group, title, xlabel) in enumerate([
        ('amplicon', 'Inside amplicons', 'Mean log$_2$ amplicon\ncopy number'),
        ('placebo', 'Outside amplicons', 'Log$_2$ sample-median\namplicon copy number')]):
    x = 21 if i == 0 else 111
    width = 65
    g = points[points.arm == group]
    r = fits[fits.arm == group].iloc[0]
    cv.text(f'{title} · {len(g)} samples', x+width/2, 103, ha='center', va='center', fontsize=8.5, fontweight='bold')
    ax = cv.ax(x, 111, width, 47)
    ns.clean(ax)
    xs = np.linspace(g.log_cn.min(), g.log_cn.max(), 50)
    ax.plot(xs, r.intercept+r.slope*xs, color=ns.shade(C['HBV'], .3), lw=1.1,
            zorder=2, gid=f'C_fit_{group}')
    scatter = ax.scatter(g.log_cn, g.log_support, s=g.marker_area, facecolors=ns.tint(C['HBV'], .55),
                         edgecolors=C['HBV'], linewidths=.75, zorder=3)
    scatter.set_gid(f'C_points_{group}')
    ax.set(xlim=(1.2, 4.6), ylim=(2.6, 10), xticks=[1.5, 2.5, 3.5, 4.5], yticks=[3, 5, 7, 9])
    ax.set_xlabel(xlabel, fontsize=8.5, labelpad=3)
    if i == 0:
        ax.set_ylabel('Mean log$_2$ junction read support', fontsize=8.5)
    for _, p in g.iterrows():
        geometry.append(dict(panel='C', arm=group, sample=p['sample'], x=float(p.log_cn), y=float(p.log_support),
                             area=float(p.marker_area), count=int(p.n_integrations), xlim=[1.2,4.6], ylim=[2.6,10]))

size_key = cv.ax(34, 177, 132, 3); size_key.set_axis_off()
size_key.set(xlim=(0,132), ylim=(-1,1))
size_key.text(64, 0, 'Integrations per sample', ha='right', va='center', fontsize=8, color=C['mid'])
nmax = int(points.n_integrations.max())
for i, n in enumerate([1, 5, 10]):
    size_key.scatter([70+i*16], [0], s=8+46*n/nmax, facecolors=ns.tint(C['HBV'], .55),
                     edgecolors=C['HBV'], linewidths=.75, clip_on=False)
    size_key.text(75+i*16, 0, str(n), va='center', fontsize=8, color=C['mid'])

qa = ns.export(cv, HERE/'rendered', 'Yan_Fig6', extra=dict(
    plotted_geometry=geometry, panel_A_intervals=2, panel_A_interval_type='95% spatial reference ranges',
    panel_B_values=4, panel_B_baseline=0, panel_B_intervals_drawn=False,
    panel_B_circular_upper_reference_bound='no finite upper bound; retained as NaN in the original source',
    panel_C_points={k:int((points.arm==k).sum()) for k in ['amplicon','placebo']},
    panel_C_marker_area=f'8 + 46 * n_integrations / {nmax}; shared by both scatter panels and size legend',
    panel_C_weighting='Unweighted OLS on source-sample means; marker area does not weight the fit',
    panel_C_uncertainty='Slope 95% t confidence intervals are in the figure legend and prepared source table',
    panel_C_fit_table='plotted_data/Figure_6C_fits_plotted_values.csv', new_statistical_tests=False))
print({k:qa[k] for k in ['min_font_pt','text_overlaps','text_outside_canvas','label_warnings']})
