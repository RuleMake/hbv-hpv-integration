"""Plot categorical enrichment profiles, reference ranges and calibration summaries."""
from pathlib import Path
import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
import numpy as np
import pandas as pd
import nstyle as ns
from nstyle import C

ROOT=Path(__file__).resolve().parent
P=ROOT/'plotted_data'
OUT=ROOT/'rendered'
cv=ns.Canvas(180,201)
fig=cv.fig
cv.text('Figure 2',178,1.8,fontsize=9,fontweight='bold',ha='right')
geom=[]
refline=dict(color=C['mid'],lw=.7,ls=(0,(2.5,2)),zorder=1)

# A. Three categorical profiles share the same physical scale and complete ranges.
cv.letter('A',2,1.8,'Duplication enrichment')
a=pd.read_csv(P/'Figure_2A_plotted_values.csv')
axesA=[]
profile_paths=[]
# Categorical spacing: a small extra gutter before the wider final label.
# Labels, estimates and ranges remain centred on the same category positions.
category_x=np.array([0.,1.,2.,3.,4.09])
for vi,v in enumerate(['HPV','HBV','EBV']):
    left=19+54*vi
    cv.text(v,left+24.5,9.4,ha='center',fontsize=9,fontweight='bold')
    ax=cv.ax(left,23,49,33)
    axesA.append(ax)
    ns.clean(ax,left=vi==0)
    ax.set(xlim=(-.4,4.4),ylim=(.8,1.6),xticks=category_x,yticks=[.8,1,1.2,1.4,1.6])
    ax.set_xticklabels(['Indep.','Cluster','Indep.','Cluster','Refinedᵃ'],fontsize=8)
    if vi==0:
        ax.set_ylabel('Observed / expected overlap')
    else:
        ax.tick_params(labelleft=False,length=0,axis='y')
    ax.axhline(1,color=C['soft'],lw=.7,zorder=0)
    for yy in [1.2,1.4,1.6]:
        ax.axhline(yy,color=ns.tint(C['rule'],.35),lw=.4,zorder=0)
    ax.axvline(1.5,color=ns.tint(C['rule'],.15),lw=.55,zorder=0)
    rows=a[a.virus==v].sort_values('scheme')
    values=rows.observed_expected.to_numpy()
    ax.plot(category_x[:2],values[:2],color=C[v],lw=1.4,zorder=3)
    ax.plot(category_x[1:3],values[1:3],color=C[v],lw=1.,ls=(0,(2.2,1.8)),zorder=3)
    ax.plot(category_x[2:],values[2:],color=C[v],lw=1.4,zorder=3)
    profile_paths.append(dict(virus=v,x=category_x.tolist(),scheme_order=list(range(1,6)),values=values.tolist(),
        solid_segments=[[0,1],[2,3,4]],dashed_segment=[1,2],smoothed=False))
    for r in rows.itertuples(index=False):
        x=category_x[r.scheme-1]
        interval=ax.errorbar(x,r.observed_expected,
            yerr=[[r.observed_expected-r.lower],[r.upper-r.observed_expected]],
            fmt='none',ecolor=ns.tint(C[v],.23),elinewidth=.8,capsize=2.2,capthick=.7,zorder=2)
        interval.lines[2][0].set_gid(f'A_range_{v}_{r.scheme}')
        point,=ax.plot(x,r.observed_expected,marker=ns.MARK[v],ms=4.8,mec='white',mew=.45,
                mfc=C[v],ls='none',zorder=4)
        point.set_gid(f'A_point_{v}_{r.scheme}')
        geom.append(dict(panel='A',virus=v,scheme=r.scheme,estimate=r.observed_expected,
                         lower=r.lower,upper=r.upper,axis_min=.8,axis_max=1.6))
    def xmm(x): return left+(x+.4)/4.8*49
    for x0,x1,label in [(-.3,1.3,'Records'),(1.7,4.39,'Loci')]:
        cv.text(label,xmm((x0+x1)/2),17.1,ha='center',fontsize=8)
        cv.rule(xmm(x0),xmm(x1),21,lw=.5)
cv.text('Reference model',97.5,65,ha='center',fontsize=8.5)

# B. Descriptive pooled rates use a neutral numeric matrix, with one shared scale.
cv.letter('B',2,80,'Null calibration')
b=pd.read_csv(P/'Figure_2B_plotted_values.csv')
generators=['Block translation','Track rotation']
arms=[('Records','Independent'),('Records','Cluster-preserving'),
      ('Loci','Independent'),('Loci','Cluster-preserving')]
heat_cmap=LinearSegmentedColormap.from_list('calibration_grey',['#F4F1EC','#4F5D61'])
norm=Normalize(0,30)
axB=cv.ax(22,98,66,21)
ns.clean(axB,left=False,bottom=False)
axB.set(xlim=(-.5,3.5),ylim=(1.5,-.5),yticks=[0,1],
        yticklabels=['Block\ntranslation','Track\nrotation'],xticks=range(4))
axB.tick_params(axis='both',length=0,pad=4)
axB.set_xticklabels(['Independent','Cluster-\npreserving']*2,fontsize=8,linespacing=1.18)
cell_colours=[]
for gi,g in enumerate(generators):
    for ai,(u,n) in enumerate(arms):
        r=b[(b.generator==g)&(b.counting_unit==u)&(b.spatial_reference==n)].iloc[0]
        val=r.rejection_rate_pct
        fc=heat_cmap(norm(val))
        rgb=np.asarray(fc[:3])
        lin=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
        lum=float(lin@np.array([.2126,.7152,.0722]))
        text_col='white' if (1.05/(lum+.05)) > ((lum+.05)/.05) else C['ink']
        axB.add_patch(Rectangle((ai-.48,gi-.47),.96,.94,facecolor=fc,edgecolor='none'))
        axB.text(ai,gi,f'{val:.1f}',ha='center',va='center',fontsize=8.5,color=text_col)
        cell_colours.append(dict(value=val,facecolor=list(fc),text_color=text_col))
        geom.append(dict(panel='B',generator=g,unit=u,spatial_reference=n,estimate=val,
                         axis_min=0.,axis_max=30.))
for xc,label in [(38.5,'Records'),(71.5,'Loci')]:
    cv.text(label,xc,91,ha='center',fontsize=8.5)
cbax=cv.ax(22,132,66,2.6)
cb=mpl.colorbar.ColorbarBase(cbax,cmap=heat_cmap,norm=norm,orientation='horizontal',ticks=[0,5,15,30])
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8,width=.55,length=2,pad=2)
cb.set_label('Rejection rate (%)',fontsize=8.5,labelpad=3)

# C. Wilson confidence intervals are whiskers; all four settings remain distinct.
cv.letter('C',96,80,'Calibration by setting')
c=pd.read_csv(P/'Figure_2C_plotted_values.csv')
settings=[('HBV',100000),('HPV',100000),('HBV',1000000),('HPV',1000000)]
key=cv.ax(111,86,65,7)
key.set_axis_off()
legC=key.legend(handles=[
    mpl.lines.Line2D([],[],marker=ns.MARK['HBV'],ls='none',ms=3.7,mfc=C['HBV'],mec=C['HBV'],label='HBV'),
    mpl.lines.Line2D([],[],marker=ns.MARK['HPV'],ls='none',ms=3.7,mfc=C['HPV'],mec=C['HPV'],label='HPV'),
    mpl.lines.Line2D([],[],marker='o',ls='none',ms=3.7,mfc=C['ink'],mec=C['ink'],label='±100 kb'),
    mpl.lines.Line2D([],[],marker='o',ls='none',ms=3.7,mfc='white',mec=C['ink'],label='±1 Mb')],
    loc='center left',ncol=2,fontsize=8,columnspacing=1.5,handlelength=1,handletextpad=.45)
axesC=[]
for gi,g in enumerate(generators):
    ax=cv.ax(112,98+51*gi,64,34)
    axesC.append(ax)
    ns.clean(ax)
    for ai,(u,n) in enumerate(arms):
        for si,(v,w) in enumerate(settings):
            r=c[(c.generator==g)&(c.counting_unit==u)&(c.spatial_reference==n)&(c.virus==v)&(c.window_bp==w)].iloc[0]
            x=ai+(si-1.5)*.19
            ax.errorbar(x,r.rejection_rate_pct,yerr=[[r.rejection_rate_pct-r.lower_pct],[r.upper_pct-r.rejection_rate_pct]],
                color=C[v],elinewidth=.85,capsize=1.5,capthick=.65,marker=ns.MARK[v],ms=3.6,
                mfc=C[v] if w==100000 else 'white',mec=C[v],mew=.85,ls='none',zorder=3)
            geom.append(dict(panel='C',virus=v,window_bp=w,generator=g,unit=u,spatial_reference=n,
                         estimate=r.rejection_rate_pct,lower=r.lower_pct,upper=r.upper_pct,
                         axis_min=-2.,axis_max=48.))
    ax.axhline(5,**refline)
    ax.set(xlim=(-.6,3.6),ylim=(-2,48),yticks=[0,10,20,30,40],xticks=range(4))
    ax.set_xticklabels(['Independent','Cluster-\npreserving']*2 if gi else ['']*4,
                       fontsize=8,linespacing=1.18)
    ax.set_title(g,loc='left',fontsize=8,color=C['mid'],pad=4)
fig.text(99/180,1-139/201,'Rejection rate (%)',fontsize=8.5,ha='center',va='center',rotation=90)
for u,x in [('Records',128.76),('Loci',159.24)]:
    cv.text(u,x,195,ha='center',fontsize=8.5)

# D. Catalogue means shown as bars with a zero baseline.
cv.letter('D',2,148,'Coordinates per cluster')
d=pd.read_csv(P/'Figure_2D_plotted_values.csv')
axD=cv.ax(22,155,66,28)
ns.clean(axD,left=False)
d_positions=[0.,.9,1.8,3.7,4.6]
for i,r in enumerate(d.itertuples(index=False)):
    y=d_positions[i]
    bar=axD.barh(y,r.mean_coordinates_per_cluster,height=.42,color=ns.tint(C[r.virus],.25),edgecolor='none',zorder=3)[0]
    bar.set_gid(f'D_bar_{i}')
    axD.text(r.mean_coordinates_per_cluster+.08,y,f'{r.mean_coordinates_per_cluster:.2f}',fontsize=8,va='center')
    geom.append(dict(panel='D',estimate=r.mean_coordinates_per_cluster,axis_min=0.,axis_max=2.85))
axD.set(ylim=(5.1,-1.3),xlim=(0,2.85),yticks=d_positions,yticklabels=d.catalogue_label,
        xticks=[0,1,2],xlabel='Mean coordinates per cluster')
axD.tick_params(axis='y',length=0,pad=3)
axD.text(0,-.8,'Primary',fontsize=8,color=C['mid'],va='bottom')
axD.text(0,3.15,'External HPV',fontsize=8,color=C['mid'],va='bottom')
lower_bound=axD.axvline(1,color=C['soft'],lw=.6,ls=(0,(1,2)),zorder=4)
lower_bound.set_gid('D_lower_bound')

fig.canvas.draw()
renderer=fig.canvas.get_renderer()
legend_clearance={'C':(legC.get_window_extent(renderer).y0-axesC[0].title.get_window_extent(renderer).y1)/fig.dpi*25.4}
qa=ns.export(cv,OUT,'Yan_Fig2',extra=dict(
    panel_sources=[f'plotted_data/Figure_2{p}_plotted_values.csv' for p in 'ABCD'],
    plotted_geometry=geom,legend_clearance_mm=legend_clearance,
    estimates_by_panel={'A':15,'B':8,'C':32,'D':5},
    spatial_reference_range_style='discrete vertical capped ranges at all five categorical schemes',
    confidence_interval_style='vertical capped whiskers',
    pooled_rate_style='numeric matrix; continuous neutral scale 0 to 30 percent',
    catalogue_mean_style='bars with zero baseline; labelled primary/external groups; dotted structural lower bound at one',
    pooled_rates_are_descriptive=True,per_setting_intervals='95% binomial Wilson; 150 simulations',
    heatmap_cells=cell_colours,
    alignment_mm={'B_top':98,'C_upper_top':98,'D_bottom':183,'C_lower_bottom':183,
                  'B_left':22,'D_left':22,'B_right':88,'D_right':88},
    profile_paths=profile_paths,shared_enrichment_ylim=[.8,1.6],
    profile_axes_mm=[dict(left=19+54*i,width=49,height=33,ylim=[.8,1.6]) for i in range(3)],
    catalogue_scale=[0,2.85],
    catalogue_group_labels=['Primary','External HPV'],catalogue_row_positions=d_positions,
    catalogue_structural_lower_bound=1.,
    refined_footnote='a: joint matching and peak exclusion use a smaller eligible coordinate set',
    palette_source='../shared_style/nstyle.py',ebv_colour=C['EBV'],
    new_statistical_tests=False))
print({k:qa[k] for k in ['min_font_pt','text_overlaps','text_outside_canvas','label_warnings','legend_clearance_mm']})

