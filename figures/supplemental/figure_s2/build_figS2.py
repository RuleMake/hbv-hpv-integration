"""Supplemental Figure S2: calibration, pool availability and fold spread."""
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import nstyle as ns
from nstyle import C

ROOT=Path(__file__).resolve().parent
data=ROOT/'plotted_data'
cal=pd.read_csv(data/'S2ABC_calibration.csv')
reject=pd.read_csv(data/'S2D_rejection_intervals.csv')
pool=pd.read_csv(data/'S2E_pool_diagnostics.csv')
fold=pd.read_csv(data/'S2F_fold_summaries.csv')
cv=ns.Canvas(180,263)
cv.text('Supplemental Figure S2',178,1,ha='right',fontsize=9,fontweight='bold')
geometry=[]
rules=['published','dedup','corrected']

# A-C share scales; window facets keep each axes to two virus curves.
for i,(rule,counting,null) in enumerate(zip(rules,['Records','Loci','Loci'],['Independent null','Independent null','Cluster-preserving null'])):
    x=2+60*i
    cv.letter('ABC'[i],x,10,counting)
    cv.text(null,x+4.6,16,color=C['mid'])
    for wi,window in enumerate([100000,1000000]):
        top=24+27*wi
        ax=cv.ax(14+60*i,top,43,22);ns.clean(ax)
        ax.plot([0,.55],[0,55],color=C['soft'],lw=.75,ls=':',zorder=1)
        ax.text(.045,.93,'±100 kb' if wi==0 else '±1 Mb',transform=ax.transAxes,
                va='top',fontsize=8,color=C['mid'])
        for virus in ['HBV','HPV']:
            r=cal[cal.rule.eq(rule)&cal.virus.eq(virus)&cal.window_bp.eq(window)].sort_values('alpha')
            ln,=ax.plot(r.alpha,100*r.rejection_fraction,color=C[virus],ls='-' if window==100000 else '--',
                marker=ns.MARK[virus],ms=3.4,mfc=C[virus] if window==100000 else 'white',mew=.85,lw=1.05,zorder=3)
            ln.set_gid(f'{"ABC"[i]}_{virus}_{window}')
            geometry.append(dict(panel='ABC'[i],rule=rule,virus=virus,window_bp=window,
                x=r.alpha.tolist(),y=(100*r.rejection_fraction).tolist(),axis_mm=[14+60*i,top,43,22],xlim=[-.015,.56],ylim=[-2,65]))
        ax.set(xlim=(-.015,.56),ylim=(-2,65),xticks=[0,.1,.3,.5],yticks=[0,30,60])
        if wi==0:ax.tick_params(axis='x',labelbottom=False)
        else:ax.set_xlabel('Nominal α',labelpad=2.5)
cv.fig.text(2/180,1-48.5/263,'Rejection rate (%)',rotation=90,ha='center',va='center',fontsize=8.5)
legend=cv.ax(4,87,173,7);legend.set_axis_off()
legend.legend(handles=[Line2D([],[],color=C[v],marker=ns.MARK[v],ls='none',ms=4,label=v) for v in ['HBV','HPV']]+[
    Line2D([],[],color=C['mid'],marker='o',ms=3.2,mfc=C['mid'],label='±100 kb'),
    Line2D([],[],color=C['mid'],ls='--',marker='o',ms=3.2,mfc='white',label='±1 Mb'),
    Line2D([],[],color=C['soft'],ls=':',label='Expected under null')],
    loc='center',ncol=5,handlelength=1.7,columnspacing=1.35)

# D: common-scale comparisons with separate magnified views of the locus schemes.
cv.letter('D',2,100,'Rejection rates and confidence intervals')
settings=[('HBV',100000),('HPV',100000),('HBV',1000000),('HPV',1000000)]
for ai,alpha in enumerate([.05,.01]):
    main_x=13+90*ai
    zoom_x=43+90*ai
    cv.text(f'α = {alpha:.2f}',main_x,111,fontsize=8.5)
    cv.text('Locus detail',zoom_x,111,fontsize=8,color=C['mid'])
    ax=cv.ax(main_x,119,74,43);ns.clean(ax)
    zoom=cv.ax(zoom_x,119,42,25);ns.clean(zoom)
    zoom.set_zorder(4)
    ax.add_patch(Rectangle((.54,-1.3),1.9,11.8,fc='none',ec=C['rule'],lw=.6,zorder=0))
    for j,rule in enumerate(rules):
        for si,(virus,window) in enumerate(settings):
            r=reject[reject.rule.eq(rule)&reject.virus.eq(virus)&reject.window_bp.eq(window)&np.isclose(reject.alpha,alpha)].iloc[0]
            x=j+(si-1.5)*.18
            for target,panel,rect,ylim in [(ax,'D',[main_x,119,74,43],[-2,46])]+(
                    [(zoom,'D_zoom',[zoom_x,119,42,25],[-.5,10.5])] if rule!='published' else []):
                ln,=target.plot([x,x],[r.wilson_lo*100,r.wilson_hi*100],color=C[virus],lw=.9,zorder=3)
                ln.set_gid(f'{panel}_CI_{ai}_{rule}_{virus}_{window}')
                mark,=target.plot(x,r.rate*100,marker=ns.MARK[virus],ms=3.6,
                    mfc=C[virus] if window==100000 else 'white',mec=C[virus],mew=.85,ls='none',zorder=4)
                mark.set_gid(f'{panel}_rate_{ai}_{rule}_{virus}_{window}')
                geometry.append(dict(panel=panel,alpha=alpha,rule=rule,virus=virus,window_bp=window,
                    value=r.rate*100,lo=r.wilson_lo*100,hi=r.wilson_hi*100,x=x,
                    axis_mm=rect,xlim=[-.5,2.5] if panel=='D' else [.5,2.5],ylim=ylim,
                    duplicate_display=panel=='D_zoom',interval='95% Wilson CI'))
    ax.set(xlim=(-.5,2.5),ylim=(-2,46),xticks=[0,1,2],yticks=[0,10,20,30,40])
    ax.set_xticklabels(['Records\nIndependent','Loci\nIndependent','Loci\nCluster'],linespacing=1.15)
    if ai==0:ax.set_ylabel('Rejection rate (%)',labelpad=3)
    zoom.set(xlim=(.5,2.5),ylim=(-.5,10.5),xticks=[1,2],yticks=[0,5,10])
    zoom.set_xticklabels(['Independent','Cluster'])

# E uses the actual window half-width on a labelled log axis.
cv.letter('E',2,180,'Matched-null pool')
ax=cv.ax(15,195,40,53);ns.clean(ax)
for virus in ['HBV','HPV','EBV']:
    r=pool[pool.virus.eq(virus)].sort_values('window_bp')
    ln,=ax.plot(r.window_bp,r.clean_pool_pct,color=C[virus],marker=ns.MARK[virus],ms=4,lw=1.05,mew=.8)
    ln.set_gid(f'E_pool_{virus}')
    geometry.append(dict(panel='E',virus=virus,x=r.window_bp.tolist(),y=r.clean_pool_pct.tolist(),
        axis_mm=[15,195,40,53],xscale='log10',xlim=[80000,6500000],ylim=[-3,105]))
ax.set_xscale('log');ax.minorticks_off()
ax.set(xlim=(80000,6500000),ylim=(-3,105),xticks=[100000,1000000,5000000],yticks=[0,25,50,75,100])
ax.set_xticklabels(['100 kb','1 Mb','5 Mb'])
ax.set_xlabel('Window half-width',labelpad=3)
ax.set_ylabel('Pool outside observed windows (%)',labelpad=3)
ebv_last=pool[pool.virus.eq('EBV')].sort_values('window_bp').iloc[-1]
ax.annotate('EBV',xy=(ebv_last.window_bp,ebv_last.clean_pool_pct),xytext=(0,-15),textcoords='offset points',ha='right',color=C['EBV'],fontsize=8)

# F contains saved quantiles and full ranges, not estimated densities or CIs.
cv.letter('F',63,180,'Null folds and observed enrichment')
legend=cv.ax(63,186,115,8);legend.set_axis_off()
legend.legend(handles=[Patch(fc=C['wash2'],ec=C['mid'],label='5th–95th percentile'),
    Line2D([],[],color=C['mid'],label='Full range'),
    Line2D([],[],ls='none',marker='|',color=C['mid'],label='Median')],
    loc='center',ncol=3,handlelength=1.25,columnspacing=1.2)
cv.text('Null',130.5,198,ha='center',fontsize=8,color=C['mid'])
cv.text('Observed',178,198,ha='right',fontsize=8,color=C['mid'])
ax=cv.ax(104,204,53,44);ns.clean(ax,left=False);ns.vref(ax,1,lw=.7)
order=[(rule,virus,window) for rule in ['published','corrected'] for virus in ['HBV','HPV'] for window in [100000,1000000]]
for i,(rule,virus,window) in enumerate(order):
    r=fold[fold.rule.eq(rule)&fold.virus.eq(virus)&fold.window_bp.eq(window)].iloc[0]
    prefix=f'F_{rule}_{virus}_{window}'
    ln,=ax.plot([r.fold_min,r.fold_max],[i,i],color=C[virus],lw=.85);ln.set_gid(prefix+'_range')
    ax.plot([r.fold_min,r.fold_max],[i,i],ls='none',marker='|',ms=4,color=C[virus])
    patch=Rectangle((r.fold_p05,i-.23),r.fold_p95-r.fold_p05,.46,fc=ns.tint(C[virus],.67),ec=C[virus],lw=.65,zorder=3)
    patch.set_gid(prefix+'_percentiles');ax.add_patch(patch)
    ln,=ax.plot([r.fold_median]*2,[i-.27,i+.27],color=C[virus],lw=1.25,zorder=4);ln.set_gid(prefix+'_median')
    cv.text('100 kb' if window==100000 else '1 Mb',100,204+(i+.5)/8*44,ha='right',va='center',fontsize=8)
    geometry.append(dict(panel='F',rule=rule,virus=virus,window_bp=window,row=i,
        minimum=r.fold_min,p05=r.fold_p05,median=r.fold_median,p95=r.fold_p95,maximum=r.fold_max,
        axis_mm=[104,204,53,44],xlim=[.90,1.15],ylim=[7.5,-.5]))
for _,r in fold[fold.rule.str.startswith('observed')].iterrows():
    i=order.index(('published',r.virus,int(r.window_bp)))
    text=cv.text(f'{r.fold_median:.2f}',178,204+(i+.5)/8*44,ha='right',va='center',fontsize=8.5,fontweight='bold',color=C['ink'])
    text.set_gid(f'F_observed_{r.virus}')
    geometry.append(dict(panel='F_observed',rule=r.rule,virus=r.virus,window_bp=int(r.window_bp),value=r.fold_median,row=i))
ax.axhline(3.5,color=C['rule'],lw=.55)
ax.set(xlim=(.90,1.15),ylim=(7.5,-.5),xticks=[.9,1,1.1,1.15],yticks=[])
ax.set_xlabel('Observed / expected overlap',labelpad=3)
for gi,label in enumerate(['Records\nIndep.','Loci\nCluster']):
    cv.text(label,63,204+(gi*4+2)/8*44,va='center',linespacing=1.25,color=C['mid'])
    for vi,virus in enumerate(['HBV','HPV']):
        cv.text(virus,79,204+(gi*4+vi*2+1)/8*44,va='center',fontsize=8)

qa=ns.export(cv,ROOT/'rendered','Supplemental_Fig_S2',extra=dict(plotted_geometry=geometry,
    calibration_points=48,wilson_intervals=24,wilson_interval_displays=40,zoom_duplicate_intervals=16,pool_points=9,fold_ranges=8,observed_folds=2,
    interval_semantics={'D':'95% Wilson binomial confidence intervals','F':'5th–95th percentile and full range of null replicates, not confidence intervals'},
    transformations=['A-C fractions converted to percent','D fraction to percent','E log10 window half-width axis; no smoothing','F saved summary values'],
    new_statistical_tests=False))
print({k:qa[k] for k in ['min_font_pt','text_overlaps','text_outside_canvas','label_warnings']})
