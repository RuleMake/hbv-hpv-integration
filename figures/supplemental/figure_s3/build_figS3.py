"""Supplemental Figure S3: complete class atlas, power scale and cluster contribution."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import nstyle as ns
from nstyle import C
ROOT=Path(__file__).resolve().parent
a,b,c,d=[pd.read_csv(ROOT/f'plotted_data/S3{p}_display.csv',float_precision='round_trip') for p in 'ABCD']
VIR=['HPV','HBV','EBV'];CLS=['DEL','DUP','INV','TRA'];EST=['bc_uniqSV_blockNull','f_clusterSignTest']
cv=ns.Canvas(180,232);geometry=[]
cv.text('Supplemental Figure S3',178,1,ha='right',fontsize=9,fontweight='bold')
cv.letter('A',2,10,'Class enrichment at ±100 kb')
for i,v in enumerate(VIR):
    x=14+60*i;ax=cv.ax(x,29,42,46);ns.clean(ax)
    ns.href(ax,1)
    lim=(.70,1.65) if v=='EBV' else (.88,1.30)
    for j,cl in enumerate(CLS):
        for k,est in enumerate(EST):
            r=a[a.virus.eq(v)&a.svclass.eq(cl)&a.estimator.eq(est)].iloc[0]
            xx=j+[-.14,.14][k];gid=f'A_{v}_{cl}_{k}'
            ln,=ax.plot([xx,xx],[r.ci_lo,r.ci_hi],color=C[v],lw=.95);ln.set_gid(gid+'_range')
            p,=ax.plot(xx,r.fold,ls='none',marker=ns.MARK[v],ms=4,mfc=C[v] if k==0 else 'white',mec=C[v],mew=.85);p.set_gid(gid+'_point')
            if r.passes_two_sided_bonferroni_24:
                t=ax.text(xx,r.ci_hi+.012,'*',ha='center',va='bottom',fontweight='bold');t.set_gid(gid+'_star')
            geometry.append(dict(panel='A',gid=gid,virus=v,svclass=cl,estimator=est,x=xx,value=r.fold,lo=r.ci_lo,hi=r.ci_hi,axis_mm=[x,29,42,46],xlim=[-.5,3.5],ylim=list(lim)))
    ax.set(xlim=(-.5,3.5),ylim=lim,xticks=range(4),xticklabels=CLS,
           yticks=[.8,1,1.2,1.4,1.6] if v=='EBV' else [.9,1,1.1,1.2,1.3])
    cv.text(v,x,21,fontsize=9,fontweight='bold')
    cv.text(f'{int(a[a.virus.eq(v)].n_coords.iloc[0]):,} coordinates',x+42,21,ha='right',color=C['mid'])
cv.fig.text(2/180,1-52/232,'Observed / matched null',rotation=90,ha='center',va='center',fontsize=8.5)
leg=cv.ax(10,85,165,7);leg.set_axis_off()
leg.legend(handles=[Line2D([],[],ls='none',marker='o',color=C['mid'],label='Record count'),
    Line2D([],[],ls='none',marker='o',mfc='white',mec=C['mid'],label='Cluster sign'),
    Line2D([],[],color=C['mid'],label='95% spatial reference range')],loc='center',ncol=3,handlelength=1.3,columnspacing=1.2)

cv.letter('B',2,101,'HPV–HBV enrichment contrast')
ax=cv.ax(27,115,45,54);ns.clean(ax,left=False);ns.vref(ax,0)
keys=CLS+['(DEL+DUP)-(TRA+INV)']
for i,key in enumerate(keys):
    r=b[b.svclass.eq(key)].iloc[0]
    patch=Rectangle((-r.mde,i-.24),2*r.mde,.48,fc=C['wash2'],ec='none',zorder=1);patch.set_gid(f'B_MDE_{i}');ax.add_patch(patch)
    p,=ax.plot(r.log_ratio,i,ls='none',marker='o',ms=4.2,color=C['ink'],zorder=3);p.set_gid(f'B_point_{i}')
    cv.text('<0.001' if r.emp_p<.001 else f'{r.emp_p:.3f}',88,115+(i+.55)/5.1*54,ha='right',va='center',color=C['mid'])
    geometry.append(dict(panel='B',gid=f'B_point_{i}',svclass=key,value=r.log_ratio,mde=r.mde,row=i,axis_mm=[27,115,45,54],xlim=[-.17,.255],ylim=[4.55,-.55]))
ax.set(ylim=(4.55,-.55),yticks=range(5),yticklabels=CLS+['Class\ncontrast'],xlim=(-.17,.255),xticks=[-.1,0,.1,.2])
ax.set_xlabel('Log(HPV / HBV)',labelpad=4)
cv.text('P',88,109,ha='right',fontstyle='italic',color=C['mid'])
cv.text('Approx. MDE',27,185,color=C['mid'])
cv.fig.patches.append(Rectangle((16/180,1-187.4/232),7/180,2.4/232,transform=cv.fig.transFigure,fc=C['wash2'],ec='none'))

cv.letter('C',95,101,'Effect size and sensitivity')
ax=cv.ax(117,115,59,54);ns.clean(ax);ax.set_xscale('log');ax.set_yscale('log')
limx=(.039,.385);limy=(.0048,.228);xx=np.geomspace(*limx,200)
ax.fill_between(xx,limy[0],np.clip(xx,*limy),color=C['wash2'],zorder=0);ax.plot(xx,xx,color=C['soft'],lw=.8)
for i,r in c.iterrows():
    p,=ax.plot(r.mde_half_width,r.effect_half_width,ls='none',marker=ns.MARK[r.virus],ms=3.2,mfc=C[r.virus] if r.estimator==EST[0] else 'white',mec=C[r.virus],mew=.8);p.set_gid(f'C_point_{i}')
    geometry.append(dict(panel='C',gid=f'C_point_{i}',virus=r.virus,svclass=r.svclass,estimator=r.estimator,x=r.mde_half_width,y=r.effect_half_width,axis_mm=[117,115,59,54],xlim=list(limx),ylim=list(limy),scale='log10 x and y'))
ax.text(.18,.007,'Below MDE',ha='center',va='center',fontsize=8,color=C['mid'])
ax.set_xlim(*limx);ax.set_ylim(*limy)
for axis,ticks in [(ax.xaxis,[.05,.1,.2]),(ax.yaxis,[.005,.01,.02,.05,.1,.2])]:
    axis.set_major_locator(mpl.ticker.FixedLocator(ticks));axis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p:f'{x:g}'));axis.set_minor_locator(mpl.ticker.NullLocator())
ax.set_xlabel('Approximate MDE',labelpad=4);ax.set_ylabel('Absolute deviation from 1',labelpad=4)
cv.text('Line: deviation = MDE',117,185,color=C['mid'])

cv.letter('D',2,197,'Share contributed by the ten leading clusters (%)')
ax=cv.ax(29,210,147,17);ax.set_axis_off();ax.set_xlim(-.5,3.5);ax.set_ylim(2.5,-.5)
for j,cl in enumerate(CLS):cv.text(cl,29+(j+.5)/4*147,205,ha='center',fontsize=8,fontweight='bold')
for i,v in enumerate(VIR):
    cv.text(v,3,210+(i+.5)/3*17,va='center',fontsize=8)
    for j,cl in enumerate(CLS):
        r=d[d.virus.eq(v)&d.svclass.eq(cl)].iloc[0];defined=np.isfinite(r.top10_share_pct)
        color=ns.tint(C['mid'],.95-.45*r.top10_share_pct/100) if defined else C['wash']
        p=Rectangle((j-.48,i-.44),.96,.88,fc=color,ec='white',lw=.8);p.set_gid(f'D_cell_{v}_{cl}');ax.add_patch(p)
        txt=ax.text(j,i,f'{r.top10_share_pct:.1f}' if defined else '—',ha='center',va='center',color=C['ink'] if defined else C['mid']);txt.set_gid(f'D_value_{v}_{cl}')
        geometry.append(dict(panel='D',gid=f'D_cell_{v}_{cl}',virus=v,svclass=cl,value=float(r.top10_share_pct) if defined else None,defined=bool(defined)))
qa=ns.export(cv,ROOT/'rendered','Supplemental_Fig_S3',extra=dict(status='saved_result_redraw',plotted_geometry=geometry,
    spatial_ranges=24,bonferroni_24_stars=int(a.passes_two_sided_bonferroni_24.sum()),contrast_points=5,contrast_MDE_windows=5,effect_MDE_points=24,defined_top10_shares=5,undefined_top10_shares=7,
    transformations=['A explicit two-sided p and Bonferroni-24 fields','B saved log ratios and zero-centred MDE rectangles','C log10 axes only','D exact share matrix; missing remains undefined'],new_statistical_tests=False))
print({k:qa[k] for k in ['text_overlaps','text_outside_canvas','label_warnings']})
