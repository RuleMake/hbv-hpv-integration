"""Window sensitivity and complete single-cluster influence summaries."""
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import nstyle as ns
from nstyle import C
ROOT=Path(__file__).resolve().parent
a,b,c=[pd.read_csv(ROOT/f'plotted_data/S4{p}_display.csv',float_precision='round_trip') for p in 'ABC']
VIR=['HPV','HBV','EBV'];CLS=['DEL','DUP','INV','TRA']
cv=ns.Canvas(180,283);geometry=[]
cv.text('Supplemental Figure S4',178,1,ha='right',fontsize=9,fontweight='bold')
cv.letter('A',2,10,'Enrichment at two window sizes')
for i,v in enumerate(VIR):
    x=14+60*i;ax=cv.ax(x,29,42,41);ns.clean(ax);ns.href(ax,1)
    lim=(.58,1.85) if v=='EBV' else (.92,1.52)
    for j,cl in enumerate(CLS):
        for k,window in enumerate(['100 kb','25 kb']):
            r=a[a.virus.eq(v)&a.svclass.eq(cl)&a.window.eq(window)].iloc[0]
            xx=j+[-.14,.14][k];gid=f'A_{v}_{cl}_{k}'
            line,=ax.plot([xx,xx],[r.ci_lo,r.ci_hi],color=C[v],lw=.95);line.set_gid(gid+'_range')
            mark,=ax.plot(xx,r.fold,ls='none',marker=ns.MARK[v],ms=4,mfc=C[v] if k==0 else 'white',mec=C[v],mew=.85);mark.set_gid(gid+'_point')
            geometry.append(dict(panel='A',gid=gid,virus=v,svclass=cl,window=window,x=xx,value=r.fold,lo=r.ci_lo,hi=r.ci_hi,axis_mm=[x,29,42,41],xlim=[-.5,3.5],ylim=list(lim)))
    ax.set(xlim=(-.5,3.5),ylim=lim,xticks=range(4),xticklabels=CLS,yticks=[.6,1,1.4,1.8] if v=='EBV' else [1,1.2,1.4])
    cv.text(v,x,21,fontsize=9,fontweight='bold');cv.text(f'{int(a[a.virus.eq(v)].n_coords.iloc[0]):,} coordinates',x+42,21,ha='right',color=C['mid'])
cv.fig.text(2/180,1-49.5/283,'Observed / matched null',rotation=90,ha='center',va='center',fontsize=8.5)
leg=cv.ax(10,79,165,7);leg.set_axis_off()
leg.legend(handles=[Line2D([],[],ls='none',marker='o',color=C['mid'],label='±100 kb'),
    Line2D([],[],ls='none',marker='o',mfc='white',mec=C['mid'],label='±25 kb'),
    Line2D([],[],color=C['mid'],label='95% spatial reference range')],loc='center',ncol=3,handlelength=1.3,columnspacing=1.5)

cv.letter('B',2,95,'Leave-one-cluster-out influence at ±25 kb')
order=[('HPV','DEL'),('HPV','DUP'),('HBV','TRA'),('HBV','INV')]
ax=cv.ax(26,109,150,30);ns.clean(ax);ns.href(ax,0)
for i,(v,cl) in enumerate(order):
    r=b[b.virus.eq(v)&b.svclass.eq(cl)].iloc[0]
    patch=Rectangle((i-.12,r.change_min),.24,r.change_max-r.change_min,fc=ns.tint(C[v],.63),ec=C[v],lw=.75);patch.set_gid(f'B_range_{i}');ax.add_patch(patch)
    for suffix,y in [('lo',r.change_min),('hi',r.change_max)]:
        line,=ax.plot([i-.16,i+.16],[y,y],color=C[v],lw=.85);line.set_gid(f'B_{suffix}_{i}')
    cv.text(f'{r.fold:.4f}',26+(i+.5)/4*150,153,ha='center',va='center');cv.text(f'{int(r.n_clusters):,}',26+(i+.5)/4*150,160,ha='center',va='center')
    geometry.append(dict(panel='B',gid=f'B_range_{i}',virus=v,svclass=cl,x=i,lo=r.change_min,hi=r.change_max,full_fold=r.fold,n_clusters=int(r.n_clusters),axis_mm=[26,109,150,30],xlim=[-.5,3.5],ylim=[-.025,.012],interval='min-max over all leave-one-cluster removals'))
ax.set(xlim=(-.5,3.5),ylim=(-.025,.012),xticks=range(4),xticklabels=[v+' '+cl for v,cl in order],yticks=[-.02,-.01,0,.01])
ax.set_ylabel('Change in fold',labelpad=3)
cv.text('Full-set fold',3,153,va='center',color=C['mid']);cv.text('Clusters',3,160,va='center',color=C['mid'])

cv.letter('C',2,176,'Largest cluster contributions at ±25 kb')
ax=cv.ax(65,190,107,68);ns.clean(ax,left=False)
ys=[];labels=[]
for gi,(v,cl) in enumerate([('HPV','DEL'),('HPV','DUP'),('HBV','TRA')]):
    g=c[c.virus.eq(v)&c.svclass.eq(cl)].reset_index(drop=True);assert len(g)==6 and g.excess_share_pct.notna().all()
    slot=gi*7
    if gi:ax.axhline(slot-.9,color=C['rule'],lw=.55)
    cv.text(v+'\n'+cl,3,190+(slot+3)/20*68,va='center',fontsize=8.5,fontweight='bold',linespacing=1.25)
    for i,r in g.iterrows():
        y=slot+i;ys.append(y);labels.append(r.note.split(';')[0])
        p=ax.barh(y,r.excess_share_pct,height=.52,left=0,color=ns.tint(C[v],.18),edgecolor=C[v],linewidth=.55)[0];p.set_gid(f'C_bar_{gi}_{i}')
        txt=ax.text(r.excess_share_pct+.10,y,f'{r.excess_share_pct:.2f}',va='center',ha='left');txt.set_gid(f'C_value_{gi}_{i}')
        geometry.append(dict(panel='C',gid=f'C_bar_{gi}_{i}',virus=v,svclass=cl,coordinate=r.note.split(';')[0],y=y,value=r.excess_share_pct,left=0,axis_mm=[65,190,107,68],xlim=[0,6.65],ylim=[19.5,-.5]))
ax.set(ylim=(19.5,-.5),yticks=ys,yticklabels=labels,xlim=(0,6.65),xticks=[0,2,4,6]);ax.tick_params(axis='y',pad=3)
ax.set_xlabel('Share of class excess (%)',labelpad=4)
cv.text('HBV INV',3,276,va='center',fontsize=8.5,fontweight='bold')
cv.text('Contribution shares not defined',65,276,va='center',color=C['mid'])
qa=ns.export(cv,ROOT/'rendered','Supplemental_Fig_S4',extra=dict(status='saved_result_redraw',plotted_geometry=geometry,
    spatial_reference_ranges=24,leave_one_out_ranges=4,defined_cluster_shares=18,undefined_cluster_shares_retained=6,
    transformations=['A original folds and full reference intervals','B saved leave-one-out min/max minus full-set fold','C zero-based bars from stored contribution percentages'],
    undefined_INV_shares='not displayed as zero; retained in source and plotted tables',new_statistical_tests=False))
print({k:qa[k] for k in ['text_overlaps','text_outside_canvas','label_warnings']})
