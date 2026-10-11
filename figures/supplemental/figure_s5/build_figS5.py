"""S5: four class facets, hotspot membership and two definition comparisons."""
from pathlib import Path
import pandas as pd
import matplotlib as mpl
from matplotlib.patches import Patch
import nstyle as ns
from nstyle import C

ROOT=Path(__file__).resolve().parent
d=pd.read_csv(ROOT/'plotted_data/S5_plotted_values.csv')
pa=d[d.panel=='a']; pc=d[d.panel=='c']
cv=ns.Canvas(180,218)
cv.text('Supplemental Figure S5',178,1,ha='right',fontsize=9,fontweight='bold')
cv.letter('A',2,10,'Class enrichment by hotspot status')
groups=[('HPV','DEL'),('HPV','DUP'),('HBV','INV'),('HBV','TRA')]
intervals=[]
for j,(virus,cls) in enumerate(groups):
    x=19+j*41
    ax=cv.ax(x,25,34,51)
    ns.clean(ax,left=j==0); ns.href(ax,1,ls='--')
    ax.set_yscale('log');ax.set_ylim(.90,6.1);ax.set_xlim(-.5,2.5)
    ax.set_yticks([1,1.5,2,3,5]);ax.set_yticklabels(['1','1.5','2','3','5'] if j==0 else [])
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_xticks(range(3),['All','Hotspots','Other'])
    cv.text(f'{virus}  {cls}',x+17,18.5,ha='center',fontweight='bold',color=ns.shade(C[virus],.15))
    for k,key in enumerate(['ALL','HOTSPOT','NON-HOTSPOT']):
        r=pa[(pa.virus==virus)&(pa.svclass==cls)&(pa.subset==key)].iloc[0]
        ax.fill_between([k-.10,k+.10],r.ci_lo,r.ci_hi,facecolor=ns.tint(C[virus],.84),edgecolor=ns.tint(C[virus],.25),lw=.75,zorder=2)
        ax.plot(k,r.fold,marker=ns.MARK[virus],ms=3.6,color=C[virus],ls='none',zorder=4)
        intervals.append(dict(panel='A',virus=virus,svclass=cls,subset=key,lower=float(r.ci_lo),upper=float(r.ci_hi),limits=[.90,6.1],complete=bool(.90<=r.ci_lo<=r.ci_hi<=6.1)))
    if j==0:ax.set_ylabel('Observed / matched null',labelpad=3)

cv.letter('B',2,91,'Hotspot membership')
ax=cv.ax(57,108,70,31)
ns.clean(ax,left=False);ax.set_xlim(0,8.4);ax.set_ylim(3.6,-.6)
ax.set_yticks(range(4),['Recurrence','Gene list','Recurrence','Gene list'])
ax.set_xticks([0,2,4,6,8]);ax.set_xlabel('Hotspot share (%)',labelpad=2)
cv.text('Coordinates',152,101,ha='right',color=C['mid'])
cv.text('Clusters',177,101,ha='right',color=C['mid'])
for gi,v in enumerate(['HPV','HBV']):
    cv.text(v,4,108+(gi*2+1.1)/4.2*31,va='center',fontweight='bold',color=ns.shade(C[v],.15))
    for k,key in enumerate(['recurrence','gene list']):
        i=gi*2+k;r=pc[(pc.virus==v)&(pc.subset==f'HOTSPOT ({key})')].iloc[0]
        ax.barh(i-.17,r.coordinate_share_pct,height=.28,color=C[v],edgecolor=C[v],lw=.7)
        ax.barh(i+.17,r.cluster_share_pct,height=.28,color='white',edgecolor=C[v],lw=.8)
        tick,=ax.plot(r.cluster_share_pct,i+.17,marker='|',ms=6,color=C[v],mew=1.0,ls='none',zorder=4)
        tick.set_gid(f'B_cluster_endpoint_{i}')
        y=108+(i+.6)/4.2*31
        cv.text(f'{int(r.n_coords):,} / {int(r.total_coordinates):,}',152,y,ha='right',va='center')
        cv.text(f'{int(r.n_clusters):,} / {int(r.total_clusters):,}',177,y,ha='right',va='center')
cv.fig.legend(handles=[Patch(fc=C['mid'],ec=C['mid'],label='Coordinates'),Patch(fc='white',ec=C['mid'],label='Clusters')],loc='upper left',bbox_to_anchor=(.28,1-96/218),ncol=2,handlelength=1,columnspacing=1.4)

cv.letter('C',2,155,'Large-duplication enrichment under two hotspot definitions')
keys=['ALL (joint reference)','HOTSPOT (recurrence)','NON-HOTSPOT (recurrence)','HOTSPOT (gene list)','NON-HOTSPOT (gene list)']
labels=['All coordinates','Recurrence: hotspots','Recurrence: other coordinates','Gene list: hotspots','Gene list: other coordinates']
for gi,v in enumerate(['HPV','HBV']):
    x=58+gi*62
    ax=cv.ax(x,171,55,33)
    ns.clean(ax,left=False); ns.vref(ax,1,ls='--');ax.set_xscale('log')
    ax.set_xlim(.45,5);ax.set_ylim(4.55,-.55);ax.set_xticks([.5,1,2,4]);ax.set_xticklabels(['0.5','1','2','4'])
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.set_yticks(range(5),labels if gi==0 else ['']*5)
    ax.set_xlabel('Observed / matched null',labelpad=3)
    cv.text(v,x+27.5,165,ha='center',fontweight='bold',color=ns.shade(C[v],.15))
    for i,key in enumerate(keys):
        r=pc[(pc.virus==v)&(pc.subset==key)].iloc[0]
        ns.null_band(ax,r.ci_lo,r.ci_hi,i,half=.11,colour=ns.tint(C[v],.84))
        ax.plot(r.fold,i,marker=ns.MARK[v],ms=3.3,color=C[v],ls='none',zorder=4)
        intervals.append(dict(panel='C',virus=v,subset=key,lower=float(r.ci_lo),upper=float(r.ci_hi),limits=[.45,5],complete=bool(.45<=r.ci_lo<=r.ci_hi<=5)))
qa=ns.export(cv,ROOT/'rendered','Supplemental_Fig_S5',extra={'source':'plotted_data/S5_plotted_values.csv','interval_type':'95% spatial-null reference ranges, not donor confidence intervals','intervals':intervals,'all_intervals_complete':all(x['complete'] for x in intervals),'new_statistical_tests':False})
print({k:qa[k] for k in ('text_overlaps','text_outside_canvas','all_intervals_complete')})
