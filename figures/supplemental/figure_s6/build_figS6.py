"""S6: common-scale sensitivity facets, case ratios and exclusion extent."""
from pathlib import Path
import pandas as pd
import nstyle as ns
from nstyle import C

ROOT=Path(__file__).resolve().parent
d=pd.read_csv(ROOT/'plotted_data/S6_plotted_values.csv')
pa=d[d.panel=='a'];pb=d[d.panel=='b']
counts=pd.read_csv(ROOT/'plotted_data/S6_exclusion_accounting.csv')
cv=ns.Canvas(180,214)
cv.text('Supplemental Figure S6',178,1,ha='right',fontsize=9,fontweight='bold')
cv.letter('A',2,10,'HBV class enrichment after sensitivity exclusions')
classes=['DEL','DUP','INV','TRA']
settings=[('baseline','None'),('complex-excl','Cases'),('amplicon-masked','Mask'),('amplicon-masked+complex-excl','Both')]
intervals=[]
for j,cl in enumerate(classes):
    x=19+j*41
    ax=cv.ax(x,25,34,55)
    ns.clean(ax,left=j==0);ns.href(ax,1,ls='--')
    ax.set_xlim(-.55,3.55);ax.set_ylim(.5,4.45)
    ax.set_yticks([1,2,3,4]);ax.set_yticklabels(['1','2','3','4'] if j==0 else [])
    ax.set_xticks(range(4),[t[1] for t in settings]);ax.tick_params(axis='x',length=0)
    if j==0:ax.set_ylabel('Observed / matched null',labelpad=3)
    cv.text(cl,x+17,18.5,ha='center',fontweight='bold',color=C['HBV'])
    for k,(key,short) in enumerate(settings):
        r=pa[(pa.arm==key)&(pa.svclass==cl)].iloc[0]
        ax.fill_between([k-.15,k+.15],r.ci_lo,r.ci_hi,facecolor=ns.tint(C['HBV'],.84),edgecolor=ns.tint(C['HBV'],.20),lw=.75,zorder=2)
        ax.plot(k,r.fold,'o',ms=4.3,color=C['HBV'],zorder=3)
        intervals.append(dict(panel='A',svclass=cl,setting=key,lower=float(r.ci_lo),upper=float(r.ci_hi),limits=[.5,4.45],complete=bool(.5<=r.ci_lo<=r.ci_hi<=4.45)))
    cv.text('MDE (fold)',x+12,91,ha='right',color=C['mid'])
    cv.text(ns.IT['P'],x+34,91,ha='right',color=C['mid'])
    cv.rule(x-3,x+36,95.8)
    for k,(key,short) in enumerate(settings):
        r=pa[(pa.arm==key)&(pa.svclass==cl)].iloc[0]
        y=98+5.4*k
        if j==0:cv.text(short,2,y)
        cv.text(f'{r.mde_fold:.3f}',x+12,y,ha='right')
        cv.text('<0.0001' if r.emp_p<.0001 else f'{r.emp_p:.4f}',x+34,y,ha='right')

cv.letter('B',2,133,'Nearby-SV fractions')
ax=cv.ax(19,148,54,48)
ns.clean(ax,left=False);ns.vref(ax,1,ls='--')
ax.set_xlim(.65,2.60);ax.set_ylim(3.55,-.55)
ax.set_xticks([1,1.5,2,2.5]);ax.set_yticks(range(4),['DEL','DUP','INV','All'])
ax.set_xlabel('Ratio of per-case means',labelpad=3)
cv.text(ns.IT['P']+' (one-sided)',85,142,ha='right',color=C['mid'])
for i,cl in enumerate(['DEL','DUP','INV','ALL']):
    r=pb[pb.svclass==cl].iloc[0]
    ax.errorbar(r.fold,i,xerr=[[r.fold-r.ci_lo],[r.ci_hi-r.fold]],fmt='o',ms=4.7,color=C['HBV'],elinewidth=1.2,capsize=2.5,capthick=.8)
    cv.text(f'{r.emp_p:.4f}',85,148+(i+.55)/4.1*48,va='center',ha='right')
    intervals.append(dict(panel='B',svclass=cl,lower=float(r.ci_lo),upper=float(r.ci_hi),limits=[.65,2.6],complete=bool(.65<=r.ci_lo<=r.ci_hi<=2.6)))

cv.letter('C',94,133,'Exclusion extent')
ax=cv.ax(127,148,22,48)
ns.clean(ax,left=False)
labels=['Annotated calls','Affected cases','DEL','DUP','INV','TRA','All classes']
ax.set_yticks(range(7),labels);ax.set_ylim(6.55,-.55);ax.set_xlim(0,100)
ax.set_xticks([0,50,100]);ax.set_xlabel('Percent',labelpad=3)
cv.text('n / N',178,142,ha='right',color=C['mid'])
for i,r in enumerate(counts.itertuples()):
    col=C['mid'] if i<2 else C['HBV']
    ax.barh(i,r.percent,height=.48,color=col,lw=0)
    ax.text(r.percent+2,i,f'{r.percent:.1f}',va='center',ha='left',fontsize=8,color=C['ink'])
    cv.text(f'{int(r.count):,} / {int(r.total):,}',178,148+(i+.55)/7.1*48,ha='right',va='center')
qa=ns.export(cv,ROOT/'rendered','Supplemental_Fig_S6',extra={'source':'plotted_data/S6_plotted_values.csv','intervals':intervals,'all_intervals_complete':all(x['complete'] for x in intervals),'panel_A_intervals':'95% spatial-null reference ranges','panel_B_intervals':'95% bootstrap confidence intervals','panel_A_MDE':'Stored approximate fold threshold at two-sided alpha=.05 and 80% power; not an effect bound','new_statistical_tests':False})
print({k:qa[k] for k in ('text_overlaps','text_outside_canvas','all_intervals_complete')})
