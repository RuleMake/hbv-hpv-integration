"""Genomic matching and donor-group comparisons, using saved results only."""
from pathlib import Path
import pandas as pd
import matplotlib as mpl
import nstyle as ns
from nstyle import C

HERE=Path(__file__).resolve().parent
P=HERE/'plotted_data'
rt=pd.read_csv(P/'Figure_5A_plotted_values.csv')
gene=pd.read_csv(P/'Figure_5B_plotted_values.csv')
tis=pd.read_csv(P/'Figure_5C_plotted_values.csv')
cv=ns.Canvas(180,183)
cv.text('Figure 5',178,2,ha='right',fontsize=9,fontweight='bold')
geom=[]

cv.letter('A',2,2,'Replication-timing references')
cv.text('HBV − HPV composition contrast',52,14,ha='center',va='center',fontsize=8.5)
tracks=[('consensus','15-cell-line\nconsensus'),('helas3','HeLa-S3'),('hepg2','HepG2')]
aa=cv.ax(21,25,62,41)
ns.clean(aa)
for i,(track,label) in enumerate(tracks):
    r=rt[rt.track==track].iloc[0]
    bars=aa.bar(i,r.estimate_pp,width=.48,color=ns.tint(C['mid'],.30),ec=C['mid'],lw=.4,zorder=2)
    bars[0].set_gid(f'A_bar_{track}')
    er=aa.errorbar(i,r.estimate_pp,yerr=[[r.estimate_pp-r.lower_pp],[r.upper_pp-r.estimate_pp]],
                  fmt='none',ecolor=C['ink'],elinewidth=.85,capsize=2.3,capthick=.7,zorder=3)
    er.lines[2][0].set_gid(f'A_range_{track}')
    geom.append(dict(panel='A',track=track,estimate=float(r.estimate_pp),lower=float(r.lower_pp),
                     upper=float(r.upper_pp),axis_min=0,axis_max=2.9))
aa.set(xlim=(-.6,2.6),ylim=(0,2.9),xticks=range(3),xticklabels=[x[1] for x in tracks],yticks=[0,1,2])
aa.set_ylabel('Contrast\n(percentage points)',fontsize=8.5)

cv.letter('B',94,2,'HPV genomic matching')
cv.text('2,568 coordinates · 1,953 clusters',142,11,ha='center',va='center',fontsize=8,color=C['mid'])
key=cv.ax(106,16,70,6);key.set_axis_off()
key.legend(handles=ns.legend_items([('Event incidence',C['HPV'],'D',True),('Cluster sign',C['HPV'],'D',False)]),
           loc='center',ncol=2,fontsize=8,handlelength=.8,columnspacing=1)
arms=[('joint_geometry','Joint\ngenomic'),('gene_status','Gene\nstatus'),('gene_length','Gene\nlength')]
ab=cv.ax(111,25,65,41);ns.clean(ab);ns.href(ab,1,color=C['soft'],lw=.65)
for i,(arm,label) in enumerate(arms):
    for j,est in enumerate(['cluster_event_incidence','calibrated_cluster_sign']):
        r=gene[(gene.arm==arm)&(gene.estimator==est)].iloc[0]
        x=i+(-.14 if j==0 else .14)
        er=ab.errorbar(x,r.estimate,yerr=[[r.estimate-r.lower],[r.upper-r.estimate]],
                      fmt='none',ecolor=C['HPV'],elinewidth=.8,capsize=2,capthick=.7,zorder=2)
        er.lines[2][0].set_gid(f'B_range_{arm}_{j}')
        ab.plot(x,r.estimate,marker='D',ms=4.5,mfc=C['HPV'] if j==0 else 'white',
                mec=C['HPV'],mew=1,ls='none',zorder=3,gid=f'B_point_{arm}_{j}')
        geom.append(dict(panel='B',arm=arm,estimator=est,estimate=float(r.estimate),lower=float(r.lower),
                         upper=float(r.upper),axis_min=.98,axis_max=1.16))
ab.set(xlim=(-.6,2.6),ylim=(.98,1.16),xticks=range(3),xticklabels=[a[1] for a in arms],yticks=[1,1.05,1.10,1.15])
ab.set_yticklabels(['1.00','1.05','1.10','1.15'])
ab.set_ylabel('Observed / matched null',fontsize=8.5)

cv.letter('C',2,84,'Large-duplication proportions by integration status')
cv.text('Donors: detected / not detected',2,99,fontsize=8,color=C['mid'],va='center')
cv.text('Mean proportion',98,99,fontsize=8.5,ha='center',va='center')
cv.text('Ratio (95% CI)',156,99,fontsize=8.5,ha='center',va='center')
ROWS=[('HBV_pooled_donor_mean','HBV',['HBV · liver']),
      ('HPV_pooled_donor_mean','HPV',['HPV · cervix and','head and neck']),
      ('HPV_common_support_histology_standardized','HPV',['HPV · squamous only','Histology standardized'])]
means=cv.ax(74,108,48,59)
ns.clean(means,left=False)
means.set(xlim=(0,50),ylim=(2.6,-.6),yticks=[],xticks=[0,25,50])
means.set_xlabel('Donor-average proportion (%)',fontsize=8.5)
ratio=cv.ax(136,108,40,59)
ns.clean(ratio,left=False);ns.vref(ratio,1,color=C['soft'])
ratio.set(xlim=(.5,1.8),ylim=(2.6,-.6),yticks=[],xticks=[.5,1,1.5])
ratio.set_xticklabels(['0.5','1.0','1.5'])
ratio.set_xlabel('Detected / not detected',fontsize=8.5)
for i,(key,virus,labels) in enumerate(ROWS):
    r=tis[tis.analysis==key].iloc[0]
    row_y=108+(i+.6)/3.2*59
    title_y=row_y-(4.3 if len(labels)==2 else 2.1)
    for j,label in enumerate(labels):
        is_group_name = j == 0 or key == 'HPV_pooled_donor_mean'
        cv.text(label,2,title_y+j*3.8,fontsize=8.5,color=C[virus] if is_group_name else C['mid'],
                fontweight='bold' if is_group_name else 'normal',va='center')
    cv.text(f'{int(r.n_detected)} / {int(r.n_comparison)}',2,
            title_y+len(labels)*3.8,fontsize=8,color=C['mid'],va='center')
    for status,off,value in [('Detected',-.16,r.mean_detected_pct),('Not detected',.16,r.mean_comparison_pct)]:
        y=i+off
        color=ns.tint(C[virus],.26) if status=='Detected' else ns.tint(C['mid'],.73)
        edge=C[virus] if status=='Detected' else ns.tint(C['mid'],.30)
        bars=means.barh(y,value,height=.23,color=color,ec=edge,lw=.4,zorder=2)
        token='detected' if status=='Detected' else 'comparison'
        bars[0].set_gid(f'C_mean_{key}_{token}')
        means.text(value+1.15,y,f'{value:.1f}',ha='left',va='center',fontsize=8,color=C['ink'])
        ym=108+(y+.6)/3.2*59
        cv.text(status,71,ym,ha='right',va='center',fontsize=8,color=C['mid'])
        geom.append(dict(panel='C_mean',analysis=key,status=token,estimate=float(value),axis_min=0,axis_max=50))
    band=ratio.hlines(i,r.lower,r.upper,color=ns.tint(C[virus],.58),lw=4.6,zorder=2)
    band.set_capstyle('butt');band.set_gid(f'C_CI_{key}')
    ratio.plot(r.estimate,i,marker=ns.MARK[virus],ms=5.2,mfc=C[virus],mec='white',mew=.5,
               ls='none',zorder=3,gid=f'C_ratio_{key}')
    ratio.text(r.estimate,i-.30,f'{r.estimate:.3f}',ha='center',va='center',fontsize=8,
               color=C[virus],fontweight='bold',bbox=dict(facecolor='white',edgecolor='none',pad=.7))
    geom.append(dict(panel='C_ratio',analysis=key,estimate=float(r.estimate),lower=float(r.lower),
                     upper=float(r.upper),axis_min=.5,axis_max=1.8))
qa=ns.export(cv,HERE/'rendered','Yan_Fig5',extra=dict(
    plotted_geometry=geom,estimates_by_panel={'A':3,'B':6,'C_means':6,'C_ratios':3},
    interval_counts={'A':3,'B':6,'C_ratios':3},
    panels_A_B_interval='95% inverse spatial reference ranges; thin capped lines',
    panel_C_interval='95% donor-bootstrap confidence intervals; thick colour bands',
    mean_intervals_displayed=False,panel_A_display_scale='100 × composition-score contrast; percentage points',
    panel_C_mean_scale='100 × donor-average or histology-standardized proportions; percent',
    new_statistical_tests=False))
print({k:qa[k] for k in ['min_font_pt','text_overlaps','text_outside_canvas','label_warnings']})
