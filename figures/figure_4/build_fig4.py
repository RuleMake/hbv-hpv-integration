"""Plot class enrichment, composition shifts and between-catalogue contrasts.

Class-proportion and composition-score shifts are displayed in percentage
points using the saved estimates and reference ranges.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import nstyle as ns
from nstyle import C

HERE=Path(__file__).resolve().parent
P=HERE/'plotted_data'
atlas=pd.read_csv(P/'Figure_4A_plotted_values.csv')
detail=pd.read_csv(P/'Figure_4B_detail_plotted_values.csv')
total=pd.read_csv(P/'Figure_4B_total_plotted_values.csv')
contrast=pd.read_csv(P/'Figure_4C_plotted_values.csv')
cv=ns.Canvas(180,210)
cv.text('Figure 4',178,2,ha='right',fontsize=9,fontweight='bold')
geometry=[]
classes=['DEL','DUP','INV','TRA']
viruses=['HPV','HBV','EBV']
cv.letter('A',2,2,'Structural-variant enrichment (±100 kb)')
key=cv.ax(74,9,102,6);key.set_axis_off()
key.legend(handles=ns.legend_items([('Call count',C['mid'],'o',True),('Cluster sign',C['mid'],'o',False)]),
           loc='center right',ncol=2,fontsize=8,handlelength=1,columnspacing=1.5)
for vi,virus in enumerate(viruses):
    left=20+56*vi
    ax=cv.ax(left,28,44,43)
    ns.clean(ax)
    ns.href(ax,1,color=C['soft'],lw=.65)
    limits=(.88,1.30) if virus!='EBV' else (.70,1.70)
    ax.set(xlim=(-.6,3.6),ylim=limits,xticks=range(4),xticklabels=classes)
    ax.set_yticks([.9,1,1.1,1.2,1.3] if virus!='EBV' else [.75,1,1.25,1.5,1.7])
    ax.set_yticklabels(['0.9','1.0','1.1','1.2','1.3'] if virus!='EBV' else ['0.75','1.00','1.25','1.50','1.70'])
    if vi==0: ax.set_ylabel('Observed / matched null',fontsize=8.5)
    cv.text(virus,left+22,18,ha='center',va='center',fontsize=9,fontweight='bold',color=C[virus])
    n=int(atlas[atlas.virus==virus].n_clusters.iloc[0])
    sub=f'{n:,} clusters'+(' · wider scale' if virus=='EBV' else '')
    cv.text(sub,left+22,22.7,ha='center',va='center',fontsize=8,color=C['mid'])
    for ci,cls in enumerate(classes):
        for ei,est in enumerate(['bc_uniqSV_blockNull','f_clusterSignTest']):
            r=atlas[(atlas.virus==virus)&(atlas.svclass==cls)&(atlas.estimator==est)].iloc[0]
            x=ci+(-.16 if ei==0 else .16)
            lo,hi=float(r.fold_null_interval_low),float(r.fold_null_interval_high)
            er=ax.errorbar(x,r.fold,yerr=[[r.fold-lo],[hi-r.fold]],fmt='none',
                          ecolor=C[virus],elinewidth=.75,capsize=1.6,capthick=.65,zorder=2)
            er.lines[2][0].set_gid(f'A_range_{virus}_{cls}_{ei}')
            ax.plot(x,r.fold,marker=ns.MARK[virus],ms=4.4,mfc=C[virus] if ei==0 else 'white',
                    mec=C[virus],mew=.9,ls='none',zorder=3,gid=f'A_point_{virus}_{cls}_{ei}')
            sig=float(r.emp_p_two_sided_bonferroni_24)<.05
            if sig:
                # Offset in physical points above the taller interval of this
                # class, including the companion cluster-sign range.
                class_top=float(atlas[(atlas.virus==virus)&(atlas.svclass==cls)].fold_null_interval_high.max())
                ax.annotate(f'{r.fold:.3f}*',xy=(x,class_top),xytext=(0,9),
                            textcoords='offset points',ha='center',va='bottom',fontsize=8,
                            color=C[virus],fontweight='bold',
                            gid=f'A_label_{virus}_{cls}')
            geometry.append(dict(panel='A',virus=virus,svclass=cls,estimator=est,estimate=float(r.fold),
                                 lower=lo,upper=hi,axis_min=limits[0],axis_max=limits[1],significant=sig))
cv.text('* Bonferroni-adjusted '+ns.IT['P']+' < 0.05',176,78.5,ha='right',va='center',fontsize=8,color=C['mid'])

# B: class-specific composition shifts connect the atlas with the composite axis.
cv.letter('B',2,86,'Local class composition (±1 Mb)')
key=cv.ax(120,90.5,56,6);key.set_axis_off()
key.legend(handles=[mpl.patches.Patch(fc=ns.tint(C[v],.35),ec=C[v],lw=.4,label=v) for v in ['HPV','HBV']],
           loc='center right',ncol=2,fontsize=8,handlelength=1.1,columnspacing=1)
cv.text('Class-fraction shifts',61,100,ha='center',va='center',fontsize=8.5)
cv.text('Composition-score shift',145,100,ha='center',va='center',fontsize=8.5)
cv.text('DEL + DUP',132.5,104.4,ha='center',va='center',fontsize=8,color=C['mid'])
cv.text('TRA + INV',161.5,104.4,ha='center',va='center',fontsize=8,color=C['mid'])
bd=cv.ax(22,107,78,35)
ns.clean(bd);ns.href(bd,0,color=C['soft'])
for vi,virus in enumerate(['HPV','HBV']):
    for ci,cls in enumerate(classes):
        r=detail[(detail.virus==virus)&(detail.svclass==cls)].iloc[0]
        x=ci+(-.17 if vi==0 else .17)
        bars=bd.bar(x,r.estimate_pp,width=.28,color=ns.tint(C[virus],.35),edgecolor=C[virus],linewidth=.35,zorder=2)
        bars[0].set_gid(f'B_detail_bar_{virus}_{cls}')
        er=bd.errorbar(x,r.estimate_pp,yerr=[[r.estimate_pp-r.lower_pp],[r.upper_pp-r.estimate_pp]],
                      fmt='none',ecolor=C[virus],elinewidth=.75,capsize=1.7,capthick=.65,zorder=3)
        er.lines[2][0].set_gid(f'B_detail_range_{virus}_{cls}')
        geometry.append(dict(panel='B_detail',virus=virus,svclass=cls,estimate=float(r.estimate_pp),
                             lower=float(r.lower_pp),upper=float(r.upper_pp),axis_min=-.65,axis_max=.65))
bd.set(xlim=(-.6,3.6),ylim=(-.65,.65),xticks=range(4),xticklabels=classes,yticks=[-.6,-.3,0,.3,.6])
bd.set_yticklabels(['−0.6','−0.3','0','0.3','0.6'])
bd.set_ylabel('Change (percentage points)',fontsize=8.5)
bt=cv.ax(118,107,58,35)
ns.clean(bt,left=False);ns.vref(bt,0,color=C['soft'])
for i,virus in enumerate(['HPV','HBV']):
    r=total[total.virus==virus].iloc[0]
    bar=bt.barh(i,r.estimate_pp,height=.42,color=ns.tint(C[virus],.35),ec=C[virus],lw=.35,zorder=2)
    bar[0].set_gid(f'B_total_bar_{virus}')
    er=bt.errorbar(r.estimate_pp,i,xerr=[[r.estimate_pp-r.lower_pp],[r.upper_pp-r.estimate_pp]],
                  fmt='none',ecolor=C[virus],elinewidth=.9,capsize=2,capthick=.7,zorder=3)
    er.lines[2][0].set_gid(f'B_total_range_{virus}')
    bt.text(r.estimate_pp,i-.31,f'{r.estimate_pp:+.3f}'.replace('-','−'),ha='center',va='bottom',
            fontsize=8,color=C[virus],fontweight='bold')
    geometry.append(dict(panel='B_total',virus=virus,estimate=float(r.estimate_pp),
                         lower=float(r.lower_pp),upper=float(r.upper_pp),axis_min=-1.75,axis_max=1.75))
bt.set(xlim=(-1.75,1.75),ylim=(1.6,-.7),yticks=[0,1],yticklabels=['HPV','HBV'],xticks=[-1,0,1])
bt.set_xticklabels(['−1','0','1'])
bt.set_xlabel('Shift (percentage points)',fontsize=8.5)

# C: all four saved contrast settings, zero baseline and complete reference ranges.
cv.letter('C',2,157,'HBV − HPV composition contrast (±1 Mb)')
cx=[0,1,2.6,3.6]
rows=[('paired_sample','stratum'),('paired_sample','density'),('donor_signature','stratum'),('donor_signature','density')]
ac=cv.ax(27,174,149,25)
ns.clean(ac);ns.href(ac,0,color=C['soft'])
for i,((variant,arm),x) in enumerate(zip(rows,cx)):
    r=contrast[(contrast.variant==variant)&(contrast.arm==arm)].iloc[0]
    bar=ac.bar(x,r.estimate_pp,width=.52,color=ns.tint(C['mid'],.22 if arm=='stratum' else .63),
               ec=C['mid'],lw=.4,zorder=2)
    bar[0].set_gid(f'C_bar_{variant}_{arm}')
    er=ac.errorbar(x,r.estimate_pp,yerr=[[r.estimate_pp-r.lower_pp],[r.upper_pp-r.estimate_pp]],
                  fmt='none',ecolor=C['ink'],elinewidth=.8,capsize=2.2,capthick=.7,zorder=3)
    er.lines[2][0].set_gid(f'C_range_{variant}_{arm}')
    geometry.append(dict(panel='C',variant=variant,arm=arm,estimate=float(r.estimate_pp),
                         lower=float(r.lower_pp),upper=float(r.upper_pp),axis_min=0,axis_max=2.6))
ac.set(xlim=(-.65,4.25),ylim=(0,2.6),xticks=cx,xticklabels=['Genomic matching','Density matching','Genomic matching','Density matching'],yticks=[0,1,2])
ac.set_ylabel('Contrast\n(percentage points)',fontsize=8.5)
for x,label in [(.5,'Paired calls'),(3.1,'Grouped calls per donor')]:
    xf=27+(x+.65)/4.9*149
    cv.text(label,xf,168.5,ha='center',va='center',fontsize=8.5)
qa=ns.export(cv,HERE/'rendered','Yan_Fig4',extra=dict(
    plotted_geometry=geometry,estimates_by_panel={'A':24,'B_detail':8,'B_total':2,'C':4},
    interval_type='95% inverse spatial reference ranges; not donor-sampling confidence intervals',
    significance_family_size=24,significance_field='emp_p_two_sided_bonferroni_24',
    composition_display_scale='100 × saved proportion or score displacement; percentage points',
    panel_A_background_donors=1902,panel_B_C_background_donors=1849,
    new_statistical_tests=False))
print({k:qa[k] for k in ['min_font_pt','text_overlaps','text_outside_canvas','label_warnings']})
