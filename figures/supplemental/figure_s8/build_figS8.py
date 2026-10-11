"""Supplemental Figure S8: HTLV-1 catalogue, class estimates and robustness."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,Rectangle

import nstyle as ns
from nstyle import C,FS
HERE=Path(__file__).resolve().parent;SRC=HERE/"plotted_data"/"S8_source_data.csv";OUT=HERE/"rendered"
d=pd.read_csv(SRC);aud=d[d.panel=="a"].reset_index(drop=True);stud=d[(d.panel=="b")&d.pmid.notna()];bal=d[(d.panel=="b")&d.covariate.notna()];atl=d[d.panel=="c"];loso=d[d.panel=="d"]
arms=[("HBV","HBV","HBV"),("HPV","HPV","HPV"),("HTLV_all","HTLV-1 pooled","HTLV-1"),("HTLV_ATL","ATL pooled","HTLV-1"),("HTLV_ATL_s35654946","ATL single study","HTLV-1"),("HTLV_nonmalignant","Non-malignant pooled","HTLV-1"),("HTLV_nonmalig_s35654946","Non-malignant single study","HTLV-1"),("HTLV_unknown","Unknown clinical state","HTLV-1")]
classes=["DEL","DUP","INV","TRA"];ests=["bc_uniqSV_blockNull","f_clusterSignTest"]
cv=ns.Canvas(180,279)
cv.text("Supplemental Figure S8",178,1,ha="right",fontsize=9,fontweight="bold")
cv.letter("A",2,10,"HTLV-1 catalogue assembly")
ax=cv.ax(57,24,29,32);ns.clean(ax,left=False);ax.set_ylim(5.5,-.5);ax.set_yticks([]);ax.set_xlim(0,60000);ax.set_xticks([0,30000,60000]);ax.set_xticklabels(["0","30","60"]);ax.set_xlabel("Records (thousands)",labelpad=3)
names=["VISDB records","Non-cell-line records","Autosomes and X","Valid positions","Mapped to GRCh37","Unique coordinates"]
for i,r in aud.iterrows():
    yy=24+(i+.5)/6*32
    ax.barh(i,r.n_records,color=C["HTLV-1"],height=.5)
    cv.text(names[i],3,yy,va="center");cv.text(f"{int(r.n_records):,}",54,yy,ha="right",va="center")
cv.letter("B",92,10,"Study contributions")
cv.text("Study (PMID)",113,19,ha="right",color=C["mid"])
ax=cv.ax(118,24,59,32);ns.clean(ax,left=False);pmids=sorted(stud.pmid.unique());ax.set_ylim(6.5,-.5);ax.set_yticks([]);ax.set_xlim(0,16000);ax.set_xticks([0,5000,10000,15000]);ax.set_xticklabels(["0","5","10","15"]);ax.set_xlabel("Coordinates (thousands)",labelpad=3)
clinical=[("HTLV_ATL","Leukaemia",C["HTLV-1"]),("HTLV_nonmalignant","Non-malignant","#789486"),("HTLV_unknown","Unknown","#C4C7C4")]
for i,pmid in enumerate(pmids):
    left=0
    for key,lab,col in clinical:
        n=float(stud[(stud.pmid==pmid)&(stud.arm==key)].n_coords.iloc[0]);ax.barh(i,n,left=left,height=.56,color=col,ec="white",lw=.35);left+=n
    cv.text(str(int(pmid)),113,24+(i+.5)/7*32,ha="right",va="center")
cv.fig.legend(handles=[Patch(fc=col,label=lab) for key,lab,col in clinical],loc="upper right",bbox_to_anchor=(.99,1-65/279),ncol=3,handlelength=.8,columnspacing=1)
cv.letter("C",2,76,"Class enrichment by clinical group")

for j,cl in enumerate(classes):
    ax=cv.ax(49+j*32.4,89,28,50);ns.clean(ax,left=False);ns.vref(ax,1)
    ax.axhline(3.665,color=C['rule'],lw=.5,zorder=0)
    for i,(key,lab,v) in enumerate(arms):
        for k,est in enumerate(ests):
            r=atl[(atl.arm==key)&(atl.svclass==cl)&(atl.estimator==est)].iloc[0];yy=i*2.2+k*.73
            m="x" if not np.isfinite(r.mde_fold) else ns.MARK[v]
            ax.plot(r.fold,yy,ls="none",marker=m,ms=3.1,mfc=C[v] if bool(r.clears_mde) else "white",mec=C[v],mew=.8)
            if v=='HTLV-1' and bool(r.survives_bonferroni): ax.text(r.fold+.013,yy,"*",va="center",fontsize=8,fontweight="bold")
        if j==0: cv.text(lab,3,89+(i*2.2+.365+.6)/17.13*50,va="center",color=C[v],fontweight="bold" if i<2 else "normal")
    ax.set_ylim(16.53,-.6);ax.set_xlim(.97,1.205);ax.set_yticks([]);ax.set_xticks([1,1.1,1.2]);ax.set_title(cl,fontweight="bold",pad=3)
cv.text("Observed / matched null",112,144,ha="center")
cv.text("Filled: deviation ≥ MDE    Open: deviation < MDE    Cross: MDE unavailable",3,152,color=C["mid"])
cv.text("Upper point: count    Lower point: cluster sign",3,157,color=C["mid"])
cv.text("* Adjusted P < 0.05: 48 HTLV-1 tests",178,157,ha="right",color=C["mid"])
cv.letter("D",2,165,"Leave-one-study-out sensitivity")
refarms=[(k,lab) for k,lab,v in arms if v=="HTLV-1" and loso[(loso.arm==k)&loso.fold.notna()].shape[0]>0]
for j,cl in enumerate(classes):
    ax=cv.ax(49+j*32.4,179,28,26);ns.clean(ax,left=False);ns.vref(ax,1)
    for i,(key,lab) in enumerate(refarms):
        g=loso[(loso.arm==key)&(loso.svclass==cl)&loso.fold.notna()]
        lo,hi=g.fold.min(),g.fold.max();baseline=float(g.all_studies_fold.iloc[0])
        ax.plot(g.fold,np.repeat(i,len(g)),ls="none",marker="o",ms=2,color=C["soft"],zorder=2);ax.plot([lo,hi],[i,i],color=C["HTLV-1"],lw=1.0);ax.plot([lo,hi],[i,i],ls="none",marker="|",ms=4,color=C["HTLV-1"])
        ax.plot(baseline,i,ls="none",marker="^",ms=4,color=C["HTLV-1"])
        if j==0: cv.text(lab,3,179+(i+.5)/4*26,va="center")
    ax.set_ylim(3.5,-.5);ax.set_yticks([]);ax.set_xlim(.98,1.056);ax.set_xticks([.99,1.02,1.05]);ax.set_title(cl,fontweight="bold",pad=3)
cv.text("Observed / matched null",112,210,ha="center")
cv.letter("E",2,220,"Covariate balance")
cv.text("Cross: before matching; open marker: after",178,225.5,ha="right",color=C["mid"])
covs=[("rt_consensus","Replication\ntiming",.60),("grad","Local\ngradient",.2),("gene_density","Gene\ndensity",.2),("in_cfs","Fragile\nsite",.2),("mapab","Mappability",.2)]
for j,(cov,lab,hi) in enumerate(covs):
    ax=cv.ax(49+j*26,239,23,32);ns.clean(ax,left=False);ax.axvline(.1,color=C["soft"],ls="--",lw=.6)
    for i,(key,labrow,v) in enumerate(arms+[("EBV","EBV","EBV")]):
        g=bal[(bal.arm==key)&(bal.covariate==cov)]
        if key=="EBV" and cov=="mapab": g=bal[(bal.arm==key)&(bal.covariate=="mappability")]
        if len(g):
            r=g.iloc[0];after=abs(r.smd_matched)
            if np.isfinite(r.smd_unmatched):
                before=abs(r.smd_unmatched);ax.plot([before,after],[i,i],color=C["rule"],lw=.7);ax.plot(before,i,ls="none",marker="x",ms=2.7,mew=.8,color=C["soft"])
            ax.plot(after,i,ls="none",marker=ns.MARK[v],ms=2.9,mfc="white",mec=C[v],mew=.8)
            if key=='EBV':
                ax.annotate(f'{after:.4f}',(after,i),xytext=(4,0),textcoords='offset points',ha='left',va='center',fontsize=8,color=C['EBV'],clip_on=False)
        else:
            # Show missingness in a separate categorical strip.
            na=Rectangle((0,i-.48),1,.96,transform=ax.get_yaxis_transform(),facecolor='#F1F3F2',edgecolor='none',zorder=5)
            na.set_gid(f'E_unavailable_{cov}');ax.add_patch(na)
            ax.text(.5,i,'N/A',transform=ax.get_yaxis_transform(),ha='center',va='center',color=C['mid'],fontsize=8,zorder=6)
        if j==0: cv.text(labrow,3,239+(i+.5)/9*32,va="center",color=C[v] if v!="EBV" else C["mid"])
    ax.set_ylim(8.5,-.5);ax.set_xlim(-.02,hi+.02);ax.set_yticks([]);ax.set_xticks([0,.3,.6] if hi>.3 else [0,.1,.2]);ax.set_title(lab,fontweight="bold",pad=3,fontsize=8)
cv.text("Absolute standardized mean difference",112,276,ha="center")
qa=ns.export(cv,OUT,"Supplemental_Fig_S8",extra={"panel_sources":[SRC.name],"atlas_cells":64,"LOSO_refits_retained_in_ranges":int(loso.fold.notna().sum()),"LOSO_ranges":"min–max over all study exclusions; triangles are all-studies baseline","MDE_fill":"absolute deviation at least approximate MDE, not significance","representation":"C/D supporting interval-and-endpoint","new_statistical_tests":False})
print(qa["text_overlaps"],qa["text_outside_canvas"])

from package_utils import finalize
assert len(atl)==64 and int(loso.fold.notna().sum())==64
assert atl.fold.min()>.97 and atl.fold.max()<1.205
assert loso.fold.dropna().between(.98,1.056).all()
family=atl[atl.arm.str.startswith('HTLV')]
assert len(family)==48 and not atl[atl.arm.isin(['HBV','HPV'])].survives_bonferroni.any()
ebv=bal[bal.arm=='EBV']
assert len(ebv)==1 and abs(float(ebv.smd_matched.iloc[0])-.1264)<1e-12
finalize(HERE,"Supplemental_Fig_S8",qa,dict(atlas_estimates=64,undefined_atlas_MDE=int(atl.mde_fold.isna().sum()),Bonferroni_family_size=len(family),HBV_HPV_outside_Bonferroni_family=True,EBV_missing_covariate_strips=4,EBV_postmatching_mappability=float(ebv.smd_matched.iloc[0]),LOSO_refits=64,LOSO_missing_rows=int(loso.fold.isna().sum()),LOSO_spans=16,covariate_rows=len(bal),atlas_all_inside_axes=True,LOSO_all_inside_axes=True,new_statistical_tests=False))
