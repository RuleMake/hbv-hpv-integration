"""Supplemental Figure S9: boundary geometry, composition and window profiles."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.lines import Line2D

import nstyle as ns
from nstyle import C,FS
HERE=Path(__file__).resolve().parent;SRC=HERE/"plotted_data"/"S9_source_data.csv";OUT=HERE/"rendered"
s=pd.read_csv(SRC);geom=s[s.panel=="a"];comp=s[s.panel=="b"];cont=s[s.panel=="c"];prof=s[s.panel=="d"]
arms=[("HBV","HBV","HBV"),("HPV","HPV","HPV"),("HTLV_all","HTLV-1 pooled","HTLV-1"),("HTLV_ATL","ATL pooled","HTLV-1"),("HTLV_ATL_s35654946","ATL single study","HTLV-1"),("HTLV_nonmalignant","Non-malignant pooled","HTLV-1"),("HTLV_nonmalig_s35654946","Non-malignant single study","HTLV-1"),("HTLV_unknown","Unknown clinical state","HTLV-1")]
nulls=[("stratum","Matched null"),("density","Density-matched null")]
cv=ns.Canvas(180,288)
cv.text("Supplemental Figure S9",178,1,ha="right",fontsize=9,fontweight="bold")

def marker(ax,x,y,virus,filled,shape=None,ms=3.6):
    ax.plot(x,y,ls="none",marker=shape or ns.MARK[virus],ms=ms,mfc=C[virus] if filled else "white",mec=C[virus],mew=.8)

def rows(ax,n,lim,ticks):
    ns.clean(ax,left=False);ns.vref(ax,0);ax.set_ylim(n-.5,-.5);ax.set_yticks([]);ax.set_xlim(*lim);ax.set_xticks(ticks);ax.set_xticklabels(["0" if t==0 else f"{t:.2f}".replace("-","−") for t in ticks])

cv.letter("A",2,10,"Boundary displacement after endpoint exclusion")
cv.text("MDE at\n0 kb",60,16,ha="right",color=C["mid"],linespacing=1.2)
for j,rad in enumerate([0,1,10,50]):
    ax=cv.ax(65+j*29,23,26,34);rows(ax,8,(-.03,.07),[-.02,0,.06]);ax.set_title(f"{rad} kb",fontweight="bold",pad=3)
    for i,(key,label,v) in enumerate(arms):
        r=geom[(geom.arm==key)&(geom.exclusion_kb==rad)].iloc[0];marker(ax,r.displacement,i,v,bool(r.clears_mde))
        if j==0:
            cv.text(label,3,23+(i+.5)/8*34,va="center",color=C[v] if i<2 else C["ink"])
            cv.text(f"{r.mde_displacement:.4f}",60,23+(i+.5)/8*34,ha="right",va="center",color=C["mid"])
cv.text("Boundary − interior displacement",121,62,ha="center")
cv.fig.legend(handles=[Line2D([],[],ls="none",marker="o",mfc=C["mid"],mec=C["mid"],ms=4,label="Absolute deviation ≥ approximate MDE"),Line2D([],[],ls="none",marker="o",mfc="white",mec=C["mid"],ms=4,label="Absolute deviation < approximate MDE")],loc="upper center",bbox_to_anchor=(.5,1-67/288),ncol=2,handlelength=1,columnspacing=1.8)
cv.letter("B",2,78,"Local class composition at ±1 Mb")
for j,(null,label) in enumerate(nulls):
    ax=cv.ax(53+j*66,90,57,30);rows(ax,8,(-.011,.011),[-.01,0,.01]);ax.set_title(label,fontweight="bold",pad=3)
    for i,(key,lab,v) in enumerate(arms):
        r=comp[(comp.arm==key)&(comp["null"]==null)].iloc[0];marker(ax,r.excess,i,v,bool(r.clears_mde))
        if j==0: cv.text(lab,3,90+(i+.5)/8*30,va="center",color=C[v] if i<2 else C["ink"])
cv.text("Observed − null composition score",113,125,ha="center")
cv.letter("C",2,135,"Contrasts between composition displacements")
short={"HPV":"HPV","HTLV_all":"HTLV-1 pooled","HTLV_ATL":"ATL pooled","HTLV_ATL_s35654946":"ATL single","HTLV_nonmalignant":"Non-malignant pooled","HTLV_nonmalig_s35654946":"Non-malignant single","HTLV_unknown":"Unknown clinical state"}
keys=[]
for v in ["HBV","HPV"]:
    partners=(["HPV"] if v=="HBV" else [])+[k for k,lab,vir in arms if vir=="HTLV-1"]
    keys.extend([(f"{v}-vs-{key}",v,short[key]) for key in partners])
for j,(null,label) in enumerate(nulls):
    ax=cv.ax(53+j*66,148,57,48);rows(ax,13,(-.016,.016),[-.01,0,.01]);ax.set_title(label,fontweight="bold",pad=3)
    ax.axhline(6.5,color=C["rule"],lw=.6)
    for i,(key,left,right) in enumerate(keys):
        r=cont[(cont.contrast==key)&(cont["null"]==null)].iloc[0]
        ax.plot(r.excess,i,ls="none",marker="s",ms=3.1,mfc=C["ink"] if bool(r.clears_mde) else "white",mec=C["ink"],mew=.8)
        if j==0: cv.text(left+" − "+right,3,148+(i+.5)/13*48,va="center")
cv.text("Difference in composition displacements",113,201,ha="center")
cv.letter("D",2,211,"Class enrichment across neighbourhood sizes")
cv.text("Solid: observed    Dashed: upper MDE reference",178,216.5,ha="right",fontsize=8,color=C["mid"])
classes=["DEL","DUP","INV","TRA"]
lims={"DEL":(.96,1.62),"DUP":(.98,1.35),"INV":(.94,1.36),"TRA":(.93,1.28)}
ticks={"DEL":[1,1.3,1.6],"DUP":[1,1.15,1.3],"INV":[1,1.15,1.3],"TRA":[.95,1.1,1.25]}
for j,cl in enumerate(classes):
    x=14+j*43
    au=cv.ax(x,224,33,24);ns.clean(au);ns.href(au,1)
    al=cv.ax(x,260,33,18);ns.clean(al);ns.href(al,1)
    for key,v,ax in [("HBV","HBV",au),("HPV","HPV",au),("HTLV_all","HTLV-1",al)]:
        g=prof[(prof.svclass==cl)&(prof.arm==key)].sort_values("window_bp")
        xs=g.window_bp/1000
        line,=ax.plot(xs,g.fold,color=C[v],lw=.85,zorder=3);line.set_gid(f'D_observed_{key}_{cl}')
        for xx,r in zip(xs,g.itertuples()): marker(ax,xx,r.fold,v,bool(r.clears_mde),ms=2.8)
        if key=="HTLV_all":
            line,=ax.plot(xs,g.mde_fold,color=C["soft"],ls="--",lw=.7,zorder=2)
        else:
            line,=ax.plot(xs,g.mde_fold,color=ns.tint(C[v],.26),ls=(0,(3,2)),lw=.75,zorder=2)
        line.set_gid(f'D_upper_MDE_{key}_{cl}')
    for ax in [au,al]: ax.set_xscale("log");ax.set_xlim(.75,700);ax.set_xticks([1,10,100,500]);ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    au.set_ylim(*lims[cl]);au.set_yticks(ticks[cl]);au.set_xticklabels([]);au.set_title(cl,fontweight="bold",pad=3)
    al.set_ylim(.986,1.1);al.set_yticks([1,1.1]);al.set_xticklabels(["1","10","100","500"])
    if j==0:
        au.set_ylabel("Fold",labelpad=3)
        cv.text("HTLV-1",1,254,color=C["HTLV-1"],fontsize=8)
cv.text("Window half-width (kb; logarithmic scale)",96,284,ha="center")
qa=ns.export(cv,OUT,"Supplemental_Fig_S9",extra={"panel_sources":[SRC.name],"rows_used":{"A":len(geom),"B":len(comp),"C":len(cont),"D":len(prof)},"A_units":"raw boundary-minus-interior displacement; MDE at 0 kb printed","MDE_fill":"two-sided absolute-deviation power-scale comparison only, not significance","intervals":"none; all dashed D lines are upper approximate MDE power references","HBV_HPV_upper_MDE_values":64,"HTLV_upper_MDE_values":32,"representation":"A as main Figure 3; B–D supporting interval-and-endpoint","new_statistical_tests":False})
print(qa["text_overlaps"],qa["text_outside_canvas"])


from package_utils import finalize
assert [len(geom),len(comp),len(cont),len(prof)]==[32,16,26,96]
assert geom.displacement.between(-.03,.07).all()
assert comp.excess.between(-.011,.011).all()
assert cont.excess.between(-.016,.016).all()
for cl in classes:
    g=prof[(prof.svclass==cl)&prof.arm.isin(['HBV','HPV'])]
    assert g.fold.between(*lims[cl]).all(), cl
    assert g.mde_fold.between(*lims[cl]).all(), cl
h=prof[prof.arm=='HTLV_all']
assert h.fold.between(.986,1.1).all() and h.mde_fold.between(.986,1.1).all()
assert prof.window_bp.min()==1000 and prof.window_bp.max()==500000
finalize(HERE,'Supplemental_Fig_S9',qa,dict(boundary_estimates=32,composition_estimates=16,composition_contrasts=26,window_profile_estimates=96,HBV_HPV_MDE_references=64,HTLV_MDE_references=32,all_estimates_inside_axes=True,all_plotted_MDE_references_inside_axes=True,intervals='No statistical intervals saved or plotted for this figure; MDE is a power scale.',new_statistical_tests=False))
