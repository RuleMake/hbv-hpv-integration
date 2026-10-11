"""Boundary geometry: definitions, composition and exclusion sensitivity."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.patches import Rectangle
import nstyle as ns
from nstyle import C

HERE=Path(__file__).resolve().parent
DATA=HERE/'plotted_data'
enrichment=pd.read_csv(DATA/'Figure_3A_plotted_values.csv',float_precision='round_trip').set_index('region')
sensitivity=pd.read_csv(DATA/'Figure_3B_plotted_values.csv',float_precision='round_trip').sort_values('excl_kb')
displacement=pd.read_csv(DATA/'Figure_3C_plotted_values.csv',float_precision='round_trip')
composition=pd.read_csv(DATA/'Figure_3_composition_plotted_values.csv',float_precision='round_trip')
relative_enrichment=pd.read_csv(DATA/'Figure_3E_display_values.csv',float_precision='round_trip').sort_values('excl_kb')
VIRUSES=[('HPV',C['HPV'],'D'),('HBV',C['HBV'],'o'),('HTLV-1',C['HTLV-1'],'^')]
RADII=[0,1,10,50]
XX=np.arange(4.)
cv=ns.Canvas(180,234)
cv.text('Figure 3',178,2,ha='right',fontsize=9,fontweight='bold')
geometry=[]

# A: geometric definition; positions are illustrative, not genomic coordinates.
cv.letter('A',2,2,'Boundary and interior')
sk=cv.canvas_ax(4,12,85,34)
left,right,flank=2.,72.,15.
ylong,yshort,half=12.5,28.,2.7
sk.plot([left,left,right,right],[ylong-8.4,ylong-10,ylong-10,ylong-8.4],color=C['mid'],lw=.65)
sk.text((left+right)/2,ylong-11.4,'Duplication >200 kb',ha='center',va='center',fontsize=8.5)
sk.add_patch(Rectangle((left,ylong-half),right-left,half*2,fc=ns.tint(C['HPV'],.91),ec='none'))
for x in [left,right-flank]:
    sk.add_patch(Rectangle((x,ylong-half),flank,half*2,fc=ns.tint(C['HPV'],.44),ec='none'))
for x0,x1 in [(left,left+flank),(right-flank,right)]:
    sk.annotate('',(x0,ylong-5),(x1,ylong-5),arrowprops=dict(arrowstyle='<->',lw=.6,color=C['mid'],shrinkA=0,shrinkB=0))
    sk.text((x0+x1)/2,ylong-6.9,'100 kb',ha='center',va='center',fontsize=8,color=C['mid'])
sk.plot([left,right],[ylong,ylong],color=C['HPV'],lw=.65)
for x in [left,right]: sk.plot([x,x],[ylong-half-.6,ylong+half+.6],color=C['HPV'],lw=.8)
for x,label in [(left+flank/2,'Boundary'),((left+right)/2,'Interior'),(right-flank/2,'Boundary')]:
    sk.text(x,ylong+half+2.6,label,ha='center',va='center',fontsize=8,color=C['mid'])
short_right=left+26
sk.add_patch(Rectangle((left,yshort-half),short_right-left,half*2,fc=ns.tint(C['HPV'],.44),ec='none'))
sk.plot([left,short_right],[yshort,yshort],color=C['HPV'],lw=.65)
for x in [left,short_right]: sk.plot([x,x],[yshort-half-.6,yshort+half+.6],color=C['HPV'],lw=.8)
sk.text(short_right+3,yshort-1.8,'100–200 kb',va='center',fontsize=8,color=C['mid'])
sk.text(short_right+3,yshort+2.1,'Boundary only',va='center',fontsize=8,color=C['mid'])

# B: full HPV enrichment ranges, with coordinate–interval overlap counts.
cv.letter('B',94,2,'HPV overlap enrichment')
cv.text('Overlaps',120,14,ha='right',va='center',fontsize=8,color=C['mid'])
eb=cv.ax(124,19,52,22);ns.clean(eb,left=False);ns.vref(eb,1,color=C['soft'])
eb.set(xlim=(.98,1.47),ylim=(1.7,-.7),yticks=[],xticks=[1,1.2,1.4])
eb.set_xlabel('Observed / expected overlap',fontsize=8)
for i,region in enumerate(['Boundary','Interior']):
    r=enrichment.loc[region]
    err=eb.errorbar(r.fold,i,xerr=[[r.fold-r.fold_ci_lo],[r.fold_ci_hi-r.fold]],fmt='none',ecolor=C['HPV'],
                    elinewidth=.85,capsize=2,capthick=.7,zorder=2)
    err.lines[2][0].set_gid(f'B_range_{region}')
    eb.plot(r.fold,i,marker='D',ms=4.8,mfc=C['HPV'] if i==0 else 'white',mec=C['HPV'],mew=1,ls='none',gid=f'B_point_{region}')
    eb.annotate(f'{r.fold:.3f}',(r.fold,i),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8,color=C['HPV'],fontweight='bold')
    y=19+(i+.7)/2.4*22
    cv.text(region,105,y,ha='right',va='center',fontsize=8)
    cv.text(f'{int(r.observed):,}',120,y,ha='right',va='center',fontsize=8)
    geometry.append(dict(panel='B',region=region,estimate=float(r.fold),lower=float(r.fold_ci_lo),upper=float(r.fold_ci_hi),xlim=[.98,1.47]))

# C: baseline observed and mean-null compositions, without converting means of
# fractions into ratios of mean counts. Shared zero-to-100 scale for all bars.
cv.letter('C',2,49,'Observed and expected overlap composition')
cv.text('Boundary',74.305,60,ha='center',va='center',fontsize=8,color=C['mid'])
cv.text('Interior',128.159,60,ha='center',va='center',fontsize=8,color=C['mid'])
cp=cv.ax(66.5,64,96,36);ns.clean(cp,left=False)
cp.set(xlim=(0,100),ylim=(2.6,-.6),yticks=[],xticks=[0,25,50,75,100])
cp.set_xlabel('Coordinate–interval overlaps (%)',fontsize=8.5,labelpad=3)
for i,(virus,color,marker) in enumerate(VIRUSES):
    r=composition[composition.label==virus].iloc[0]
    row_y=64+(i+.6)/3.2*36
    cv.text(virus,15.5,row_y-1.8,fontsize=8.5,fontweight='bold',color=color,va='center')
    cv.text(f'{int(r.n_total):,} overlaps',15.5,row_y+1.9,fontsize=8,color=C['mid'],va='center')
    for status,off,col in [('Observed',-.18,'observed_boundary_share'),('Null mean',.18,'null_boundary_share')]:
        y=i+off; boundary=float(r[col])*100
        token='observed' if status=='Observed' else 'null'
        boundary_patch=cp.barh(y,boundary,height=.28,color=ns.tint(color,.38),ec=color,lw=.4,zorder=2)[0]
        boundary_patch.set_gid(f'C_boundary_{virus}_{token}')
        interior_patch=cp.barh(y,100-boundary,left=boundary,height=.28,color=ns.tint(C['mid'],.89),ec=ns.tint(C['mid'],.58),lw=.4,zorder=1)[0]
        interior_patch.set_gid(f'C_interior_{virus}_{token}')
        cp.text(boundary/2,y,f'{boundary:.1f}',va='center',ha='center',fontsize=8,color=C['ink'])
        label_y=64+(y+.6)/3.2*36
        cv.text(status,63,label_y,ha='right',va='center',fontsize=8,color=C['mid'])
        geometry.append(dict(panel='C',virus=virus,status=token,boundary_pct=boundary,interior_pct=100-boundary,
                             total_overlaps=int(r.n_total) if status=='Observed' else None,xlim=[0,100]))

# D/E: retain the existing exclusion counts and folds as two aligned views.
cv.letter('D',2,113,'Retained boundary overlaps')
cv.letter('E',94,113,'Enrichment after exclusion')
dc=cv.ax(21,125,65,29);ns.clean(dc)
bars=dc.bar(XX,sensitivity.n_edge,width=.55,color=ns.tint(C['HPV'],.45),ec=C['HPV'],lw=.4,zorder=2)
for j,(patch,value) in enumerate(zip(bars,sensitivity.n_edge)):
    patch.set_gid(f'D_count_{RADII[j]}')
    dc.annotate(f'{int(value):,}',(j,value),xytext=(0,4),textcoords='offset points',ha='center',va='bottom',fontsize=8)
    geometry.append(dict(panel='D',excl_kb=RADII[j],count=float(value),ylim=[0,12600]))
dc.set(xlim=(-.62,3.62),ylim=(0,12600),xticks=XX,xticklabels=[str(x) for x in RADII],yticks=[0,5000,10000])
dc.set_yticklabels(['0','5,000','10,000']);dc.set_ylabel('Overlaps',fontsize=8.5)
ef=cv.ax(111,125,65,29);ns.clean(ef);ns.href(ef,0,color=C['soft'])
for name,field,relative_field,filled,style in [
    ('Boundary','edge_fold','boundary_relative_change_pct',True,'-'),
    ('Interior','interior_fold','interior_relative_change_pct',False,(0,(2.4,1.8)))]:
    values=relative_enrichment[relative_field].to_numpy()
    ef.plot(XX,values,color=C['HPV'],lw=1.1 if filled else .85,ls=style,gid=f'E_profile_{name}')
    ef.plot(XX,values,ls='none',marker='D',ms=4.3,mfc=C['HPV'] if filled else 'white',mec=C['HPV'],mew=.8,gid=f'E_points_{name}')
    ef.annotate(name,(3,values[-1]),xytext=(-1,7),textcoords='offset points',ha='right',fontsize=8,
                color=C['HPV'] if filled else C['mid'],fontweight='bold' if filled else 'normal')
    for j,value in enumerate(values): geometry.append(dict(panel='E',region=name,excl_kb=RADII[j],
        fold=float(sensitivity[field].iloc[j]),relative_change_pct=float(value),ylim=[-.1,.55]))
ef.set(xlim=(-.62,3.62),ylim=(-.1,.55),xticks=XX,xticklabels=[str(x) for x in RADII],yticks=[0,.2,.4])
ef.set_yticklabels(['0','0.2','0.4'])
ef.set_ylabel('Change in enrichment\nfrom 0 kb (%)',fontsize=8)
cv.text('Endpoint-exclusion radius (kb)',98.5,161.5,ha='center',va='center',fontsize=8.5)

# F: same reference ranges and discrete MDE rectangles, expressed in percentage
# points of the boundary–interior balance score, not fraction changes.
cv.letter('F',2,168,'Boundary geometry across viral catalogues')
leg=cv.ax(21,175.5,155,6);leg.set_axis_off()
leg.legend(handles=[mpl.lines.Line2D([],[],color=C['mid'],lw=.9,marker='o',ms=3.8,label='Displacement and 95% spatial range'),
                    mpl.patches.Patch(fc=ns.tint(C['mid'],.88),ec=ns.tint(C['mid'],.55),lw=.4,label='Approx. MDE (80% power)')],
           loc='center',ncol=2,fontsize=8,handlelength=1.3,columnspacing=1.3)
ylim=(-3.1,8.1)
for vi,(virus,color,marker) in enumerate(VIRUSES):
    left=22+vi*54
    cv.text(virus,left+23,187,ha='center',va='center',color=color,fontweight='bold',fontsize=9)
    ax=cv.ax(left,193,46,29);ns.clean(ax,left=vi==0);ns.href(ax,0,color=C['soft'])
    rows=displacement[displacement.label==virus].sort_values('excl_kb')
    for j,(_,r) in enumerate(rows.iterrows()):
        point=float(r.displacement)*100; low=float(r.lower)*100; high=float(r.upper)*100; mde=float(r.mde_displacement)*100
        ax.add_patch(Rectangle((j-.36,-mde),.72,2*mde,fc=ns.tint(C['mid'],.88),ec=ns.tint(C['mid'],.55),
                               lw=.4,zorder=1,gid=f'F_MDE_{virus}_{RADII[j]}'))
        err=ax.errorbar(j,point,yerr=[[point-low],[high-point]],fmt='none',ecolor=ns.tint(color,.12),
                        elinewidth=.9,capsize=2,capthick=.7,zorder=3)
        err.lines[2][0].set_gid(f'F_range_{virus}_{RADII[j]}')
        geometry.append(dict(panel='F',virus=virus,excl_kb=RADII[j],displacement_pp=point,lower_pp=low,upper_pp=high,mde_pp=mde,ylim=list(ylim)))
    values=rows.displacement.to_numpy()*100
    ax.plot(XX,values,color=color,lw=1,zorder=3,gid=f'F_profile_{virus}')
    ax.plot(XX,values,ls='none',marker=marker,ms=4.3,mfc=color,mec='white',mew=.45,zorder=4,gid=f'F_points_{virus}')
    ax.set(xlim=(-.62,3.62),ylim=ylim,xticks=XX,xticklabels=[str(x) for x in RADII],yticks=[-2,0,2,4,6,8])
    if vi==0: ax.set_ylabel('Displacement\n(percentage points)',fontsize=8.5)
    else: ax.tick_params(labelleft=False,length=0,axis='y')
cv.text('Endpoint-exclusion radius (kb)',98.5,230,ha='center',va='center',fontsize=8.5)

qa=ns.export(cv,HERE/'rendered','Yan_Fig3',extra=dict(
    plotted_geometry=geometry,interval_counts={'B':2,'F':12},counts_retained=4,enrichment_profile_values=8,
    composition_bars=6,composition_segments=12,composition_formula='Observed=n_edge/(n_edge+n_interior); mean null=observed−displacement/2',
    composition_layout='paired bars narrowed from 123 to 96 mm',composition_axis_mm=[66.5,64,96,36],
    panel_E_formula='100 * (fold_at_radius / fold_at_0_kb - 1)',panel_E_unit='percent relative change in enrichment, not percentage points',
    panel_E_baseline_folds={'Boundary':1.3188,'Interior':1.0683},
    composition_unit='Coordinate–interval overlaps; a coordinate may contribute more than once',
    composition_null_precision='Reconstructed from the saved six-decimal displacement, validated against the four-decimal stored null mean balance',
    panel_F_units='Percentage points of the boundary–interior balance score; all saved displacement/range/MDE values multiplied by 100',
    panel_F_shared_ylim=list(ylim),mde_rectangles=12,mde_interpretation='Approximate 80% power; not confidence limits',
    new_statistical_tests=False))
print({k:qa[k] for k in ['min_font_pt','text_overlaps','text_outside_canvas','label_warnings']})
