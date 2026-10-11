"""S7: linked matching-scheme facets, displacements and descriptive reductions."""
from pathlib import Path
import pandas as pd
from matplotlib.lines import Line2D
import nstyle as ns
from nstyle import C

ROOT=Path(__file__).resolve().parent
d=pd.read_csv(ROOT/'plotted_data/S7_plotted_values.csv')
pa=d[d.panel=='a'].reset_index(drop=True)
pb=d[(d.panel=='b')&d.numeral.notna()].copy()
pc=d[d.panel=='c'].sort_values('panel_C_display_order').reset_index(drop=True)
pd2=d[d.panel=='d']
labels=['Chromosome','Replication\ntiming','Gene\ndensity','Fragile-site\nstatus','18 strata','142 strata','Joint +\ndensity']
cv=ns.Canvas(180,258)
cv.text('Supplemental Figure S7',178,1,ha='right',fontsize=9,fontweight='bold')
cv.letter('A',2,10,'HBV − HPV contrast after genomic matching')
ax=cv.ax(21,25,152,43)
ns.clean(ax);ns.href(ax,0,ls='--')
ax.set_xlim(-.5,6.5);ax.set_ylim(0,2.75);ax.set_yticks([0,1,2])
ax.set_xticks(range(7),labels);ax.tick_params(axis='x',length=0)
ax.set_ylabel('Score difference (pp)',labelpad=3)
intervals=[]
for i,r in enumerate(pa.itertuples()):
    ax.fill_between([i-.16,i+.16],r.lower_pp,r.upper_pp,facecolor=ns.tint(C['ink'],.87),edgecolor=ns.tint(C['ink'],.3),lw=.8,zorder=2)
    ax.plot(i,r.excess_pp,'o',color=C['ink'],ms=4.5,zorder=3)
    x=21+(i+.5)/7*152
    cv.text(f'{r.mde_pp:.3f}',x,86,ha='center')
    cv.text(f'{r.emp_p_two_sided:.4f}',x,93,ha='center')
    intervals.append(dict(panel='A',scheme=r.rung,lower=float(r.lower_pp),upper=float(r.upper_pp),limits=[0,2.75],complete=bool(0<=r.lower_pp<=r.upper_pp<=2.75)))
cv.text('MDE (pp)',2,86,color=C['mid']);cv.text(ns.IT['P'],2,93,color=C['mid'])

cv.letter('B',2,106,'Virus-specific displacements under the same schemes')
cv.text('HBV',21,115,color=C['HBV'],fontweight='bold')
cv.text('HPV',35,115,color=ns.shade(C['HPV'],.15),fontweight='bold')
cv.fig.legend(handles=[Line2D([],[],marker='_',ms=10,ls='none',mec=C['ink'],mew=1.1,label='Approx. MDE'),Line2D([],[],marker='o',ms=4.4,ls='none',mec=C['mid'],mfc=C['mid'],label='At or above MDE'),Line2D([],[],marker='o',ms=4.4,ls='none',mec=C['mid'],mfc='white',label='Below MDE')],loc='upper right',bbox_to_anchor=(173/180,1-115/258),ncol=3,handlelength=.9,columnspacing=1.0)
ax=cv.ax(21,124,152,35)
ns.clean(ax);ns.href(ax,0,lw=.75,color=C['mid'])
ax.set_xlim(-.5,6.5);ax.set_ylim(-1.23,1.05);ax.set_yticks([-1,-.5,0,.5,1])
ax.set_xticks(range(7),labels);ax.tick_params(axis='x',length=0)
ax.set_ylabel('Score displacement (pp)',labelpad=3)
for v,off in [('HBV',-.10),('HPV',.10)]:
    sub=pb[pb.virus==v].set_index('rung')
    for i,scheme in enumerate(pa.rung):
        r=sub.loc[scheme];val=r.excess_pp
        ax.bar(i+off,val,width=.20,color=ns.tint(C[v],.75),edgecolor=C[v],lw=.65,zorder=2)
        threshold=r.approx_mde_signed_position_pp
        ax.plot([i+off-.16,i+off+.16],[threshold,threshold],color=C[v],lw=1.15,solid_capstyle='butt',zorder=4)
        reached=str(r.clears_mde).lower()=='true'
        ax.plot(i+off,val,marker=ns.MARK[v],ms=4.5,mfc=C[v] if reached else 'white',mec=C[v],mew=.9,ls='none',zorder=5)

cv.letter('C',2,183,'Reduction in the contrast')
short={'1_replication_timing':'Replication timing','2_gene_density':'Gene density','3_fragile_site':'Fragile-site status','4_stratum_18way':'18 strata','5_joint_142way':'142 strata','6_joint_plus_density':'Joint + density','arithmetic sum of the three singles':'Sum of single factors','breakage-density step':'Density increment'}
ax=cv.ax(43,200,52,43)
ns.clean(ax,left=False);ax.set_xlim(0,42);ax.set_ylim(9.05,-.55)
ax.set_xticks([0,20,40]);ax.set_yticks(pc.panel_C_position,[short[r] for r in pc.rung]);ax.set_xlabel('Contrast reduction (%)',labelpad=3)
for i,r in enumerate(pc.itertuples()):
    derived=r.quantity_type=='arithmetic contrast'
    ypos=r.panel_C_position
    ax.barh(ypos,r.attenuation_pct,height=.52,color='white' if derived else C['mid'],edgecolor=C['mid'],hatch='///' if derived else None,lw=.75)
    ax.text(r.attenuation_pct+1.2,ypos,f'{r.attenuation_pct:.1f}',ha='left',va='center',fontsize=8)

ax.axhline(5.7,color=C['rule'],lw=.5,zorder=.5)
cv.rule(3,96,200+(5.7+.55)/9.6*43)
cv.text('Arithmetic',3,200+(6.4+.55)/9.6*43,va='center',color=C['mid'])

cv.letter('D',108,183,'Permutation streams')
for x,label in [(130,'Draws'),(155,'Contrast (pp)'),(178,'MDE (pp)')]:cv.text(label,x,194,ha='right',color=C['mid'])
cv.rule(111,178,199)
for j,(scheme,title) in enumerate([('5_joint_142way','Joint matching'),('6_joint_plus_density','Joint + density')]):
    r=pd2[pd2.rung==scheme].iloc[0];y=203+j*27
    cv.text(title,111,y,fontweight='bold')
    for k,(draws,est,mde) in enumerate([(10000,r.excess_pp,r.mde_pp),(20000,r.excess_20k_pp,r.mde_excess_20k_pp)]):
        yy=y+8+k*6.5
        cv.text(f'{draws:,}',130,yy,ha='right')
        cv.text(f'{est:.3f}',155,yy,ha='right')
        cv.text(f'{mde:.3f}',178,yy,ha='right')
qa=ns.export(cv,ROOT/'rendered','Supplemental_Fig_S7',extra={'source':'plotted_data/S7_plotted_values.csv','intervals':intervals,'all_intervals_complete':all(x['complete'] for x in intervals),'panel_A_intervals':'Observed contrast ± 1.95996398454 null SD, multiplied by 100; spatial-null normal approximation, not bootstrap CI','panel_B_threshold_encoding':'Short horizontal ticks show reconstructed approximate signed MDE; filled/open endpoint symbols retain original threshold classification; original empty MDE fields remain unchanged','panel_B_approximate_thresholds_pp':pb[['virus','rung','approx_mde_signed_position_pp']].to_dict(orient='records'),'panel_B_threshold_limits':[-1.23,1.05],'panel_B_all_thresholds_in_axes':bool(pb.approx_mde_signed_position_pp.between(-1.23,1.05).all()),'panel_C_order':pc.rung.tolist(),'panel_C_derived_quantities':['Sum of single factors','Density increment'],'panel_D_units':'Percentage points for contrast and MDE; no intervals are implied','new_statistical_tests':False})
print({k:qa[k] for k in ('text_overlaps','text_outside_canvas','all_intervals_complete')})
