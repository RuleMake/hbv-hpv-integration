"""Plot tissue-exchange and external-catalogue comparisons from saved results."""
from pathlib import Path
import matplotlib as mpl
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import nstyle as ns
from nstyle import C
from package_utils import finalize,sha

HERE=Path(__file__).resolve().parent
DATA=HERE/'plotted_data'
tissue=pd.read_csv(DATA/'Tissue_exchange.csv')
ext=pd.read_csv(DATA/'External_composition.csv')
vis=pd.read_csv(DATA/'VISDB_composition.csv')
cv=ns.Canvas(180,269)
cv.text('Supplemental Figure S10',178,1,ha='right',fontsize=9,fontweight='bold')
plotted=[]
all_reference=[]

def inverse_range(r):
    return r.excess-(r.null_ci_hi-r.null_mean),r.excess-(r.null_ci_lo-r.null_mean)

def save_row(panel,r,label,condition,estimate,lo,hi,p=np.nan,mde=np.nan,**kwargs):
    plotted.append(dict(panel=panel,catalogue=label,condition=condition,estimate=estimate,
        range_lo=lo,range_hi=hi,p=p,mde=mde,**kwargs))

cv.letter('A',2,10,'Tissue-specific duplication enrichment')
cv.text('Coordinates',3,20,color=C['mid'])
cv.text('SV tissue pool',29,20,color=C['mid'])
cv.text('Fold',176,20,color=C['mid'],ha='right')
a=cv.ax(77,28,87,42)
ns.clean(a,left=False);ns.vref(a,1)
order=[('HBV','1_own_tissue','Liver'),('HBV','2_cross_tissue','Cervical /\nhead and neck'),
       ('HPV','1_own_tissue','Cervical /\nhead and neck'),('HPV','2_cross_tissue','Liver')]
for i,(v,key,label) in enumerate(order):
    r=tissue[(tissue.virus==v)&(tissue.control==key)].iloc[0]
    own=key=='1_own_tissue'
    ns.dot_ci(a,i,r.fold,r.fold_ci_lo,r.fold_ci_hi,C[v],ns.MARK[v],4.8,own)
    prefix='3_downsampled_pan_' if own else '4_downsampled_to_cross_'
    rr=tissue[(tissue.virus==v)&tissue.control.str.startswith(prefix)]
    assert len(rr)==30
    # Offset only the categorical row coordinate; every saved reference fold is unchanged.
    y=i+.25+np.tile([-.035,0,.035],10)
    a.scatter(rr.fold,y,s=7,color=C['soft'],alpha=.8,linewidths=0,zorder=2)
    all_reference.extend(rr.fold.tolist())
    page_y=28+(i+.55)/4.25*42
    cv.text(v,3,page_y,va='center',color=C[v],fontweight='bold')
    cv.text(label,29,page_y,va='center',linespacing=1.1)
    cv.text(f'{r.fold:.3f}',176,page_y,ha='right',va='center',fontweight='normal')
    save_row('A',r,v,label.replace('\n',' '),r.fold,r.fold_ci_lo,r.fold_ci_hi,r.emp_p,
        n_coordinates=r.coords,n_clusters=r.n_clusters,n_donors=r.n_donors,n_events=r.n_events)
a.set(ylim=(3.7,-.55),xlim=(.78,1.94),yticks=[],xticks=[.8,1,1.2,1.4,1.6,1.8])
a.axhline(1.6,color=C['rule'],lw=.6)
cv.text('Observed / expected overlap',120.5,76,ha='center')
cv.fig.legend(handles=[mpl.lines.Line2D([],[],marker='o',ls='none',mfc=C['ink'],mec=C['ink'],label='Corresponding tissue'),
    mpl.lines.Line2D([],[],marker='o',ls='none',mfc='white',mec=C['ink'],label='Exchanged tissue'),
    mpl.lines.Line2D([],[],marker='o',ls='none',mfc=C['soft'],mec=C['soft'],ms=3,label='Pan-cancer reference')],
    loc='upper left',bbox_to_anchor=(.014,1-81/269),ncol=3,fontsize=8,columnspacing=1.5,handletextpad=.35)

cv.letter('B',2,94,'Class composition in external coordinate catalogues')
cv.text('Gray windows: ±MDE',178,101,ha='right',color=C['mid'])
b=cv.ax(22,107,151,42)
ns.clean(b);ns.href(b,0)
entries=[]
for frame,col,key,v,label in [(ext,'catalogue','CLCA','HBV','CLCA\nHBV'),
        (ext,'catalogue','Sasa2025','HPV','Sasa\nHPV'),
        (vis,'virus','VISDB_HBV','HBV','VISDB\nHBV'),
        (vis,'virus','VISDB_HPV','HPV','VISDB\nHPV')]:
    r=frame[(frame[col]==key)&(frame.arm=='stratum')&(frame.statistic=='class_axis_TRAINV_minus_DELDUP')].iloc[0]
    entries.append((r,v,label))
bxs=[0,1,2.6,3.6]
mdes=[]
for j,(x,(r,v,label)) in enumerate(zip(bxs,entries)):
    lo,hi=inverse_range(r)
    window=Rectangle((x-.19,-r.mde_excess),.38,2*r.mde_excess,
        facecolor='#E8EBEA',edgecolor='#C8CFCC',linewidth=.5,zorder=1)
    window.set_gid(f'B_MDE_window_{j}');b.add_patch(window)
    mdes.append(dict(catalogue=label.replace('\n',' '),x=x,lower=-float(r.mde_excess),upper=float(r.mde_excess)))
    ns.vdot_ci(b,x,r.excess,lo,hi,C[v],ns.MARK[v],3.5)
    b.lines[-2].set_gid(f'B_spatial_range_{j}');b.lines[-1].set_gid(f'B_estimate_{j}')
    for suffix,y in [('lo',lo),('hi',hi)]:
        cap,=b.plot([x-.05,x+.05],[y,y],color=C[v],lw=.8,zorder=3)
        cap.set_gid(f'B_spatial_cap_{suffix}_{j}')
    px=22+(x+.65)/4.9*151
    cv.text(f'{r.mde_excess:.4f}',px,166,ha='center',va='center',color=C['mid'])
    save_row('B',r,label.replace('\n',' '),'Covariate matched',r.excess,lo,hi,
        r.emp_p_two_sided,r.mde_excess,n_coordinates=r.n_coords,n_clusters=r.n_clusters)
b.set(xlim=(-.65,4.25),ylim=(-.084,.052),xticks=bxs,xticklabels=[e[2] for e in entries],
    yticks=[-.08,-.04,0,.04],ylabel='Observed − null score')
b.tick_params(axis='x',length=0,pad=5)
cv.text('MDE',3,166,va='center',color=C['mid'])

cv.letter('C',2,181,'HBV − HPV composition contrasts')
c=cv.ax(22,195,151,43)
ns.clean(c);ns.href(c,0)
method=[('stratum','Matched',C['ink'],'o',True),('density','+ density',C['mid'],'s',False)]
contrasts=[]
cx=[0,.5,2,2.5]
for i,(frame,col,key,label) in enumerate([(ext,'catalogue','CLCA-vs-Sasa2025','CLCA–Sasa'),
        (vis,'virus','VISDB_HBV-vs-VISDB_HPV','VISDB')]):
    for j,(arm,mname,color,mark,filled) in enumerate(method):
        r=frame[(frame[col]==key)&(frame.arm==arm)&frame.statistic.str.startswith('interaction')].iloc[0]
        x=cx[2*i+j];lo,hi=inverse_range(r)
        ns.vdot_ci(c,x,r.excess,lo,hi,color,mark,4.5,filled)
        px=22+(x+.5)/3.5*151
        cv.text(mname,px,252,ha='center',va='center',color=C['mid'])
        cv.text(f'{r.emp_p_two_sided:.4f}',px,258,ha='center',va='center')
        cv.text(f'{r.mde_excess:.4f}',px,264,ha='center',va='center',color=C['mid'])
        contrasts.append((r,color,mark,filled,lo,hi))
        save_row('C',r,label,'Covariate matched' if j==0 else 'Plus breakpoint density',
            r.excess,lo,hi,r.emp_p_two_sided,r.mde_excess)
c.set(xlim=(-.5,3),ylim=(-.006,.09),xticks=[.25,2.25],xticklabels=['CLCA–Sasa','VISDB'],
    yticks=[0,.04,.08],ylabel='HBV − HPV contrast')
c.tick_params(axis='x',length=0,pad=6)
cv.text('P',3,258,va='center')
cv.text('MDE',3,264,va='center',color=C['mid'])

# The inset repeats both VISDB ranges on an expanded scale.
ci=cv.ax(126,198,47,25)
ns.clean(ci);ns.href(ci,0)
for j,(r,col,mark,filled,lo,hi) in enumerate(contrasts[2:]):
    ns.vdot_ci(ci,j,r.excess,lo,hi,col,mark,4.0,filled)
ci.set(xlim=(-.6,1.6),ylim=(-.0045,.0045),xticks=[0,1],xticklabels=['Matched','+ density'],yticks=[-.004,0,.004])
ci.tick_params(axis='x',length=0,pad=3)
ci.set_title('VISDB, expanded',fontsize=8,fontweight='bold',pad=4)

pd.DataFrame(plotted).to_csv(DATA/'Figure_S10_plotted_values.csv',index=False)
qa=ns.export(cv,HERE/'rendered','Supplemental_Fig_S10',extra={
    'panel_sources':['Tissue_exchange.csv','External_composition.csv','VISDB_composition.csv'],
    'unique_spatial_ranges':12,'reference_pool_points':120,'VISDB_ranges_repeated_in_inset':2,
    'interval_type':'Stored fold reference ranges and additive inverse-null spatial ranges; not donor-sampling confidence intervals.',
    'B_encoding':'Colored point estimates and full spatial reference ranges; gray windows span negative to positive approximate MDE, a power scale, not a confidence interval.',
    'B_MDE_windows':mdes,
    'C_inset':'VISDB expanded scale duplicates the two unique VISDB ranges.',
    'new_statistical_tests':False})
print({k:qa[k] for k in ['min_font_pt','text_overlaps','text_outside_canvas','label_warnings']})
pp=pd.DataFrame(plotted)
assert list(pp.groupby('panel').size())==[4,4,4]
for panel,limits in [('A',(.78,1.94)),('B',(-.084,.052)),('C',(-.006,.09))]:
    sub=pp[pp.panel==panel]
    assert sub.range_lo.between(*limits).all() and sub.range_hi.between(*limits).all()
assert len(all_reference)==120 and min(all_reference)>.78 and max(all_reference)<1.94
assert all(-.0045<r[-2]<r[-1]<.0045 for r in contrasts[2:])
assert all(-.084<r['lower']<0<r['upper']<.052 for r in mdes)
finalize(HERE,'Supplemental_Fig_S10',qa,dict(unique_estimates=12,unique_spatial_ranges=12,
    reference_pool_points=120,VISDB_ranges_repeated_in_inset=2,all_ranges_and_points_inside_axes=True,
    full_MDE_windows=4,all_MDE_windows_inside_axes=True,B_shared_vertical_axis=[-.084,.052],
    table_P_values=4,table_MDE_values=8,new_statistical_tests=False))
mapping=pd.read_csv(HERE/'SOURCE_MAP.csv')
derived=mapping.copy()
derived['plotted_file']='plotted_data/Figure_S10_plotted_values.csv'
derived['plotted_sha256']=sha(DATA/'Figure_S10_plotted_values.csv')
derived['operation']='Selected saved estimates; B/C endpoints = excess minus reversed centred null percentiles. See inverse_range in build_figS10.py. No re-estimation.'
pd.concat([mapping,derived],ignore_index=True).to_csv(HERE/'SOURCE_MAP.csv',index=False)
