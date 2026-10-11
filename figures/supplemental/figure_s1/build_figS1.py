"""Plot catalogue filtering, analysis inputs and structural-variant counting units."""
from pathlib import Path
import pandas as pd
import matplotlib as mpl
import nstyle as ns
from nstyle import C

ROOT=Path(__file__).resolve().parent
flow=pd.read_csv(ROOT/'plotted_data/S1A_catalogue_flow.csv')
inputs=pd.read_csv(ROOT/'plotted_data/S1B_analysis_inputs.csv')
units=pd.read_csv(ROOT/'plotted_data/S1C_reference_units.csv').set_index('row_id')
cv=ns.Canvas(180,216)
cv.text('Supplemental Figure S1',178,1,ha='right',fontsize=9,fontweight='bold')
geometry=[]

# A. The same count scale is retained throughout; only the counted unit changes.
cv.letter('A',2,10,'Integration catalogue filtering')
cv.text('Catalogue stage',3,21,va='center',color=C['mid'],fontsize=8.5)
key=cv.ax(62,18,73,6);key.set_axis_off()
key.legend(handles=[mpl.patches.Patch(fc=ns.tint(C[v],.38),ec=C[v],lw=.4,label=v) for v in ['HPV','HBV']],
           ncol=2,loc='center',fontsize=8.5,handlelength=1.2,columnspacing=2.5)
cv.text('EBV',153,21,ha='right',va='center',fontsize=8.5,color=C['EBV'])
ebv_key=cv.canvas_ax(140.5,18,4,6)
ebv_key.add_patch(mpl.patches.Rectangle((0,1.8),3.3,2.4,fc=C['EBV'],ec=C['EBV'],lw=.4))
cv.text('Total',178,21,ha='right',va='center',fontsize=8.5)
ax=cv.ax(62,29,73,54);ns.clean(ax,left=False)
ax.set(xlim=(0,9000),ylim=(4.5,-.5),yticks=[],xticks=[0,3000,6000,9000])
ax.set_xticklabels(['0','3,000','6,000','9,000'])
ax.set_xlabel('Catalogue size',fontsize=8.5,labelpad=3)
for i,r in flow.iterrows():
    y=29+(i+.5)/5*54
    cv.text(r.stage_label,3,y,va='center',fontsize=8.5)
    left=0
    for virus in ['HPV','HBV','EBV']:
        value=int(r[virus])
        patch=ax.barh(i,value,left=left,height=.47,fc=ns.tint(C[virus],.38) if virus!='EBV' else C[virus],
                      ec='white',lw=.45,zorder=2)[0]
        patch.set_gid(f'A_count_{i+1}_{virus}')
        if virus!='EBV':
            ax.text(left+value/2,i,f'{value:,}',va='center',ha='center',fontsize=8.5,color=C['ink'])
        geometry.append(dict(panel='A',stage_id=int(r.stage_id),virus=virus,count=value,left=left,
                             unit=r.unit,axis_mm=[62,29,73,54],xlim=[0,9000],ylim=[4.5,-.5]))
        left+=value
    assert left==int(r.total)
    cv.text(f'{int(r.EBV):,}',153,y,ha='right',va='center',fontsize=8.5,color=C['EBV'])
    cv.text(f'{int(r.total):,}',178,y,ha='right',va='center',fontsize=8.5,fontweight='bold')

# B. This maps the retained sets to analyses rather than repeating panel A bars.
cv.letter('B',2,99,'Analysis input sets')
cv.text('Coordinates',3,108,va='center',fontsize=8.5,color=C['mid'])
cv.text('Total',88,108,ha='right',va='center',fontsize=8.5,color=C['mid'])
cv.text('Used for',101,108,va='center',fontsize=8.5,color=C['mid'])
cv.rule(3,178,112,lw=.5)
for y,total,label,use in [
    (119,7846,'Matched-null eligible','Reference-model comparison'),
    (132,6957,'Amplification peaks masked','Class atlas and joint-corrected\nenrichment')]:
    r=inputs[inputs.total.eq(total)].iloc[0]
    cv.text(label,3,y,va='center',fontsize=8.5)
    cv.text(f'{int(r.total):,}',88,y,ha='right',va='center',fontsize=9,fontweight='bold')
    cv.text(use,101,y,va='center',fontsize=8.5,linespacing=1.3)
    geometry.append(dict(panel='B',total=int(r.total),unit='coordinates',display_label=label))

# C. List units separately because intervals, breakends and paired calls
# represent different counting units.
cv.letter('C',2,145,'Structural-variant reference units')
cv.text('Record set',3,155,va='center',fontsize=8.5,color=C['mid'])
for label,x in [('Donors',111),('Aliquots',140),('Calls / records',178)]:
    cv.text(label,x,155,ha='right',va='center',fontsize=8.5,color=C['mid'])
cv.rule(3,178,160,lw=.5)
rows=[
    ('paired_all','All paired calls',166),
    ('paired_excluded','Integration-excluded paired calls',175),
    ('support_all','All intervals and breakends',186),
    ('support_excluded','Integration-excluded intervals and breakends',195),
    ('large_duplication_catalogue','Large-duplication source catalogue',207)]
cv.rule(3,178,181,lw=.4)
cv.rule(3,178,201,lw=.4)
for row_id,label,y in rows:
    r=units.loc[row_id]
    cv.text(label,3,y,va='center',fontsize=8.5)
    for field,x in [('donors',111),('aliquots',140),('records',178)]:
        cv.text(f'{int(r[field]):,}',x,y,ha='right',va='center',fontsize=9,
                fontweight='bold' if field=='records' and row_id.startswith('paired') else 'normal')
    geometry.append(dict(panel='C',row_id=row_id,donors=int(r.donors),aliquots=int(r.aliquots),
                         records=int(r.records),record_unit=r.record_unit))

qa=ns.export(cv,ROOT/'rendered','Supplemental_Fig_S1',extra=dict(
    plotted_geometry=geometry,flow_stage_count=5,flow_segment_count=15,analysis_input_rows=2,reference_rows=5,
    sources=[f'plotted_data/{n}' for n in ['S1A_catalogue_flow.csv','S1B_analysis_inputs.csv','S1C_reference_units.csv']],
    no_runtime_correction=True,new_statistical_tests=False))
print({k:qa[k] for k in ['min_font_pt','text_overlaps','text_outside_canvas','label_warnings']})
