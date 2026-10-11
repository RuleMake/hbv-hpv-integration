"""Prepare Figure 2 display labels and units from saved estimates."""
from pathlib import Path
import hashlib
import json
import pandas as pd

ROOT=Path(__file__).resolve().parent
(ROOT/'review').mkdir(exist_ok=True)
(ROOT/'plotted_data').mkdir(exist_ok=True)
RAW=ROOT/'source_data/main'
OUT=ROOT/'plotted_data'
OUT.mkdir(exist_ok=True)
arms={'published':('Records','Independent'), 'recordblock':('Records','Cluster-preserving'),
      'dedup':('Loci','Independent'), 'corrected':('Loci','Cluster-preserving')}
schemes={0:('Records','Independent','None'),1:('Records','Cluster-preserving','None'),
         2:('Loci','Independent','None'),3:('Loci','Cluster-preserving','None'),
         4:('Loci','Cluster-preserving','Joint matching; peak exclusion')}

a=pd.read_csv(RAW/'Figure_2A.csv')
rows=[]
for r in a.itertuples(index=False):
    unit,null,adjust=schemes[r.rung]
    rows.append(dict(figure='Figure 2',panel='A',scheme=r.rung+1,counting_unit=unit,
        spatial_reference=null,additional_adjustment=adjust,virus=r.virus,
        observed_expected=r.fold,lower=r.lo,upper=r.hi,range_type='95% spatial reference range',
        n_coordinates=r.n_coords,window_bp=100000))
pd.DataFrame(rows).to_csv(OUT/'Figure_2A_plotted_values.csv',index=False)

for panel in ['B','C']:
    df=pd.read_csv(RAW/f'Figure_2{panel}.csv')
    rows=[]
    for r in df.itertuples(index=False):
        unit,null=arms[r.arm]
        row=dict(figure='Figure 2',panel=panel,generator=r.generator,counting_unit=unit,
            spatial_reference=null,virus=r.virus,window_bp=r.window_bp,n_simulations=r.n_rep)
        rows.append(row)
    out=pd.DataFrame(rows)
    out['rejection_rate_pct']=100*df['FPR_0.05'].to_numpy()
    if panel=='C':
        out['lower_pct']=100*df['FPR_0.05_lo'].to_numpy()
        out['upper_pct']=100*df['FPR_0.05_hi'].to_numpy()
        out['interval_type']='95% binomial Wilson confidence interval'
    else:
        out['interval_type']='None; descriptive pooled rate with reused seeds across settings'
    out.to_csv(OUT/f'Figure_2{panel}_plotted_values.csv',index=False)

d=pd.read_csv(RAW/'Figure_2D.csv')
d['figure']='Figure 2'
d['panel']='D'
d['catalogue_group']=d.catalogue.map(lambda s:'Primary' if s.startswith('Primary') else 'External HPV')
d['virus']=d.catalogue.map(lambda s:s.split()[-1])
d['catalogue_label']=d.catalogue.map({'Primary HBV':'HBV','Primary HPV':'HPV','Primary EBV':'EBV',
    'Sasa HPV':'HPV (Sasa)','Fan HPV':'HPV (Fan)'})
d=d.rename(columns={'n_coords':'n_coordinates','multiplicity':'mean_coordinates_per_cluster'})
d.to_csv(OUT/'Figure_2D_plotted_values.csv',index=False)

mapping=[]
for panel in 'ABCD':
    p=RAW/f'Figure_2{panel}.csv'
    mapping.append(dict(figure='Figure 2',panel=panel,input=str(p.relative_to(ROOT)),
        plotted_table=f'plotted_data/Figure_2{panel}_plotted_values.csv',sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        transformation='Preserve saved values; use readable method names'+('; multiply rates and bounds by 100' if panel in 'BC' else '')))
pd.DataFrame(mapping).to_csv(ROOT/'SOURCE_MAP.csv',index=False)
(ROOT/'review/source_preparation.json').write_text(json.dumps(dict(source_files=mapping,new_statistical_tests=False),indent=2)+'\n',encoding='utf-8')
print('Prepared 15 enrichment estimates, 8 pooled rates, 32 calibration intervals and 5 catalogue means.')
