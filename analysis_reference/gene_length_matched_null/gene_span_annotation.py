"""Annotate distinct GENCODE gene-body spans at genomic coordinates.

Segment tables are built once and reused through vectorized searchsorted calls.
The annotation uses overlapping gene bodies rather than nearest-gene assignment."""
from pathlib import Path
import numpy as np
import pandas as pd

SCOPES=('protein_coding','all_gene')

class GeneSpanAnnotation:
    def __init__(self,genes:pd.DataFrame,chrom_sizes:dict[str,int]):
        needed={'gene_id','gene_type','chrom','start0','end0','span_bp'}
        if not needed<=set(genes):raise ValueError(f'Missing columns: {needed-set(genes)}')
        g=genes.copy()
        if not g.gene_id.is_unique:raise ValueError('gene_id must be unique; no symbol deduplication')
        for col in ['start0','end0','span_bp']:
            if g[col].isna().any() or not np.equal(g[col],np.floor(g[col])).all():raise ValueError('Nonintegral or missing gene interval')
            g[col]=g[col].astype(np.int64)
        if not ((g.start0>=0)&(g.end0>g.start0)&g.span_bp.eq(g.end0-g.start0)).all():raise ValueError('Invalid half-open gene spans')
        self.chrom_sizes={self.normalize_chrom(c):int(n) for c,n in chrom_sizes.items()}
        if not g.chrom.isin(self.chrom_sizes).all():raise ValueError('Genes outside supported chromosome table')
        if not (g.end0<=g.chrom.map(self.chrom_sizes)).all():raise ValueError('Gene past reference chromosome end')
        self.genes=g;self._tables={}
        for scope in SCOPES:
            sub=g[g.gene_type.eq('protein_coding')] if scope=='protein_coding' else g
            self._tables[scope]={c:self._build_segments(sub[sub.chrom.eq(c)]) for c in self.chrom_sizes}

    @staticmethod
    def normalize_chrom(chrom):
        s=str(chrom);s=s[3:] if s.lower().startswith('chr') else s
        return 'chr'+s.upper()

    @classmethod
    def from_csv(cls,gene_csv=None,chrom_sizes_csv=None):
        root=Path(__file__).resolve().parent/'outputs'
        g=pd.read_csv(gene_csv or root/'gencode_v19_gene_spans.csv')
        sizes=pd.read_csv(chrom_sizes_csv or root/'chrom_sizes.csv')
        return cls(g,sizes.set_index('chrom').chrom_size_bp.to_dict())

    @staticmethod
    def _build_segments(g):
        events={};logs={}
        for r in g.itertuples():
            gid=str(r.gene_id);logs[gid]=float(np.log10(r.span_bp))
            events.setdefault(int(r.start0),[[],[]])[1].append(gid)
            events.setdefault(int(r.end0),[[],[]])[0].append(gid)
        starts=[];counts=[];med=[];ids=[];logsets=[];active=set()
        # Include 0 even when the first gene begins later; final end has an empty set.
        events.setdefault(0,[[],[]])
        for x,(remove,add) in sorted(events.items()):
            active.difference_update(remove);active.update(add)
            keys=sorted(active);vals=np.array([logs[z] for z in keys],dtype=float)
            starts.append(x);counts.append(len(keys));med.append(float(np.median(vals)) if len(vals) else np.nan)
            ids.append(';'.join(keys));logsets.append(';'.join(format(v,'.17g') for v in vals))
        return dict(starts=np.asarray(starts,np.int64),counts=np.asarray(counts,np.int32),median=np.asarray(med,float),gene_ids=np.asarray(ids,object),log10_spans=np.asarray(logsets,object))

    def annotate_arrays(self,chrom,pos0,scope='protein_coding',include_details=False,compact=False):
        """Returns aligned numpy arrays; count=-1 denotes unsupported reference domain.

        count=0 means within supported reference, outside this scope's gene bodies.
        Unsupported/negative/out-of-chromosome inputs have no length, never a false zero.
        Nonintegral or missing positions raise; callers must distinguish missing upstream.
        """
        if scope not in SCOPES:raise ValueError('Unknown scope')
        ch=np.asarray(chrom,dtype=str);p=np.asarray(pos0)
        if ch.ndim!=1 or p.ndim!=1 or len(ch)!=len(p):raise ValueError('chrom,pos0 must be equally sized vectors')
        try:numeric=p.astype(np.float64)
        except (ValueError,TypeError) as e:raise ValueError('Positions must be numeric integers') from e
        if not np.isfinite(numeric).all() or not np.equal(numeric,np.floor(numeric)).all():raise ValueError('Missing/nonintegral positions')
        if (np.abs(numeric)>2**53).any():raise ValueError('Position exceeds exactly representable integer range')
        if compact and include_details:raise ValueError('compact mode returns counts/median only')
        p=numeric.astype(np.int64);n=len(p);counts=np.full(n,-1,np.int32);median=np.full(n,np.nan)
        if not compact:norm=np.empty(n,object)
        if include_details:ids=np.full(n,'',object);logs=np.full(n,'',object)
        for raw in np.unique(ch):
            c=self.normalize_chrom(raw);idx=np.flatnonzero(ch==raw)
            if not compact:norm[idx]=c
            if c not in self.chrom_sizes:continue
            good=idx[(p[idx]>=0)&(p[idx]<self.chrom_sizes[c])]
            tab=self._tables[scope][c];j=np.searchsorted(tab['starts'],p[good],side='right')-1
            counts[good]=tab['counts'][j];median[good]=tab['median'][j]
            if include_details:ids[good]=tab['gene_ids'][j];logs[good]=tab['log10_spans'][j]
        if compact:return dict(n_overlapping_genes=counts,median_log10_span=median)
        out=dict(chrom=norm,pos0=p,n_overlapping_genes=counts,median_log10_span=median,overlap_class=np.where(counts<0,'unsupported',np.where(counts==0,'0',np.where(counts==1,'1','multi'))),annotation_state=np.where(counts<0,'unsupported_reference_domain',np.where(counts==0,'outside_scope_gene_bodies','genic')))
        if include_details:out.update(overlapping_gene_ids=ids,log10_span_set_by_gene_id=logs)
        return out

    def annotate(self,chrom,pos0,scope='protein_coding',include_details=True):
        """DataFrame convenience wrapper; original point values are never shifted."""
        d=pd.DataFrame(self.annotate_arrays(chrom,pos0,scope,include_details));d.insert(2,'scope',scope)
        return d

    def annotate_boundary_sensitivity(self,chrom,pos0,scope='protein_coding'):
        """Diagnostic annotation-only +/-1 offsets; no spatial-analysis input mutation."""
        base=self.annotate(chrom,pos0,scope,True);p=base.pos0.to_numpy(np.int64)
        for shift,label in [(-1,'minus1'),(1,'plus1')]:
            x=self.annotate_arrays(chrom,p+shift,scope,True)
            base[f'{label}_n_overlapping_genes']=x['n_overlapping_genes']
            base[f'{label}_median_log10_span']=x['median_log10_span']
            base[f'{label}_gene_ids']=x['overlapping_gene_ids']
            base[f'{label}_identity_changed']=base.overlapping_gene_ids.ne(x['overlapping_gene_ids'])|base.annotation_state.ne(x['annotation_state'])
        return base
