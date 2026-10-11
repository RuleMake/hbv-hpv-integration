"""Prepare replication-timing tracks and reconstruct their cached standardization."""
from pathlib import Path
import csv, hashlib, json, shutil, sys, urllib.request

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
SOURCE = Path("upstream_project")
sys.path.insert(0, str(WORKSPACE / "analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps"))
import numpy as np
import pandas as pd


def hash_file(p, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()


def dump_csv(p, rows):
    with p.open("w", encoding="utf-8", newline="") as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    for sub in ("inputs/bigwig", "inputs/source", "inputs/official", "outputs", "runtime_wheels", "runtime_deps", "tmp"):
        (HERE/sub).mkdir(parents=True, exist_ok=True)
    sources = [(p, HERE/"inputs/bigwig"/p.name) for p in sorted((SOURCE/"data/raw/repliseq").glob("*WaveSignalRep1.bigWig"))]
    assert len(sources)==15
    for rel in ("data/interim/repliseq_bins_50kb.parquet", "src/common/repliseq.py", "src/common/config.py", "src/common/genome_features.py", "env/requirements.txt", "data/raw/okseq/PROVENANCE.md"):
        p=SOURCE/rel; sources.append((p,HERE/"inputs/source"/p.name))
    sources.append((WORKSPACE/"analysis_inputs/gene_length_matching/annotation/outputs/chrom_sizes.csv",HERE/"inputs/source/hg19_chrom_sizes.csv"))
    rows=[]
    for src,dst in sources:
        before=hash_file(src)
        if not dst.exists(): shutil.copy2(src,dst)
        copied=hash_file(dst); after=hash_file(src)
        assert before==copied==after, str(src)
        rows.append(dict(source=str(src),copy=str(dst.relative_to(HERE)),bytes=dst.stat().st_size,
                         sha256_before=before,sha256_copy=copied,sha256_after=after,
                         md5=hash_file(dst,"md5")))
    dump_csv(HERE/"input_manifest.csv",rows)
    urls={
        "ucsc_files.txt":"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/files.txt",
        "ucsc_md5sum.txt":"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/md5sum.txt",
        "ucsc_track.html":"https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeUwRepliSeq",
        "pybigtools_api.html":"https://bigtools.readthedocs.io/en/latest/api.html",
    }
    external=[]
    for name,url in urls.items():
        dst=HERE/"inputs/official"/name
        if not dst.exists():
            req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (scientific source audit)"})
            with urllib.request.urlopen(req,timeout=45) as r: dst.write_bytes(r.read())
        external.append(dict(url=url,path=str(dst.relative_to(HERE)),bytes=dst.stat().st_size,sha256=hash_file(dst)))
    dump_csv(HERE/"public_input_manifest.csv",external)
    d=pd.read_parquet(HERE/"inputs/source/repliseq_bins_50kb.parquet")
    names=[p.name.split("RepliSeq")[1].split("WaveSignal")[0] for p,_ in sources[:15]]
    z=d[names].apply(lambda s:(s-s.mean())/s.std(ddof=0))
    rebuilt=z.mean(axis=1,skipna=True)
    assert np.array_equal(rebuilt.to_numpy(),d.rt_consensus.to_numpy(),equal_nan=True)
    assert not d.duplicated(["chrom","bin_start"]).any()
    out=d[["chrom","bin_start"]].copy(); out["source_chrom_raw"]=out.chrom.astype(str)
    out["chrom"]="chr"+out.source_chrom_raw
    out["bin_end"]=out.bin_start+50000
    out["rt_consensus_legacy"]=d.rt_consensus
    out["rt_helas3_raw_cached"]=d.Helas3; out["rt_hepg2_raw_cached"]=d.Hepg2
    out["rt_helas3_z"]=z.Helas3; out["rt_hepg2_z"]=z.Hepg2
    out["n_consensus_cells"]=d[names].notna().sum(axis=1)
    out["all_three_finite"]=out[["rt_consensus_legacy","rt_helas3_z","rt_hepg2_z"]].notna().all(axis=1)
    out.to_csv(HERE/"outputs/rt_comparison_50kb.csv",index=False,float_format="%.17g")
    out.to_parquet(HERE/"outputs/rt_comparison_50kb.parquet",index=False)
    norm=[]
    for n in names:
        norm.append(dict(cell=n,n_finite=int(d[n].notna().sum()),mean_raw=float(d[n].mean()),sd_raw_ddof0=float(d[n].std(ddof=0)),
                         min_raw=float(d[n].min()),max_raw=float(d[n].max()),normalization="global finite chr1-22/X complete 50kb bins; z=(raw-mean)/sd(ddof0)"))
    dump_csv(HERE/"outputs/normalization_parameters.csv",norm)
    audit=dict(n_bins=len(d),n_cells=len(names),n_all15_finite=int(d[names].notna().all(axis=1).sum()),
               n_all3_finite=int(out.all_three_finite.sum()),consensus_exact_equal=True,max_consensus_difference=float(np.nanmax(np.abs(rebuilt-d.rt_consensus))),
               consensus_contains_target_cells=True,source_bigwig_gate="not checked by this preparation script",n_source_files=len(rows))
    (HERE/"outputs/cache_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
    print(json.dumps(audit),flush=True)


if __name__=="__main__": main()
