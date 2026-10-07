"""Build an exact GENCODE v36 transcript map and sum Salmon TPM/NumReads by gene ID."""
import argparse, gzip, re
from pathlib import Path
from collections import defaultdict
from clinical_io import rows,save,unique

def txmap(gtf,out):
    opener=gzip.open if str(gtf).endswith('.gz') else open
    records={}
    with opener(gtf,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)!=9 or fields[2]!='transcript':continue
            a=dict(re.findall(r'(\w+) "([^"]+)"',fields[8]))
            d={k:a[k] for k in ['transcript_id','gene_id','gene_name','gene_type']}
            tid=d['transcript_id']
            if tid in records and records[tid]!=d:raise ValueError('Conflicting exact transcript ID '+tid)
            records[tid]=d
    if not records:raise ValueError('No transcript annotations')
    save(out,[records[t] for t in sorted(records)])

def aggregate(manifest,mapping,quant_dir,out):
    m=unique(rows(mapping),'transcript_id');samples=rows(manifest);unique(samples,'ID');unique(samples,'Run')
    if len(samples)!=48:raise ValueError('Expected primary 48-run manifest')
    metadata={}
    for r in m.values():
        gid=r['gene_id'];d={'gene_id':gid,'gene_name':r['gene_name'],'gene_type':r['gene_type']}
        if gid in metadata and metadata[gid]!=d:raise ValueError('Conflicting gene metadata '+gid)
        metadata[gid]=d
    tpms={};counts={}
    for r in samples:
        sid=r['ID'];q=rows(Path(quant_dir)/sid/'quant.sf');unique(q,'Name')
        t=defaultdict(float);c=defaultdict(float)
        for v in q:
            if v['Name'] not in m:raise ValueError('Unmapped exact transcript '+v['Name'])
            gid=m[v['Name']]['gene_id'];t[gid]+=float(v['TPM']);c[gid]+=float(v['NumReads'])
        tpms[sid]=t;counts[sid]=c
    out=Path(out);out.mkdir(parents=True,exist_ok=True);ids=[r['ID'] for r in samples]
    save(out/'gene_metadata.tsv',[metadata[g] for g in sorted(metadata)])
    for name,matrix in [('gene_TPM.tsv.gz',tpms),('gene_NumReads.tsv.gz',counts)]:
        save(out/name,[dict(gene_id=g,**{sid:matrix[sid][g] for sid in ids}) for g in sorted(metadata)])

def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    t=s.add_parser('mapping');t.add_argument('--gtf',required=True);t.add_argument('--output',required=True)
    g=s.add_parser('aggregate');g.add_argument('--manifest',required=True);g.add_argument('--mapping',required=True);g.add_argument('--quant-dir',required=True);g.add_argument('--output-dir',required=True)
    a=p.parse_args()
    if a.command=='mapping':txmap(a.gtf,a.output)
    else:aggregate(a.manifest,a.mapping,a.quant_dir,a.output_dir)
if __name__=='__main__':main()
