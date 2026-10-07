"""Portable implementation of the documented 48-patient clinical scoring formula.
This script is provided for reproduction; it was not run during release preparation.
"""
import argparse, math, statistics
from clinical_io import rows,save,unique

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['tpm','mapping','manifest','output']:p.add_argument('--'+n,required=True)
    a=p.parse_args();samples=rows(a.manifest);ids=[r['ID'] for r in samples];unique(samples,'ID')
    mapping=rows(a.mapping);unique(mapping,'Gene')
    present=[r for r in mapping if str(r['Measurable']).lower()=='true']
    if len(mapping)!=146 or len(present)!=141 or len(ids)!=48:raise ValueError('Expected 146 / 141 genes and 48 patients')
    tpm=unique(rows(a.tpm),'gene_id');zmat=[]
    for r in present:
        gid=r['GENCODE_v36_gene_id']
        if gid not in tpm:raise ValueError('Missing exact GENCODE ID '+gid)
        values=[float(tpm[gid][sid]) for sid in ids]
        if any(not math.isfinite(v) or v<0 for v in values):raise ValueError('Invalid TPM')
        logs=[math.log2(v+1) for v in values];mean=statistics.mean(logs);sd=statistics.stdev(logs)
        if sd<=0:raise ValueError('Non-variable measurable gene '+r['Gene'])
        zmat.append([(v-mean)/sd for v in logs])
    raw=[statistics.mean(g[j] for g in zmat) for j in range(48)];mean=statistics.mean(raw);sd=statistics.stdev(raw)
    save(a.output,[{'Public_ID':sid,'Program146_raw_score':raw[j],'Program146_1SD':(raw[j]-mean)/sd} for j,sid in enumerate(ids)])
if __name__=='__main__':main()
