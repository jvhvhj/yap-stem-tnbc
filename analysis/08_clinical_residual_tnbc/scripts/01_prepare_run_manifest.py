"""Validate the supplied public 48-run manifest and list paired-read download IDs."""
import argparse
from collections import Counter
from clinical_io import rows,save,unique

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    d=rows(a.manifest);unique(d,'ID');unique(d,'Run')
    if len(d)!=48 or Counter(r['Group'] for r in d)!={'D':34,'ND':14}:raise ValueError('Expected 48 patients, D34 / ND14')
    if any(r['LibraryLayout']!='PAIRED' or r['Neo_subtype']!='HR-HER2-' for r in d):raise ValueError('Unexpected sequencing layout or Neo subtype')
    save(a.output,[{'Public_ID':r['ID'],'Run':r['Run'],'Read1_name':r['Run']+'_1.fastq.gz','Read2_name':r['Run']+'_2.fastq.gz'} for r in d])
if __name__=='__main__':main()
