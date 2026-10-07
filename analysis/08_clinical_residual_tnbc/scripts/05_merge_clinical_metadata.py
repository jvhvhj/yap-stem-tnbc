"""Exact-ID clinical merge; operative sensitivity reuses scores from all 48 patients."""
import argparse
from collections import Counter
from clinical_io import rows,save,unique

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',required=True);p.add_argument('--scores',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    m=unique(rows(a.manifest),'ID');s=unique(rows(a.scores),'Public_ID')
    if set(m)!=set(s) or len(m)!=48:raise ValueError('Patient IDs do not match the primary 48')
    if Counter(r['Group'] for r in m.values())!={'D':34,'ND':14}:raise ValueError('Unexpected groups')
    data=[dict(s[sid],Outcome_group=m[sid]['Group'],D_event=int(m[sid]['Group']=='D'),RCB_class=m[sid]['RCB_class'],Neo_subtype=m[sid]['Neo_subtype'],Operative_subtype=m[sid]['Op_subtype'],Included_operative_TNBC_sensitivity=m[sid]['Op_subtype']=='HR-HER2-') for sid in s]
    if sum(r['Included_operative_TNBC_sensitivity'] for r in data)!=45:raise ValueError('Expected sensitivity n45')
    save(a.output,data)
if __name__=='__main__':main()
