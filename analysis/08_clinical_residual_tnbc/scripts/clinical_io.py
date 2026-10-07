"""Shared TSV utilities for the residual-TNBC workflow (Python standard library)."""
import csv, gzip
from pathlib import Path

def rows(path):
    path=Path(path)
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def save(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not data:raise ValueError('No rows to export')
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'wt',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)

def unique(data,key):
    result={r[key]:r for r in data}
    if len(result)!=len(data):raise ValueError('Duplicate '+key)
    return result
