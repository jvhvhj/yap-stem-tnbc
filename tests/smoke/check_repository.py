#!/usr/bin/env python3
"""Read-only static checks. Never executes a manuscript analysis."""
import argparse
import ast
import csv
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript', default='Rscript', help='Rscript executable for parse-only checks')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[2]
    signature={}
    for name,n in [('YAP17',17),('Stem21',21),('Program146',146)]:
        data=rows(root/'config/signatures'/f'{name}.tsv')
        genes=[x['gene'] for x in data]
        if len(genes)!=n or len(set(genes))!=n:
            raise RuntimeError('Invalid membership count: '+name)
        signature[name]=set(genes)
    if signature['Program146'] & (signature['YAP17']|signature['Stem21']):
        raise RuntimeError('Score-definition genes overlap Program146')
    for p in (root/'config').glob('*.tsv'):
        if not rows(p): raise RuntimeError('Empty configuration: '+p.name)
    source=rows(root/'docs/script_provenance.tsv')
    for record in source:
        p=root/record['release_script']
        observed=hashlib.sha256(p.read_bytes()).hexdigest()
        if observed!=record['release_SHA256']:
            raise RuntimeError('Release checksum mismatch: '+record['release_script'])
    python_files=list((root/'analysis').rglob('*.py'))+list((root/'tests').rglob('*.py'))
    for p in python_files:
        ast.parse(p.read_text(encoding='utf-8-sig'), filename=str(p.relative_to(root)))
    rscript=shutil.which(args.rscript)
    if rscript is None and Path(args.rscript).is_file(): rscript=args.rscript
    if rscript is None:
        print('R parse checks NOT RUN: Rscript executable unavailable.',file=sys.stderr)
        return 2
    result=subprocess.run([rscript,'tests/smoke/parse_R_scripts.R','.'],cwd=root,capture_output=True,text=True,encoding='utf-8',errors='replace')
    if result.returncode:
        print(result.stdout,end='');print(result.stderr,end='',file=sys.stderr)
        return result.returncode
    print(result.stdout.strip())
    print(f'Python syntax: {len(python_files)} files. Release script checksums: {len(source)} unchanged.')
    print('Signatures: 17 / 21 / 146; score/program overlap: 0. Configuration files readable.')
    print('Static checks completed. No biological results reproduced and no figures generated.')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
