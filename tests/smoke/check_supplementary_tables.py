#!/usr/bin/env python3
"""Read-only publication workbook/reference checks using the standard library."""
import csv
from pathlib import Path
import re
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL = '{http://schemas.openxmlformats.org/package/2006/relationships}'
RID = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'


def workbook(path):
    """Extract typed cell values without opening Excel or evaluating formulas."""
    with zipfile.ZipFile(path) as z:
        strings = []
        if 'xl/sharedStrings.xml' in z.namelist():
            strings = [''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        rels = {t.get('Id'): t.get('Target') for t in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        result = {}
        for sheet in ET.fromstring(z.read('xl/workbook.xml')).find('s:sheets', NS):
            target = rels[sheet.get(RID)]
            target = target.lstrip('/') if target.startswith('/') else 'xl/' + target
            cells = {}
            for c in ET.fromstring(z.read(target)).findall('.//s:sheetData/s:row/s:c', NS):
                if c.find('s:f', NS) is not None:
                    raise AssertionError('Unexpected formula in publication results: ' + path.name)
                v = c.find('s:v', NS)
                text = v.text if v is not None else None
                if c.get('t') == 's': value = strings[int(text)]
                elif c.get('t') == 'inlineStr': value = ''.join(c.find('s:is', NS).itertext())
                elif c.get('t') == 'b': value = text == '1'
                elif c.get('t') in ('str', 'e'): value = text
                elif text is None: value = None
                else: value = float(text)
                cells[c.get('r')] = value
            result[sheet.get('name')] = cells
        return result


def records(cells):
    def col(a):
        n = 0
        for ch in a: n = n * 26 + ord(ch) - 64
        return n
    header = {col(re.match('[A-Z]+', a)[0]): v for a, v in cells.items() if a.endswith('4') and re.fullmatch('[A-Z]+4', a)}
    rows = {}
    for address, value in cells.items():
        m = re.fullmatch('([A-Z]+)([0-9]+)', address)
        r, c = int(m[2]), col(m[1])
        if r >= 5 and c in header: rows.setdefault(r, {})[header[c]] = value
    return [r for _, r in sorted(rows.items()) if any(v is not None for v in r.values())]


def main():
    with (ROOT / 'supp/table_code_map.tsv').open(encoding='utf-8-sig', newline='') as f:
        mapping = list(csv.DictReader(f, delimiter='\t'))
    forbidden = re.compile(r'\b(?:Earlier|PASS|gate|closure|audit)\b|FINAL_PAPER|ARCHIVED|QC_PASS|RUN_PRIMARY|sha256|[FC]:[\\/]', re.I)
    books = {}
    for row in mapping:
        path = ROOT / 'supp' / row['filename']
        assert path.is_file(), path
        if '_CLEAN.xlsx' not in path.name: continue
        book = workbook(path)
        assert list(book) == row['description'].split('; '), path.name
        for name, cells in book.items():
            for address, value in cells.items():
                if isinstance(value, str): assert not forbidden.search(value), (path.name, name, address, value)
        books[row['manuscript_table_id']] = book
    assert len(books) == 4
    s2 = books['Supplementary Table S2']
    assert len(records(s2['2_Program146'])) == 146
    broad = records(s2['Primary_six_programs'])
    assert len({r['Program'] for r in broad}) == 6
    compact = records(s2['Compact_functional_modules'])
    assert len({r['Program'] for r in compact if r['Analysis_role'] == 'Spatial primary module'}) == 5
    assert all(r['Version'] == 'score-gene-excluded' for r in compact if r['Analysis_role'] == 'Spatial primary module')
    s4 = books['Supplementary Table S4']
    assert len(records(s4['4a_spatial_program_patient'])) == 110
    assert len(records(s4['4b_spatial_program_summary'])) == 5
    assert len(records(s4['7_BSW2_methods_coverage'])) == 9
    s5 = books['Supplementary Table S5']
    assert len(records(s5['2_Patient_scores'])) == 48
    genes = records(s5['3_Program146_mapping'])
    assert len(genes) == 146 and sum(r['Measurable'] is True for r in genes) == 141
    print('Publication references and four workbook structures verified. No statistical analysis executed.')


if __name__ == '__main__': main()
