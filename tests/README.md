# Static tests

From the repository root:

```bash
python tests/smoke/check_repository.py --help
python tests/smoke/check_repository.py --rscript Rscript
python tests/smoke/check_supplementary_tables.py
```

These tests read fixed memberships, configuration and source checksums, parse Python without importing it, and parse R without evaluating it. The tests use Python's standard library and base R. They do not require Seurat, DESeq2, CNA callers or expression matrices.

The supplementary-table check verifies workbook filenames, worksheet references,
publication terminology and expected result-table sizes. It reads XLSX XML
directly and does not evaluate formulas or recompute any result.

They do not test a biological effect, generate a figure or establish full numerical equivalence. R must be supplied for a complete syntax check; an unavailable R executable produces a non-success exit code rather than a claimed success.
