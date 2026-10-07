# Static tests

From the repository root:

```bash
python tests/smoke/check_repository.py --help
python tests/smoke/check_repository.py --rscript Rscript
```

These tests read fixed memberships, configuration and source checksums, parse Python without importing it, and parse R without evaluating it. The tests use Python's standard library and base R. They do not require Seurat, DESeq2, CNA callers or expression matrices.

They do not test a biological effect, generate a figure or establish full numerical equivalence. R must be supplied for a complete syntax check; an unavailable R executable produces a non-success exit code rather than a claimed success.
