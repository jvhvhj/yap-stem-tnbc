# Fixed gene definitions

| File | Members | Meaning |
| --- | ---: | --- |
| `YAP17.tsv` | 17 | Primary YAP-associated transcriptional score |
| `Stem21.tsv` | 21 | Primary stemness-associated transcriptional score |
| `Program146.tsv` | 146 | Score-independent High-associated transcriptional program |

Each TSV preserves the source gene symbols and order. The two score sets jointly contain 38 distinct genes. None of these 38 genes belongs to Program146. Program146 membership matches both the existing patient-direction gene list and the gene flags in the 17,597-gene formal result table.

## Scoring

The primary component scores use the arithmetic mean of measurable fixed genes in the appropriate gene-expression layer, not a low-dimensional integrated embedding. In Wu discovery, the component values come from the RNA/data log-normalized expression layer. The stored discovery Joint axis is the sum of the two component z scores over the fixed 10,836-cell population. Pooled tertiles define its Low, Intermediate and High states.

Different validation/spatial source scripts contain explicitly scoped normalization and transfer definitions. Their constants and scale must be recovered from those specific sources. This snapshot does not impose a newly unified scaling rule or regenerate missing constants. A rank correlation is invariant to a positive scalar multiplication, but a numeric score value or a transferred threshold is not interchangeable across scales.

## Program146 derivation

The final recorded criterion is patient-blocked High-versus-Other BH FDR < 0.05, unshrunk log2 fold change >= 0.5, and positive within-patient direction in at least 6 of 8 patients. The 38 score-defining genes are excluded before the final testable universe. The existing formal table contains 17,597 testable genes and 362 FDR-significant genes.

The corrected July 27 original derivation has been recovered and archived. The separate August 8 regeneration source matches its historically recorded checksum; it is not relabelled as the original. A [portable derivation](../../analysis/04_program146/README.md) preserves the original primary calculations. Existing-output checks and input-only preflight do not establish full numerical rerun equivalence.

`Program146_representative_genes.tsv` contains only ATF3, PPP1R15A, MCL1 and ZFP36L1, all verified members. CDKN1A, PERP and S100A10 are not Program146 members. This restricted representative list does not replace the full definition.

## Provenance and transfer rules

YAP17/Stem21 membership was reconciled against the July 19 signature provenance records; Program146 was reconciled against the July 28 direction-sensitivity records and the formal gene-level result table. These dates identify the source records, not a new scientific definition.

No re-selection, fuzzy substitution or subtype-specific optimization is permitted. Missing genes remain missing. Previously approved, unambiguous aliases must be documented by the corresponding dataset's existing mapping table. Full membership files are definitions, not a claim that every external dataset measures all members.
