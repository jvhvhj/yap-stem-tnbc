# study-defined gene definitions

| File | Members | Meaning |
| --- | ---: | --- |
| `YAP17.tsv` | 17 | Primary YAP-associated transcriptional score |
| `Stem21.tsv` | 21 | Primary stemness-associated transcriptional score |
| `Program146.tsv` | 146 | score-gene-excluded High-associated transcriptional program |

Each TSV preserves the source gene symbols and order. The two score sets jointly contain 38 distinct genes. None of these 38 genes belongs to Program146. Program146 membership matches both the existing patient-direction gene list and the gene flags in the 17,597-gene formal result table.

## Scoring

The primary component scores use the arithmetic mean of measurable study-defined genes in the appropriate gene-expression layer, not a low-dimensional integrated embedding. In Wu discovery, the component values come from the RNA/data log-normalized expression layer. The stored discovery Joint axis is the sum of the two component z scores over the fixed 10,836-cell population. Pooled tertiles define its Low, Intermediate and High states.

Different validation/spatial source scripts contain explicitly scoped normalization and transfer definitions. Their constants and scale must be recovered from those specific sources. This snapshot does not impose a newly unified scaling rule or regenerate missing constants. A rank correlation is invariant to a positive scalar multiplication, but a numeric score value or a transferred threshold is not interchangeable across scales.

## Program146 derivation

The final recorded criterion is patient-blocked High-versus-Other BH FDR < 0.05, unshrunk log2 fold change >= 0.5, and positive within-patient direction in at least 6 of 8 patients. The 38 score-defining genes are excluded before the final testable universe. The existing formal table contains 17,597 testable genes and 362 FDR-significant genes.

The corrected July 27 original derivation has been recovered and archived. The separate August 8 regeneration source matches its historically recorded checksum; it is not relabelled as the original. A [portable derivation](../../analysis/04_program146/README.md) preserves the original primary calculations. Existing-output checks and input-only preflight do not establish full numerical rerun equivalence.

`Program146_representative_genes.tsv` contains only ATF3, PPP1R15A, MCL1 and ZFP36L1, all verified members. CDKN1A, PERP and S100A10 are not Program146 members. This restricted representative list does not replace the full definition.

## Provenance and transfer rules

## Auxiliary resources

`YAP_TEAD_auxiliary9.tsv` defines the **9-gene YAP/TEAD-associated auxiliary gene set**. `epithelial_plasticity_auxiliary12.tsv` defines the **12-gene epithelial stemness/plasticity auxiliary gene set**. Both are literature-informed gene sets defined for this study and used only for marker-set/scoring-definition sensitivity. They are not independently validated signatures and do not replace YAP17 or Stem21. Representative literature supports biological rationale, not the exact complete membership of these study-defined sets.

`cell_cycle_covariates.tsv` records the exact Seurat `cc.genes.updated.2019` S (43) and G2M (54) memberships and cohort coverage. Scores use arithmetic means of measurable log-normalized expression, with no matched-control subtraction or discrete phase assignment. Missing genes are omitted, not set to zero. Wu covers S43/43 and G2M52/54 (PIMREG and JPT1 absent); ARTEMIS covers S43/43 and G2M54/54. Biological source: [Tirosh et al. 2016](https://doi.org/10.1126/science.aad0501).

`current_functional_program_membership.tsv` separates original membership, YAP17/Stem21 exclusion, primary retained membership and Program146-overlap-excluded sensitivity membership. The four Hallmarks are TNFα–NF-κB, EMT, Hypoxia and Apoptosis (MSigDB Hallmark v2024.Hs). Survival–stress and adhesion–remodeling are gene sets defined for this study, without an external signature-source claim. These broad six-program memberships are distinct from the compact study-defined spatial modules in S2 `Compact_functional_modules`. The spatial primary versions exclude YAP17/Stem21 genes but retain Program146 overlap. Spatial UPR is the study-defined eight-gene module, not the complete Hallmark UPR set.

YAP17/Stem21 membership was reconciled against the July 19 signature provenance records; Program146 was reconciled against the July 28 direction-sensitivity records and the formal gene-level result table. These dates identify the source records, not a new scientific definition.

No re-selection, fuzzy substitution or subtype-specific optimization is permitted. Missing genes remain missing. Previously approved, unambiguous aliases must be documented by the corresponding dataset's existing mapping table. Full membership files are definitions, not a claim that every external dataset measures all members.
