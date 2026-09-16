# Historical derivation sources

| File | Scientific status | Historical source timestamp (UTC) | SHA-256 |
| --- | --- | --- | --- |
| `01_wu_pseudobulk_and_reconciliation_original_2026-07-27.R` | AUTHORITATIVE_ORIGINAL_DERIVATION | 2026-07-27 13:55:03.200 | `eca08e90dd0726acae9f8fc2fd0ed03b6a3dd1fa2e939785bb3122832614f1e7` |
| `Figure3_I_effect_consistency_regeneration_2026-08-08.R` | Later regeneration/verification, not the original derivation | 2026-08-08 13:29:07.820 | `a27d5557359ed4c8857bd576c7c7ce7dc2d7a6a1ab5d7b11452cdd0aac757a4f` |

The corrected July 27 source was recovered by exact replay of its archived creation patch (13:54:23.145 UTC) and the subsequent brace correction. No scientific statement was reconstructed from a result table. Its original path was `F:/桌面/AHIPPO-YAP/0727_analysis_closure/scripts/01_wu_pseudobulk_and_reconciliation.R`. The recovered hash identifies those exact bytes; a contemporary July 27 checksum was not available for an independent comparison.

The August 8 source came from historical `tmp/Figure3_I_effect_consistency.R`. Its checksum matches the checksum recorded in the formal downstream result provenance. It consumes an established primary program and regenerates/verifies results; it is not the original discovery source.

**These files are preserved for provenance and should not be executed directly because they contain historical absolute paths and additional historical branches.** They are intentional exceptions to the portable-script path policy. No third-party package source is included. Internal execution-history transcripts are not distributed here.

The July 27 historical wrapper returned a nonzero exit code despite writing core outputs and printing a completion message. Downstream read/write records and surviving outputs close source provenance, but this is not a claim that the whole historical job or its optional enrichment branches exited cleanly. Packaging includes only parse checks and comparisons with existing data, not a new full biological run.

Use the separately documented [portable derivation](../README.md) for a future isolated reproduction. Its full execution still requires explicit input and environment validation.
