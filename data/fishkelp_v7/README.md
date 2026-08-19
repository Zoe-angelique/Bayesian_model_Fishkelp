# FishKelp v7 data for the Blamey and Bolton fish

This directory contains a portable snapshot generated from FishKelp commit
`8c46ebb0f3d31857e468376778833e933f6e4f5d` on 19 August 2026.

## Contents

| File | Rows | Purpose |
|---|---:|---|
| `blamey_bolton_2018_fish.csv` | 30 | Table 2 fish, published and accepted names, fishery and status fields, scope flags, and record counts |
| `fishkelp_model_records_v7.csv` | 1,598 | Model-ready FishKelp records from all sources for the paper fish |
| `fishkelp_analysis_records_v7.csv` | 1,617 | Detailed pre-model records, including 15 out-of-scope records for *Chrysoblephus anglicus* |
| `fishkelp_posterior_estimates_v7.csv` | 30 | Posterior summaries; five unestimated fish retain blank estimates and explicit status fields |
| `fishkelp_species_frame_v7.csv` | 30 | Review scope, accepted names, and search names |
| `fishkelp_species_traits_v7.csv` | 30 | Trait predictors used by the v7 model where available |
| `model/tdr_model_data_v7.rds` | complete model | Full model-ready object required for hierarchical fitting |
| `model/tdr_convergence_v7.csv` | 1 | Production convergence result |
| `SOURCE_PROVENANCE.csv` | 8 | Source paths, MD5 checksums, and FishKelp commit |

The detailed analysis file contains paraphrased coding rationales. It does not
contain the source PDFs and is not a substitute for consulting them.

## Scope of the paper table

Blamey and Bolton's Table 2 covers species associated with South African kelp
forest and temperate reef ecosystems. Its columns report species, fishery, and
exploitation status. The table does not distinguish reef association from a
kelp-specific relationship. It defines the fish set here, while the FishKelp
evidence records determine the trophic and habitat rungs.

## Missing posterior estimates

The `model_status_v7` field distinguishes three reasons why a fish lacks an
estimate.

- `outside_rfbf_scope` identifies a species that occurs in the species frame
  but falls outside the review boundary.
- `not_in_rfbf_species_frame` identifies a fish listed by the paper but absent
  from the v7 review frame.
- `no_model_ready_records` identifies an in-scope fish for which no records
  occur in the fitted data.

Blank posterior fields must remain missing. They do not represent zero
evidentiary reach.

## Rebuilding the package

Run this command from the repository root.

```sh
Rscript R/package_fishkelp_v7_data.R /path/to/FishKelp
```

The full run instructions for Zoë are in
`Bayesian_framework/FishKelp_full_species_suite.qmd`.
