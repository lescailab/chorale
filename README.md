
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chorale

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/lescailab/chorale/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lescailab/chorale/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/lescailab/chorale/graph/badge.svg)](https://app.codecov.io/gh/lescailab/chorale)
<!-- badges: end -->

CHORALE matches factors across omics modalities measured on different
individuals. It does not invent sample pairs and does not claim to
reconstruct a unique hidden biological state.

The production rule is hierarchical:

1.  ICA extracts factors separately in each modality, without phenotype
    or covariate weights.
2.  A multivariable model estimates each factor’s adjusted phenotype
    effect.
3.  Phenotype defines which cross-modal matches are eligible.
4.  Every other eligible shared covariate can refine those candidates,
    with equal weight per covariate block.
5.  Secondary agreement can never rescue absent or incompatible
    phenotype evidence.

This gives phenotype genuine priority without an arbitrary numerical
weight.

## What a result says

Each match is classified as `resolved`, `ambiguous`,
`phenotype_unsupported`, or `incompatible`. The output includes adjusted
effects, uncertainty, factor orientation, candidate sets, phenotype and
secondary margins, and reasons why covariates were excluded.

Because the cohorts are unpaired, CHORALE reports a range of compatible
cross-modal correlations. For example, `[0.35, 0.70]` says every
compatible value is positive but the data cannot select one point.
`[-0.95, 0.96]` says the relationship remains undetermined. Bootstrap
output is called a sensitivity envelope until confidence coverage is
established.

## Status

This branch is a production candidate for external review, not a stable
or CRAN release. Validation includes planted-signal simulations,
complete-null calibration, destroyed-pairing recovery, deterministic
tests, package checks, and site builds. Remaining limitations are
reported in the [methods
documentation](https://lescailab.github.io/chorale/methods).

## Installation

``` r
# install.packages("devtools")
devtools::install_github("lescailab/chorale")
```

For a reproducible development environment:

``` bash
conda env create -f environment.yml
conda activate chorale_development
```

The repository-owned environment uses R 4.5 and Node 22 and lists direct
dependencies without machine-specific build hashes.

## Minimal example

``` r
library(chorale)

sim <- chorale_simulate(
  n_modalities = 3,
  n_features = 120,
  n_shared_factors = 2,
  n_private_factors = 1,
  n_strains = 4,
  n_per_cell = 3,
  seed = 1
)

containers <- Map(chorale_load, sim$modalities, sim$col_data)
fit <- chorale_fit(containers, n_factors = c(3, 3, 3), n_init = 3,
                   n_ambiguity_boot = 19) # short documentation run
chorale_programmes(fit, significant_only = FALSE)
#>   programme n_modalities                         modalities   modality   factor
#> 4        P2            3 modality_1, modality_2, modality_3 modality_1 factor_2
#> 5        P2            3 modality_1, modality_2, modality_3 modality_2 factor_2
#> 6        P2            3 modality_1, modality_2, modality_3 modality_3 factor_2
#> 1        P1            3 modality_1, modality_2, modality_3 modality_1 factor_1
#> 2        P1            3 modality_1, modality_2, modality_3 modality_2 factor_1
#> 3        P1            3 modality_1, modality_2, modality_3 modality_3 factor_3
#> 7        P3            3 modality_1, modality_2, modality_3 modality_1 factor_3
#> 8        P3            3 modality_1, modality_2, modality_3 modality_2 factor_3
#> 9        P3            3 modality_1, modality_2, modality_3 modality_3 factor_1
#>   joint_statistic     joint_p supported     resolution_status
#> 4         15.5145 0.004975124      TRUE             ambiguous
#> 5         15.5145 0.004975124      TRUE             ambiguous
#> 6         15.5145 0.004975124      TRUE             ambiguous
#> 1          1.8293 0.398009950     FALSE phenotype_unsupported
#> 2          1.8293 0.398009950     FALSE phenotype_unsupported
#> 3          1.8293 0.398009950     FALSE phenotype_unsupported
#> 7          0.0659 1.000000000     FALSE phenotype_unsupported
#> 8          0.0659 1.000000000     FALSE phenotype_unsupported
#> 9          0.0659 1.000000000     FALSE phenotype_unsupported
#>                               secondary_evidence phenotype_column
#> 4 selected among phenotype-compatible candidates        phenotype
#> 5 selected among phenotype-compatible candidates        phenotype
#> 6 selected among phenotype-compatible candidates        phenotype
#> 1 selected among phenotype-compatible candidates        phenotype
#> 2 selected among phenotype-compatible candidates        phenotype
#> 3 selected among phenotype-compatible candidates        phenotype
#> 7                                     conflicted        phenotype
#> 8                                     conflicted        phenotype
#> 9                                     conflicted        phenotype
#>   phenotype_reference pure_features all_pure
#> 4             control          TRUE     TRUE
#> 5             control          TRUE     TRUE
#> 6             control          TRUE     TRUE
#> 1             control          TRUE     TRUE
#> 2             control          TRUE     TRUE
#> 3             control          TRUE     TRUE
#> 7             control          TRUE     TRUE
#> 8             control          TRUE     TRUE
#> 9             control          TRUE     TRUE
```

Use `chorale_control()` to set the phenotype column, reference level,
exchangeability blocks, support threshold, and ambiguity bootstrap.
There is no phenotype-weight parameter.

## Documentation

The repository vignettes are the canonical source for the [documentation
site](https://lescailab.github.io/chorale):

- [How the estimator
  works](https://lescailab.github.io/chorale/how-it-works)
- [Outputs and
  interpretation](https://lescailab.github.io/chorale/outputs)
- [Statistical methods](https://lescailab.github.io/chorale/methods)
- [Tutorial](https://lescailab.github.io/chorale/tutorial)
