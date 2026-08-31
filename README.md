
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chorale

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/lescailab/chorale/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lescailab/chorale/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/lescailab/chorale/graph/badge.svg)](https://app.codecov.io/gh/lescailab/chorale)
<!-- badges: end -->

CHORALE analyses omics modalities measured on different individuals by
first placing every sample in the same space of named biological
concepts.

The modalities share no sample and need share no feature. Their
connection is a fixed vocabulary: pathways, cell-type signatures, or
other feature sets chosen before the results are inspected. Once each
modality has been scored on that vocabulary, CHORALE supports two
complementary analyses.

- **Concept evidence** estimates the adjusted phenotype effect
  separately in each modality and combines those effects by
  inverse-variance weighting.
- **Joint-state analysis** stacks all samples by their concept scores
  and fits a low-rank representation, then tests and transfers the
  resulting directions.

Related concepts can also be tested as families, while variation outside
the vocabulary is retained as free dimensions.

The shared preparation runs in a fixed order:

1.  The vocabulary is fixed and its coverage reported, per concept and
    per modality.
2.  Every modality is scored on it, without the design being consulted.
3.  What the vocabulary does not explain is kept as free dimensions,
    selected by reproducibility.
4.  The concept branch regresses each score on the design within each
    modality and combines the adjusted effects.
5.  The joint branch standardises scores within modality, optionally
    removes modality-local nuisance effects, and fits the observed
    entries of the stacked score matrix.
6.  Freedman–Lane residual permutations hold the design fixed while
    calibrating concept, family, joint-component, and free-dimension
    evidence.

## What a result says

Each concept carries its adjusted phenotype effect in every modality
that expresses it, the combined statistic, agreement diagnostics,
family-wise and false-discovery control, and the statistic that survives
adjustment for overlapping neighbours. A high heterogeneity p-value
means that excess disagreement was not detected, which with few
modalities the test has little power to do.

A joint component is a fitted direction over concepts. Its sign and
label are arbitrary, so it is interpreted through its loadings,
calibrated phenotype evidence, per-modality effects, and
leave-one-modality-out transfer. Coordinated variation no concept
explains is reported separately as free dimensions.

## Status

This branch is a production candidate for external review. Validation
includes planted-concept simulations, complete-null calibration,
recovery measured against a withheld pairing, deterministic tests,
package checks, and site builds. Conditions of applicability are stated
in the [methods
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

The package requires R 4.5. The repository-owned environment uses R 4.5
and Node 22 and lists direct dependencies without machine-specific build
hashes.

## Minimal example

``` r
library(chorale)

# Two cohorts on disjoint individuals, with one named concept planted in both.
fx <- chorale_concept_example(seed = 1)

fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
                           n_permutations = 199, n_init = 3)
fit$evidence$joint[, c("concept", "n_modalities", "joint_z",
                       "sign_agreement", "q_value")]
#>   concept n_modalities joint_z sign_agreement q_value
#> 1 planted            2 26.9901            1.0  0.0000
#> 2   quiet            2  0.2657            1.0  0.9112
#> 3   other            2  0.1165            0.5  0.9112

# The same encoding can be analysed jointly. Component count and nuisance
# adjustment are decisions of the call.
state <- chorale_joint_state(fit$encoding, n_components = 1)
chorale_joint_evidence(state, n_permutations = 199)$components
#>   component           term  share    effect        se        z p_value p_family
#> 1  joint_01 phenotype=case 0.3562 0.3604175 0.1871729 1.925586    0.07     0.07
#>   q_value significant family_significant
#> 1  0.0653       FALSE              FALSE
```

Use `chorale_control()` to set the phenotype column, reference level,
exchangeability blocks, the reproducibility threshold the free
dimensions are selected at, and the error rate concepts are called at.

## Documentation

The repository vignettes are the canonical source for the [documentation
site](https://lescailab.github.io/chorale):

- [Why CHORALE exists](https://lescailab.github.io/chorale/why)
- [How the estimator
  works](https://lescailab.github.io/chorale/how-it-works)
- [The input data](https://lescailab.github.io/chorale/input-format)
- [Tutorial](https://lescailab.github.io/chorale/tutorial)
- [Outputs and how to read
  them](https://lescailab.github.io/chorale/outputs)
- [Statistical methods](https://lescailab.github.io/chorale/methods)
- [Simulation and
  validation](https://lescailab.github.io/chorale/simulation)
