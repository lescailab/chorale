
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chorale

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/lescailab/chorale/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lescailab/chorale/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/lescailab/chorale/graph/badge.svg)](https://app.codecov.io/gh/lescailab/chorale)
<!-- badges: end -->

chorale estimates a shared latent biological state across omics
modalities measured on **disjoint sets of animals**, where no subject
contributes to more than one modality, and expresses that state in
biological terms.

The design takes its starting point from work showing that, in unpaired
multi-domain linear latent models, a shared latent block can be
identified from marginal distributions alone where the modality-specific
components are non-exchangeable and non-Gaussian (Sturma et al., NeurIPS
2023; Timilsina, Shrestha & Fu, NeurIPS 2024). Identification fails
where those modality-specific structures are made statistically alike,
so modality heterogeneity is treated as the identifying resource and
cross-modality harmonisation is excluded.

What the package implements is not those papers’ estimator. They match
recovered component distributions; chorale compares how factors respond
to the experimental design, and corroborates the result on a curated
pathway vocabulary. Its error control is therefore empirical,
established by rerunning the whole procedure under a null, and is not
inherited from those theorems.

Three properties hold throughout.

**One assignment, not many.** The correspondence between modalities is
solved once over the whole collection rather than a pair at a time, so a
programme cannot be assembled by chaining pairwise decisions through a
single mistaken link.

**Two channels that can disagree.** Factors are compared on what they do
to the design and, separately, on which curated pathways they implicate.
The factors are fitted without ever seeing those pathways, so agreement
on biology corroborates the design result instead of restating it. A
lipidome joins that comparison through its classes, in the same
vocabulary the genes use.

**Every result carries the check that would have caught it.**
Permutation calibration, annotation-matched nulls, modality shuffle, and
stability across random initialisations run alongside the estimate, and
no output is reportable without them.

## Status

Under development. The estimator, the identified sets, the controls and
the report generator all run end to end on simulated and on retrieved
data. What remains open is validation rather than implementation: the
calibration of the whole selection pipeline, the destroy-the-pairing
benchmark against data whose cross-modality truth is known, and
replication in an independent cohort. Until those are in place the
outputs say whether the approach is workable on a given pair of cohorts,
not what the biology is.

| Function                        | Status      |
|---------------------------------|-------------|
| `chorale_load()`                | Implemented |
| `chorale_map()`                 | Implemented |
| `chorale_simulate()`            | Implemented |
| `chorale_python_setup()`        | Implemented |
| `chorale_fit()`                 | Implemented |
| `chorale_programmes()`          | Implemented |
| `chorale_leave_one_out()`       | Implemented |
| `chorale_pathway_evidence()`    | Implemented |
| `chorale_metabolite_pathways()` | Implemented |
| `chorale_bound()`               | Implemented |
| `chorale_null()`                | Implemented |
| `chorale_report()`              | Implemented |

## Documentation

Full documentation, including a tutorial that runs end to end on
simulated data, is published at <https://lescailab.github.io/chorale>.
The [methods section](https://lescailab.github.io/chorale/methods) gives
every equation, where it is implemented, and the literature it rests on.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("lescailab/chorale")
```

Bioconductor dependencies resolve through `BiocManager`. Comparator
methods listed in `Suggests` are optional: their absence degrades a
comparison rather than breaking the package.

## Example

Generate disjoint cohorts over a shared feature space with known ground
truth, then load one modality into the common container:

``` r
library(chorale)

sim <- chorale_simulate(
  n_modalities = 3,
  n_features = 200,
  n_shared_factors = 4,
  n_strains = 6,
  seed = 1
)

# No animal appears in more than one modality.
length(intersect(sim$col_data[[1]]$sample_id, sim$col_data[[2]]$sample_id))
#> [1] 0

se <- chorale_load(sim$modalities[[1]], sim$col_data[[1]])
se
#> class: SummarizedExperiment 
#> dim: 200 144 
#> metadata(0):
#> assays(1): counts
#> rownames(200): modality1_feature00001 modality1_feature00002 ...
#>   modality1_feature00199 modality1_feature00200
#> rowData names(0):
#> colnames(144): modality1_sample0001 modality1_sample0002 ...
#>   modality1_sample0143 modality1_sample0144
#> colData names(9): sample_id cohort ... region batch
```

Harmonise feature identifiers to mouse Entrez, keeping one-to-many
mappings with fractional weight rather than resolving them arbitrarily:

``` r
chorale_map(c("Bdnf", "Trem2", "App"), from = "SYMBOL")
#> 
#> 'select()' returned 1:1 mapping between keys and columns
#>      id ENTREZID weight
#> 1  Bdnf    12064      1
#> 2 Trem2    83433      1
#> 3   App    11820      1
```

## Scope

The estimand is the conditional law of the latent state given the design
covariates, together with the modality-specific measurement operators.
Individual-level cross-modality inference is undefined, because each
animal is measured in exactly one modality. The joint distribution
across modalities is not identified nonparametrically, so
partial-identification bounds on the cross-modality coupling are
reported alongside any point estimate.

## References

- Sturma N, Squires C, Drton M, Uhler C. Unpaired Multi-Domain Causal
  Representation Learning. *NeurIPS* 2023.
- Timilsina S, Shrestha S, Fu X. Identifiable Shared Component Analysis
  of Unpaired Multimodal Mixtures. *NeurIPS* 2024. arXiv:2409.19422.
- Mao W, Zaslavsky E, Hartmann BM, Sealfon SC, Chikina M. Pathway-level
  information extractor (PLIER) for gene expression data. *Nat Methods*
  16:607-610, 2019.
- Taroni JN, Grayson PC, Hu Q et al. MultiPLIER: a transfer learning
  framework for transcriptomics reveals systemic features of rare
  disease. *Cell Systems* 8(5):380-394, 2019.
