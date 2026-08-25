# chorale 0.0.0.9000

First development version. The package estimates a shared latent biological
state across omics modalities measured on disjoint sets of individuals.

## Estimator

* `chorale_concept_fit()` runs the estimator end to end: build the vocabulary,
  encode each modality on it, and test the evidence.
* `chorale_concepts()` builds the vocabulary of named biological concepts a
  collection can be scored on and reports, per concept and modality, how many
  features carry it.
* `chorale_encode()` scores every modality on the vocabulary and keeps what the
  vocabulary does not explain as free dimensions.
* `chorale_concept_evidence()` tests whether each concept separates cases from
  controls in every modality that expresses it and combines modalities by
  inverse-variance weighting, with family-wise and false-discovery control.
* `chorale_concept_families()` and `chorale_family_evidence()` test concepts in
  overlapping groups as well as one at a time.
* `chorale_concept_specificity()` tests each shared covariate in the
  phenotype's place.

## Joint state

* `chorale_joint_state()` estimates one latent state from all modalities at
  once by stacking their concept scores and decomposing them together.
* `chorale_joint_evidence()` tests whether a component moves with the
  phenotype, adjusting for shared covariates and modality.
* `chorale_joint_concepts()` names the concepts whose loadings define each
  component.
* `chorale_joint_transfer()` holds each modality out, fits on the rest, and
  tests whether the projected held-out samples separate cases from controls.

## Free dimensions and factor selection

* `chorale_free_dimensions()` reports coordinated variation no concept explains,
  tested against permuted phenotypes.
* `chorale_n_factors()` and `chorale_select_factors()` choose how many free
  dimensions a modality supports.
* `chorale_ica()` runs the factorisation in the basis the samples support.

## Design, benchmarking and reporting

* `chorale_check_design()` reports what a collection can anchor on before it is
  fitted.
* `chorale_destroy_pairing()` measures recovery of a withheld pairing.
* `chorale_plant()`, `chorale_score_concepts()` and
  `chorale_validate_concepts()` plant named concepts and measure their
  recovery; `chorale_concept_example()` builds example collections.
* `chorale_null()` calibrates the smallest attainable p-value for each control.
* `chorale_report()` writes the concept, free-dimension and added-value tables,
  the joint-state section, and the fitted object as `fit.rds`.
* `chorale_control()` holds the run settings.
