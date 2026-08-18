---
layout: ../layouts/Base.astro
title: "Why CHORALE exists"
description: "What integration across disjoint cohorts requires, what existing methods provide, and what is missing."
---

# Why CHORALE exists

<p class="lede">What integration across disjoint cohorts requires, what existing methods provide, and what is missing.</p>

## Why disjoint samples are the common case

Multi-omics integration is usually described for a design in which every individual is measured on
every platform. That design is expensive, and it is rare in animal models, where tissue is consumed
by the assay: a brain used for transcriptomics is not available for proteomics, and an animal
sacrificed at one age cannot be sampled at another. The consequence is that most published omics
resources for a given disease model describe the same biology in different animals.

Integration across such resources is worth attempting because the modalities are complementary.
Transcript abundance reports regulation, protein abundance reports what the regulation produced
after turnover, and metabolite or lipid abundance reports the biochemical consequence. A disease
process can be visible in the combination while no single layer carries enough of it to be called.

## What existing methods do

**Matched-sample factor models.** MOFA and its successors recover latent factors shared across
assays by exploiting the fact that the same individual appears in each (Argelaguet et al. 2018).
The individual is the key that joins the tables. Without it these methods have nothing to join on.

**Methods that learn from bridge data.** INTEND trains a cross-omic predictor on paired samples and
then applies it to unpaired ones (Itai et al. 2023). It is effective and it was benchmarked
honestly, by withholding known pairing across eleven cancer datasets. It requires paired data to
exist for the tissue and platform in question, which for most animal-model resources it does not.

**Methods that join on shared features.** mosaicMPI factorises heterogeneous datasets separately
and connects the resulting programmes through correlations of their loadings on features the
datasets share (Verhey et al. 2024). This requires a common feature space, which transcripts and
lipids do not have.

**Methods that anchor on a shared label.** Propensity score alignment uses shared treatment or
perturbation information to align unpaired modalities (Xi et al. 2024); a dual-channel transformer integrates unpaired
transcriptomic and epigenomic data by cross-modal reconstruction (Liao et al. 2026); UnCOT-AD couples representations across
three unpaired omics layers in Alzheimer's disease (Abir et al. 2025). Anchoring on an experimental label is
therefore an established idea rather than a new one.

**Identification theory for unpaired data.** Sturma et al. (2023) show that a shared latent block
is identifiable from unpaired multi-domain data when the modality-specific components are
non-Gaussian and pairwise different, and give a finite-sample procedure that runs independent
component analysis separately and matches the recovered component distributions with
Kolmogorov-Smirnov tests. Timilsina et al. (2024) identify shared components under a modality
variability assumption using distributional discrepancies inside a joint estimator.

## What is missing

Four gaps run through that body of work.

**The unidentified quantity is rarely named.** With disjoint samples, the joint distribution of two
modalities is not identified: the data fix each marginal and say nothing about how they are
coupled within an individual. Methods that construct pseudo-pairs produce a number for that
coupling. The number is an artefact of the construction, and nothing in the data distinguishes it
from the alternatives.

**Correspondence is decided pair by pair.** Where more than two modalities are involved, the usual
route is to compare each pair and reconcile the results. Independent pairwise decisions need not
agree with one another, and reconciling them by merging propagates a single wrong pairing through
the whole result.

**Biology enters as annotation rather than as evidence.** Programmes are typically recovered on
statistical grounds and named afterwards from curated sets. The naming then cannot corroborate the
recovery, because it played no part in it and because any factor yields some enrichment.

**Metabolomes and lipidomes sit outside the shared vocabulary.** Genes and proteins reach curated
pathway sets through their identifiers. Lipid species do not, so a lipidome is usually analysed in
a vocabulary of its own and compared with the others by hand.

## What CHORALE does about them

**It reports the identified set instead of inventing a coupling.** Where the data cannot determine
how strongly two programmes move together within an individual, the output is the range of
correlations the data cannot exclude, narrowed as far as the experimental design allows. Where that
range stays wide, its width is the result.

**It solves the correspondence once, over all modalities together.** Every pairwise agreement enters
one matrix and the assignment of factors to programmes is recovered from it in a single step, so
the correspondences agree around every cycle and no programme is assembled by chaining.

**It corroborates on biology through a separate channel.** Factors are fitted without ever seeing
the curated sets. Whether the modalities implicate the same biology is then asked as a second,
independent question, against a null that holds set sizes and feature frequencies fixed. The two
channels are never combined into one score, because keeping them apart is what shows that the
agreement on biology was not built in.

**It brings a lipidome into the same vocabulary.** A table shipped with the package maps lipid
classes onto the pathway sets the genes are described in, so a lipid species and a gene are placed
in one set. The transcriptome then contributes the enzymes of a pathway and the lipidome the
molecules that pathway acts on.

**It reports failure as an outcome.** A run in which no programme survives its null says so.

## Where the value lies

The combination is what is distinctive: an assignment solved jointly, two evidence channels that
can disagree, identified sets in place of a fabricated coupling, calibration that reruns the whole
selection procedure, and outputs a reader can inspect rather than a score from a black box. Each
element has precedent; their assembly into an auditable workflow for disjoint bulk-omics cohorts
does not, to the extent a search of the applied literature reveals.
