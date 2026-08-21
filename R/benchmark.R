#' Recover a correspondence that is known, with the pairing withheld
#'
#' Every other check asks whether the estimator's output beats a null. This asks
#' the harder question: on data where the true cross-modality correspondence is
#' known, does the estimator find it? Two modalities measured on the *same*
#' individuals carry that answer. The pairing is discarded, the collection is
#' encoded as though the samples were disjoint, and the recovered correspondence
#' is then compared with the pairing that was withheld.
#'
#' The vocabulary asserts that a concept in one modality is the same concept in
#' the other, and it asserts it without seeing a single paired sample. The
#' withheld pairing is what can check that: a concept's score in one modality is
#' compared with its score in the other across the individuals both were
#' measured on.
#' The comparison is signed, because a concept score has a fixed orientation and
#' two modalities ranking the same state in opposite directions have not
#' recovered it.
#'
#' Three quantities frame the result. The **paired upper bound** is the agreement
#' the best alignment the pairing admits reaches, which is what any method
#' knowing the pairing could do. The **random baseline** is the agreement of
#' one-to-one relabellings drawn at random. The recovery is reported between
#' them, and a result near the baseline is a failure however small its p-value.
#'
#' Random alignment is a baseline and not a bound. A vocabulary whose names do
#' not track the same biology in the two modalities places below it, and the
#' placement reports that as a negative value rather than clipping it.
#'
#' @param paired_a,paired_b Feature-by-sample matrices for two modalities
#'   measured on the same individuals, with the same sample identifiers in the
#'   same order.
#' @param design A design table for those individuals, carrying `sample_id` and
#'   the covariates the effects are adjusted for.
#' @param sets A named list of concepts, as returned by [chorale_genesets()],
#'   giving the vocabulary recovery is scored in.
#' @param n_random Random alignments forming the baseline.
#' @param feature_space Passed to [chorale_concepts()], for a modality whose
#'   features reach the vocabulary through lipid class rather than through gene
#'   identifiers.
#' @param min_features Features of a modality a concept must reach to be scored
#'   in it.
#' @param n_free Free dimensions per modality when the collection is encoded.
#' @param n_init Initialisations per encoding.
#' @param seed Integer seed.
#'
#' @returns A list with `concept_summary`, a one-row data frame carrying the
#'   recovered agreement against the paired upper bound and the random baseline;
#'   `concept_recovery`, one row per concept; and `concept_encoding`, the
#'   encoding the comparison was computed from.
#' @export
#' @examples
#' fx <- chorale_concept_example(n_samples = 40, seed = 1)
#' # Give the two modalities the same individuals, so a pairing exists to withhold.
#' a <- SummarizedExperiment::assay(fx$containers$A)
#' b <- SummarizedExperiment::assay(fx$containers$B)
#' colnames(b) <- colnames(a)
#' design <- as.data.frame(
#'   SummarizedExperiment::colData(fx$containers$A))
#' chorale_destroy_pairing(a, b, design, sets = fx$sets,
#'                         n_random = 20)$concept_summary
chorale_destroy_pairing <- function(paired_a, paired_b, design, sets,
                                    n_random = 200L, feature_space = NULL,
                                    min_features = 5L, n_free = 0L,
                                    n_init = 5L, seed = 1L) {
  shared <- intersect(colnames(paired_a), colnames(paired_b))
  shared <- intersect(shared, as.character(design$sample_id))
  if (length(shared) < 20) {
    rlang::abort(paste0(
      "The two modalities share ", length(shared), " samples. A pairing that ",
      "can be withheld needs at least twenty."
    ))
  }
  a <- paired_a[, shared, drop = FALSE]
  b <- paired_b[, shared, drop = FALSE]
  d <- design[match(shared, design$sample_id), , drop = FALSE]

  # The encoder never sees that the samples correspond: each modality is given
  # its own identifiers, as a disjoint pair of cohorts would be.
  da <- d
  db <- d
  da$sample_id <- paste0("a_", seq_along(shared))
  db$sample_id <- paste0("b_", seq_along(shared))
  colnames(a) <- da$sample_id
  colnames(b) <- db$sample_id
  containers <- list(A = chorale_load(a, da), B = chorale_load(b, db))

  concept <- chorale_concept_pairing(containers, sets, feature_space,
                                     min_features, n_free, n_init, n_random,
                                     seed)
  list(concept_summary = concept$summary,
       concept_recovery = concept$recovery,
       concept_encoding = concept$encoding)
}

#' Score the correspondence the vocabulary asserts against the withheld pairing
#'
#' The estimator is given the two modalities as though they were disjoint, so
#' the concept scores are computed without any knowledge that a sample in one is
#' a sample in the other. The pairing is then used, and only then, to ask
#' whether a concept scored in one modality tracks the same concept scored in
#' the other on the same individual.
#'
#' The correspondence the vocabulary asserts is the identity: concept for
#' concept, by name. It is placed between the best alignment the pairing admits,
#' which is what any method knowing the pairing could have found, and alignment
#' at random.
#'
#' @keywords internal
#' @noRd
chorale_concept_pairing <- function(containers, sets, feature_space,
                                    min_features, n_free, n_init, n_random,
                                    seed) {
  concepts <- chorale_concepts(containers, sets, feature_space = feature_space,
                               min_features = min_features,
                               min_modalities = 2L)
  shared <- Reduce(intersect, lapply(concepts$membership, colnames))
  if (length(shared) < 2) {
    return(list(
      summary = data.frame(
        n_concepts = length(shared),
        paired_upper_bound = NA_real_, recovered_agreement = NA_real_,
        random_baseline = NA_real_, placement_between_bounds = NA_real_,
        fraction_partner_correct = NA_real_,
        reason = "fewer than two concepts reach both modalities",
        stringsAsFactors = FALSE),
      recovery = data.frame(), encoding = NULL))
  }

  encoding <- chorale_encode(containers, concepts, n_free = n_free,
                             n_init = n_init, seed = seed)
  sa <- encoding$encodings$A$concept_scores[, shared, drop = FALSE]
  sb <- encoding$encodings$B$concept_scores[, shared, drop = FALSE]

  # The pairing enters here and nowhere else: the two matrices are in the same
  # sample order because the withheld pairing put them there.
  #
  # The correlation is signed. A concept score has a fixed orientation, since a
  # larger score means a larger weighted mean of the features that carry the
  # concept, so two modalities ranking the same biological state in opposite
  # directions have not recovered it. Sign indeterminacy belongs to the free
  # dimensions, which are not what is being compared here.
  agreement <- suppressWarnings(stats::cor(sa, sb))
  agreement[!is.finite(agreement)] <- 0

  # The assignment solver takes non-negative costs. Adding a constant to every
  # entry adds the same amount to every complete assignment, so the optimal
  # assignment is the one the signed agreement implies.
  best <- chorale_assign(agreement + 1)
  paired_bound <- mean(agreement[cbind(seq_along(shared), best)])
  self <- diag(agreement)
  recovered_agreement <- mean(self)
  correct <- mean(best == seq_along(shared))

  # A random alignment is a one-to-one relabelling of the concepts, which is the
  # thing the identity is being compared against. Drawing with replacement would
  # compare against something no alignment could be.
  set.seed(seed)
  random <- vapply(seq_len(n_random), function(i) {
    mean(agreement[cbind(seq_along(shared), sample(seq_along(shared)))])
  }, numeric(1))
  random_baseline <- mean(random)

  # A baseline, not a bound: a vocabulary whose names do not track the same
  # biology in the two modalities can place below it, and the placement says so
  # by going negative.
  placement <- if (paired_bound > random_baseline) {
    (recovered_agreement - random_baseline) / (paired_bound - random_baseline)
  } else {
    NA_real_
  }

  recovery <- data.frame(
    concept = shared,
    self_agreement = round(self, 4),
    best_partner = shared[best],
    best_agreement = round(agreement[cbind(seq_along(shared), best)], 4),
    partner_correct = best == seq_along(shared),
    stringsAsFactors = FALSE
  )
  recovery <- recovery[order(-recovery$self_agreement), , drop = FALSE]
  rownames(recovery) <- NULL

  list(
    summary = data.frame(
      n_concepts = length(shared),
      paired_upper_bound = round(paired_bound, 4),
      recovered_agreement = round(recovered_agreement, 4),
      random_baseline = round(random_baseline, 4),
      placement_between_bounds = round(placement, 3),
      fraction_partner_correct = round(correct, 3),
      reason = NA_character_,
      stringsAsFactors = FALSE),
    recovery = recovery,
    encoding = encoding
  )
}
