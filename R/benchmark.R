#' Recover a correspondence that is known, with the pairing withheld
#'
#' Every other check asks whether the estimator's output beats a null. This asks
#' the harder question: on data where the true cross-modality correspondence is
#' known, does the estimator find it? Two modalities measured on the *same*
#' individuals carry that answer. The pairing is discarded, the estimator is run
#' as though the samples were disjoint, and the recovered correspondence is then
#' compared with the pairing that was withheld.
#'
#' The truth is a correspondence between factors, not between samples: a factor
#' recovered in one modality and a factor recovered in the other are the same
#' programme when their scores agree across the animals both were measured on.
#' That agreement is computable only because the pairing exists, and it is
#' withheld from the estimator throughout.
#'
#' In concept space the question is sharper. The vocabulary asserts that a
#' concept in one modality is the same concept in the other, and that assertion
#' is made without seeing a single paired sample. The withheld pairing is what
#' can check it: a concept's score in one modality is compared with its score in
#' the other across the animals both were measured on, and the comparison the
#' vocabulary asserts is placed against the best alignment any method could have
#' found and against alignment at random. The comparison is signed, because a
#' concept score has a fixed orientation and two modalities ranking the same
#' state in opposite directions have not recovered it.
#'
#' Three quantities frame the result, in the manner the plan sets out. The
#' **paired benchmark** is the correspondence obtained by using the pairing,
#' which is the best any method could do. **Random alignment** is the
#' correspondence obtained by pairing factors at random, which is the worst.
#' The estimator's recovery is reported between them, and a result near random
#' alignment is a failure however small its p-value.
#'
#' In concept space random alignment is a baseline rather than a bound. The
#' identity the vocabulary asserts can do worse than relabelling the concepts at
#' random, which is what a vocabulary whose names do not track the same biology
#' in both modalities looks like, and the placement reports it as a negative
#' value rather than clipping it.
#'
#' @param paired_a,paired_b Feature-by-sample matrices for two modalities
#'   measured on the same individuals, with the same sample identifiers in the
#'   same order.
#' @param design A design table for those individuals, carrying `sample_id` and
#'   the covariates the estimator anchors on.
#' @param n_factors Factors per modality, or `"auto"`.
#' @param n_init Initialisations per fit.
#' @param n_random Random alignments forming the lower bound.
#' @param sets A named list of concepts, as returned by [chorale_genesets()].
#'   Supplying it adds the concept-space benchmark; without it only the factor
#'   correspondence is scored.
#' @param spaces Which spaces to score, `"factor"`, `"concept"`, or both.
#'   Defaults to whichever `sets` makes available.
#' @param feature_space Passed to [chorale_concepts()], for a modality whose
#'   features reach the vocabulary through lipid class rather than through gene
#'   identifiers.
#' @param min_features Features of a modality a concept must reach to be scored
#'   in it.
#' @param n_free Free dimensions per modality when the collection is encoded.
#' @param seed Integer seed.
#'
#' @returns A list with `truth`, the factor correspondence the pairing implies;
#'   `recovered`, the correspondence the estimator found without it; and
#'   `summary`, a one-row data frame carrying the recovered agreement against
#'   the paired upper bound and the random lower bound, and the fraction of
#'   factors whose partner was recovered correctly. Where the concept space is
#'   scored, `concept_summary` carries the same three quantities for the
#'   vocabulary and `concept_recovery` one row per concept.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' # Give the two modalities the same animals, so a pairing exists to withhold.
#' a <- sim$modalities[[1]]
#' b <- sim$modalities[[2]][, seq_len(ncol(a))]
#' colnames(b) <- colnames(a)
#' chorale_destroy_pairing(a, b, sim$col_data[[1]], n_factors = 3, n_init = 2,
#'                         n_random = 20)$summary
chorale_destroy_pairing <- function(paired_a, paired_b, design,
                                    n_factors = "auto", n_init = 5L,
                                    n_random = 200L, sets = NULL,
                                    spaces = NULL, feature_space = NULL,
                                    min_features = 5L, n_free = 0L,
                                    seed = 1L) {
  spaces <- spaces %||% (if (is.null(sets)) "factor" else c("factor", "concept"))
  spaces <- match.arg(spaces, c("factor", "concept"), several.ok = TRUE)
  if ("concept" %in% spaces && is.null(sets)) {
    rlang::abort("Scoring the concept space needs `sets`, the vocabulary to score in.")
  }
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

  # The estimator never sees that the samples correspond: each modality is given
  # its own identifiers, as a disjoint pair of cohorts would be.
  da <- d
  db <- d
  da$sample_id <- paste0("a_", seq_along(shared))
  db$sample_id <- paste0("b_", seq_along(shared))
  colnames(a) <- da$sample_id
  colnames(b) <- db$sample_id
  containers <- list(A = chorale_load(a, da), B = chorale_load(b, db))

  concept <- if ("concept" %in% spaces) {
    chorale_concept_pairing(containers, sets, feature_space, min_features,
                            n_free, n_init, n_random, seed)
  } else {
    NULL
  }
  if (!"factor" %in% spaces) {
    return(list(truth = data.frame(), recovered = data.frame(),
                summary = data.frame(), fit = NULL,
                concept_summary = concept$summary,
                concept_recovery = concept$recovery,
                concept_encoding = concept$encoding))
  }

  fit <- chorale_fit(containers, n_factors = n_factors, n_init = n_init,
                     n_pathway_perm = 0L, seed = seed)

  # The withheld truth: factors correspond when their scores agree across the
  # animals both were measured on.
  sa <- fit$fits$A$scores
  sb <- fit$fits$B$scores
  agreement <- abs(suppressWarnings(stats::cor(sa, sb)))
  agreement[!is.finite(agreement)] <- 0
  truth_assign <- chorale_assign(agreement)
  truth <- data.frame(
    factor_a = colnames(sa),
    factor_b = colnames(sb)[truth_assign],
    agreement = round(agreement[cbind(seq_len(nrow(agreement)),
                                      truth_assign)], 4),
    stringsAsFactors = FALSE
  )
  paired_bound <- mean(truth$agreement, na.rm = TRUE)

  # What the phenotype-led estimator supported, knowing nothing of the pairing.
  m <- fit$matches
  recovered <- if (is.data.frame(m) && nrow(m) > 0) {
    eligible <- m$supported & m$resolution_status %in% c("resolved", "ambiguous")
    m <- m[eligible, , drop = FALSE]
    data.frame(factor_a = m$factor_a, factor_b = m$factor_b,
               significant = m$significant, stringsAsFactors = FALSE)
  } else {
    data.frame(factor_a = character(0), factor_b = character(0),
               significant = logical(0), stringsAsFactors = FALSE)
  }
  recovered_agreement <- if (nrow(recovered) > 0) {
    mean(vapply(seq_len(nrow(recovered)), function(i) {
      agreement[recovered$factor_a[i], recovered$factor_b[i]]
    }, numeric(1)), na.rm = TRUE)
  } else {
    NA_real_
  }
  correct <- if (nrow(recovered) > 0) {
    key <- stats::setNames(truth$factor_b, truth$factor_a)
    mean(recovered$factor_b == key[recovered$factor_a], na.rm = TRUE)
  } else {
    NA_real_
  }

  # The lower bound: what pairing factors at random achieves on the same data.
  set.seed(seed)
  random <- vapply(seq_len(n_random), function(i) {
    j <- sample(seq_len(ncol(agreement)), nrow(agreement), replace = TRUE)
    mean(agreement[cbind(seq_len(nrow(agreement)), j)])
  }, numeric(1))
  random_bound <- mean(random)

  # Where the recovery sits between the bounds: one is the paired fit, zero is
  # random alignment.
  placement <- if (is.finite(recovered_agreement) &&
                   paired_bound > random_bound) {
    (recovered_agreement - random_bound) / (paired_bound - random_bound)
  } else {
    NA_real_
  }

  summary <- data.frame(
    n_samples = length(shared),
    n_factors_a = ncol(sa),
    n_factors_b = ncol(sb),
    paired_upper_bound = round(paired_bound, 4),
    recovered_agreement = round(recovered_agreement, 4),
    random_lower_bound = round(random_bound, 4),
    placement_between_bounds = round(placement, 3),
    fraction_partner_correct = round(correct, 3),
    n_recovered_pairs = nrow(recovered),
    stringsAsFactors = FALSE
  )
  list(truth = truth, recovered = recovered, summary = summary, fit = fit,
       concept_summary = concept$summary, concept_recovery = concept$recovery,
       concept_encoding = concept$encoding)
}

#' Score the correspondence the vocabulary asserts against the withheld pairing
#'
#' The estimator is given the two modalities as though they were disjoint, so
#' the concept scores are computed without any knowledge that a sample in one is
#' a sample in the other. The pairing is then used, and only then, to ask
#' whether a concept scored in one modality tracks the same concept scored in
#' the other on the same animal.
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
