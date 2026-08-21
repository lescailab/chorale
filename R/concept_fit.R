#' Estimate the shared state of a collection through its concepts
#'
#' The whole estimator in one call. A vocabulary of named concepts is built for
#' the collection, every modality is scored on it, what the vocabulary does not
#' explain is kept as free dimensions, and the design is asked whether each
#' concept separates cases from controls.
#'
#' The order matters and is fixed. Scoring comes before the design is consulted,
#' so the concept scores cannot have been shaped by the contrast they are then
#' tested on. The vocabulary is fixed before any of it, so which concepts are
#' tested does not depend on which of them turn out to be interesting.
#'
#' What connects the modalities is the vocabulary, not the design. Two
#' modalities measured on different animals share no sample and need share no
#' covariate beyond the phenotype: they meet because both are scored on the same
#' named concepts. The design is what supplies the evidence, one modality at a
#' time, that a concept moves with the disease.
#'
#' @param containers A named list of [chorale_load()] containers.
#' @param sets A named list of concepts, as returned by [chorale_genesets()], or
#'   a `chorale_concepts` object already built for this collection.
#' @param feature_space,feature_map Passed to [chorale_concepts()].
#' @param n_free,transform Passed to [chorale_encode()].
#' @param n_permutations Permutations calibrating the vocabulary.
#' @param min_features Features of a modality a concept must reach before that
#'   modality is counted as expressing it.
#' @param min_modalities Modalities a concept must reach to enter the
#'   vocabulary.
#' @param assay_name Assay to read from each container.
#' @param control A [chorale_control()] object.
#' @param seed Integer seed.
#' @param ... Named overrides applied to `control`.
#'
#' @returns An object of class `chorale_concept_fit` carrying `concepts`,
#'   `encoding` and `evidence`.
#' @export
#' @examples
#' fx <- chorale_concept_example(seed = 1)
#' chorale_concept_fit(fx$containers, fx$sets, n_free = 1,
#'                     n_permutations = 99, n_init = 2)
chorale_concept_fit <- function(containers, sets,
                                feature_space = NULL,
                                feature_map = NULL,
                                n_free = "auto",
                                transform = "auto",
                                n_permutations = 999L,
                                min_features = 5L,
                                min_modalities = 1L,
                                assay_name = NULL,
                                control = chorale_control(),
                                seed = 1L, ...) {
  if (!is.list(containers) || length(containers) < 1) {
    rlang::abort("`containers` must be a list of at least one modality.")
  }
  if (is.null(names(containers))) {
    names(containers) <- paste0("modality_", seq_along(containers))
  }
  control <- chorale_merge_control(control, list(...))

  concepts <- if (inherits(sets, "chorale_concepts")) {
    sets
  } else {
    chorale_concepts(containers, sets, feature_space = feature_space,
                     feature_map = feature_map, min_features = min_features,
                     min_modalities = min_modalities,
                     min_lipid_compounds = control$min_lipid_compounds,
                     min_lipid_specificity = control$min_lipid_specificity,
                     assay_name = assay_name)
  }
  if (length(concepts$vocabulary) == 0) {
    rlang::abort(paste0(
      "No concept reaches ", min_modalities, " modality with at least ",
      min_features, " features. The collection cannot be scored on this ",
      "vocabulary."), class = "chorale_empty_vocabulary")
  }

  encoding <- chorale_encode(containers, concepts, n_free = n_free,
                             transform = transform, assay_name = assay_name,
                             control = control, seed = seed)
  evidence <- chorale_concept_evidence(encoding, n_permutations = n_permutations,
                                       control = control, seed = seed)

  structure(
    list(
      modalities = names(containers),
      concepts = concepts,
      encoding = encoding,
      evidence = evidence,
      designs = encoding$designs,
      feature_space = concepts$feature_space,
      control = control,
      n_permutations = as.integer(n_permutations),
      seed = seed
    ),
    class = "chorale_concept_fit"
  )
}

#' A small collection with a concept planted in every modality
#'
#' Documentation and tests both need a collection where a named concept really
#' does separate cases from controls in more than one modality, on features the
#' modalities share by name. This builds one directly rather than deriving it
#' from a factor model, so what the concept does is stated rather than inferred.
#'
#' @param n_samples Samples per modality.
#' @param n_features Features per modality.
#' @param effect Shift, in standard deviations, applied to the members of the
#'   planted concept in cases.
#' @param seed Integer seed.
#'
#' @returns A list with `containers`, a named list of two modalities on disjoint
#'   samples; `sets`, three concepts of which the first is the planted one; and
#'   `planted`, its name.
#' @export
#' @examples
#' fx <- chorale_concept_example()
#' names(fx$sets)
chorale_concept_example <- function(n_samples = 60L, n_features = 90L,
                                    effect = 1, seed = 1L) {
  if (n_features < 45L) {
    rlang::abort("`n_features` must be at least 45 for three concepts and a tail.")
  }
  set.seed(seed)
  ids <- sprintf("feature_%04d", seq_len(n_features))
  # Three concepts of equal size separated by gaps, and a tail of features no
  # concept covers, so the vocabulary is genuinely incomplete.
  k <- floor((n_features - 15L) / 3L)
  sets <- list(planted = ids[seq_len(k)],
               quiet = ids[k + 5L + seq_len(k)],
               other = ids[2L * k + 10L + seq_len(k)])
  tail_ids <- ids[(3L * k + 11L):n_features]
  containers <- list()
  for (i in seq_len(2L)) {
    # Disjoint samples: the two modalities share no animal, which is the
    # situation the estimator exists for.
    sample_id <- sprintf("m%d_s%03d", i, seq_len(n_samples))
    phenotype <- rep(c("control", "case"), length.out = n_samples)
    sex <- rep(c("F", "F", "M", "M"), length.out = n_samples)
    x <- matrix(stats::rnorm(n_features * n_samples), nrow = n_features,
                dimnames = list(ids, sample_id))
    x[sets$planted, phenotype == "case"] <-
      x[sets$planted, phenotype == "case"] + effect
    # A modality-private direction, so the free dimensions have something to
    # carry that no concept explains.
    private <- stats::rnorm(n_samples)
    x[tail_ids, ] <- x[tail_ids, ] +
      matrix(rep(private, each = length(tail_ids)), nrow = length(tail_ids))
    containers[[i]] <- chorale_load(
      x, data.frame(sample_id = sample_id, phenotype = phenotype, sex = sex,
                    stringsAsFactors = FALSE))
  }
  names(containers) <- c("A", "B")
  list(containers = containers, sets = sets, planted = "planted")
}

#' @export
print.chorale_concept_fit <- function(x, ...) {
  cat("<chorale_concept_fit>\n")
  cat("  modalities:", paste(x$modalities, collapse = ", "), "\n")
  cat("  concepts in the vocabulary:", length(x$concepts$vocabulary),
      sprintf("(%d expressible in every modality)\n", x$concepts$n_in_all))
  v <- x$encoding$variance
  cat("  free dimensions:",
      paste(sprintf("%s %d", v$modality, v$n_free), collapse = ", "), "\n")
  j <- x$evidence$joint
  if (nrow(j) > 0) {
    cat("  concepts supported at q <=", x$evidence$alpha, ":",
        sum(j$significant), "\n")
    cat("  strongest:", j$concept[1], sprintf("(z = %.2f)\n", j$joint_z[1]))
  } else {
    cat("  no concept could be tested on this design\n")
  }
  invisible(x)
}
