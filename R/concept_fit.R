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

#' A collection with named concepts planted in it
#'
#' Documentation, tests and recovery curves all need a collection where named
#' concepts really do separate cases from controls, on features the modalities
#' share by name. This builds one directly rather than deriving it from a factor
#' model, so what each concept does is stated rather than inferred.
#'
#' Coverage is a dial rather than a fixture. A concept is defined on the whole
#' vocabulary, but only a share of its members is measured in any one modality:
#' the rest carry identifiers of that modality's own, as an assay that does not
#' quantify every gene would. Lowering `coverage` is therefore the question of
#' how far the vocabulary has to reach into a modality before a concept planted
#' in it can still be recovered.
#'
#' @param n_samples Samples per modality.
#' @param n_features Features per modality.
#' @param effect Shift, in standard deviations, applied in cases to the members
#'   of every planted concept.
#' @param n_modalities Modalities in the collection, on disjoint samples.
#' @param n_concepts Concepts in the vocabulary.
#' @param n_planted Concepts of those that carry the case-control shift. The
#'   rest are the vocabulary the estimator must leave alone.
#' @param coverage Share of each concept's members a modality measures.
#' @param seed Integer seed.
#'
#' @returns A list with `containers`, a named list of modalities on disjoint
#'   samples; `sets`, the vocabulary; and `planted`, the names of the concepts
#'   that carry the shift.
#' @export
#' @examples
#' fx <- chorale_concept_example()
#' names(fx$sets)
#' fx$planted
chorale_concept_example <- function(n_samples = 60L, n_features = 90L,
                                    effect = 1, n_modalities = 2L,
                                    n_concepts = 3L, n_planted = 1L,
                                    coverage = 1, seed = 1L) {
  if (n_planted > n_concepts) {
    rlang::abort("`n_planted` cannot exceed `n_concepts`.")
  }
  if (coverage <= 0 || coverage > 1) {
    rlang::abort("`coverage` must lie in (0, 1].")
  }
  set.seed(seed)
  ids <- sprintf("feature_%04d", seq_len(n_features))
  # Concepts of equal size separated by gaps, and a tail of features no concept
  # covers, so the vocabulary is genuinely incomplete.
  k <- floor((n_features - 5L * (n_concepts + 1L)) / n_concepts)
  if (k < 8L) {
    rlang::abort(paste0(
      "`n_features` is too small for ", n_concepts,
      " concepts; each would reach fewer than eight features."))
  }
  names_of <- chorale_example_names(n_concepts, n_planted)
  sets <- stats::setNames(lapply(seq_len(n_concepts), function(j) {
    ids[(j - 1L) * (k + 5L) + 5L + seq_len(k)]
  }), names_of)
  covered <- unlist(sets, use.names = FALSE)
  tail_ids <- setdiff(ids, covered)
  planted <- names_of[seq_len(n_planted)]

  containers <- list()
  for (i in seq_len(n_modalities)) {
    # Disjoint samples: the modalities share no animal, which is the situation
    # the estimator exists for.
    sample_id <- sprintf("m%d_s%03d", i, seq_len(n_samples))
    phenotype <- rep(c("control", "case"), length.out = n_samples)
    sex <- rep(c("F", "F", "M", "M"), length.out = n_samples)
    x <- matrix(stats::rnorm(n_features * n_samples), nrow = n_features,
                dimnames = list(ids, sample_id))
    for (cn in planted) {
      x[sets[[cn]], phenotype == "case"] <-
        x[sets[[cn]], phenotype == "case"] + effect
    }
    # A modality-private direction, so the free dimensions have something to
    # carry that no concept explains.
    if (length(tail_ids) > 0) {
      private <- stats::rnorm(n_samples)
      x[tail_ids, ] <- x[tail_ids, ] +
        matrix(rep(private, each = length(tail_ids)), nrow = length(tail_ids))
    }
    if (coverage < 1) {
      # A member the modality does not measure is a feature under an identifier
      # of that modality's own, so the concept reaches fewer features here.
      unmeasured <- unlist(lapply(sets, function(members) {
        n_drop <- length(members) - max(1L, floor(length(members) * coverage))
        if (n_drop <= 0) character() else sample(members, n_drop)
      }), use.names = FALSE)
      at <- match(unmeasured, rownames(x))
      rownames(x)[at] <- paste0("m", i, "_", unmeasured)
    }
    containers[[i]] <- chorale_load(
      x, data.frame(sample_id = sample_id, phenotype = phenotype, sex = sex,
                    stringsAsFactors = FALSE))
  }
  names(containers) <- LETTERS[seq_len(n_modalities)]
  list(containers = containers, sets = sets, planted = planted)
}

#' Names for an example vocabulary, planted concepts first
#' @keywords internal
#' @noRd
chorale_example_names <- function(n_concepts, n_planted) {
  planted <- if (n_planted >= 1L) {
    c("planted", if (n_planted > 1L) paste0("planted_", 2:n_planted))
  } else {
    character()
  }
  n_quiet <- n_concepts - n_planted
  quiet <- c("quiet", "other")[seq_len(min(2L, max(n_quiet, 0L)))]
  if (n_quiet > 2L) quiet <- c(quiet, paste0("concept_", 3:n_quiet))
  c(planted, quiet)
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
