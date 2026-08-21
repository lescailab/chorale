#' The named concepts a collection can be scored on
#'
#' Modalities measured on different individuals share no samples, so something
#' other than a sample must connect them. That something is a fixed vocabulary
#' of named biological concepts: a concept is a set of features, defined once
#' for the collection, and every modality that carries enough of those features
#' can be scored on it. A concept is the same concept in every modality, so no
#' correspondence between modalities has to be estimated.
#'
#' Features reach a concept by the route their identifiers allow. Genes and
#' proteins reach it through gene identifiers, by way of
#' [chorale_geneset_matrix()]. Lipid species reach it through their class, by
#' way of [chorale_metabolite_matrix()], since a pathway acts on a class rather
#' than on one chain length. Either way the concept keeps its name, so the two
#' modalities end up in one vocabulary rather than in two that have to be
#' reconciled afterwards.
#'
#' Coverage is reported rather than assumed. A concept expressible in every
#' modality of the collection can carry evidence across all of them; one
#' expressible in a subset carries evidence across that subset; one no modality
#' expresses carries none, and is reported as absent rather than silently
#' dropped, because an absence is a statement about the measured features.
#'
#' @param containers A named list of [chorale_load()] containers, one per
#'   modality.
#' @param sets A named list of concepts, each a character vector of member
#'   identifiers, as returned by [chorale_genesets()].
#' @param feature_space Optional named character vector saying how each
#'   modality's features reach the vocabulary: `"gene"` through gene
#'   identifiers, `"lipid"` through lipid class. Modalities not named default to
#'   `"gene"`.
#' @param feature_map Optional named list, one entry per modality, each a
#'   mapping from [chorale_map()] carrying `id`, `ENTREZID` and `weight`. Where
#'   a modality supplies one, a feature identifier matching several genes
#'   contributes proportionally rather than counting once per gene.
#' @param min_features Features of a modality a concept must reach before that
#'   modality is counted as expressing it. A concept matching one or two
#'   measured features is matched by chance rather than measured.
#' @param min_modalities Modalities a concept must reach to enter the
#'   vocabulary. The default keeps a concept reaching one modality, since it
#'   still carries evidence there; raising it to two keeps only concepts that
#'   can connect modalities.
#' @param min_lipid_compounds,min_lipid_specificity Passed to
#'   [chorale_metabolite_matrix()] for modalities whose feature space is
#'   `"lipid"`.
#' @param assay_name Assay to read from each container. Defaults to the first.
#'
#' @returns An object of class `chorale_concepts` with `vocabulary`, the
#'   concepts retained; `membership`, one feature-by-concept matrix per
#'   modality; `coverage`, one row per concept per modality carrying the number
#'   of features that modality has for it; and `summary`, one row per concept
#'   carrying how many modalities express it and whether all do.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' ids <- rownames(sim$modalities[[1]])
#' sets <- list(concept_a = ids[1:20], concept_b = ids[15:40])
#' chorale_concepts(containers, sets, min_features = 5)
chorale_concepts <- function(containers, sets,
                             feature_space = NULL,
                             feature_map = NULL,
                             min_features = 5L,
                             min_modalities = 1L,
                             min_lipid_compounds = 2L,
                             min_lipid_specificity = 0.05,
                             assay_name = NULL) {
  if (!is.list(containers) || length(containers) < 1) {
    rlang::abort("`containers` must be a list of at least one modality.")
  }
  if (is.null(names(containers))) {
    names(containers) <- paste0("modality_", seq_along(containers))
  }
  if (!is.list(sets) || is.null(names(sets)) || length(sets) == 0) {
    rlang::abort("`sets` must be a named list of concepts.")
  }
  modalities <- names(containers)
  space <- chorale_feature_space(feature_space, modalities)

  membership <- list()
  for (m in modalities) {
    se <- containers[[m]]
    an <- assay_name %||% SummarizedExperiment::assayNames(se)[1]
    ids <- rownames(SummarizedExperiment::assay(se, an))
    if (is.null(ids)) {
      rlang::abort(paste0(
        "Modality '", m, "' has no feature identifiers, so nothing can be ",
        "placed in the vocabulary."))
    }
    mm <- if (identical(unname(space[[m]]), "lipid")) {
      chorale_metabolite_matrix(
        ids, sets, min_compounds = min_lipid_compounds,
        min_specificity = min_lipid_specificity, min_features = 0L)
    } else {
      mapping <- if (!is.null(feature_map)) feature_map[[m]] else NULL
      out <- chorale_geneset_matrix(ids, sets, mapping = mapping,
                                    min_features = 0L)
      if (ncol(out) > 0) rownames(out) <- ids
      out
    }
    rownames(mm) <- ids
    membership[[m]] <- mm
  }

  # Membership is built without a threshold, so coverage records how many
  # features a modality actually has for a concept and the threshold is applied
  # afterwards. Building it thresholded would report zero for every concept the
  # threshold removed, which is a different statement. The grid is complete, so
  # a concept no modality expresses is a row carrying zero rather than a name
  # that quietly disappeared.
  coverage <- do.call(rbind, lapply(names(sets), function(cn) {
    data.frame(
      concept = cn,
      modality = modalities,
      n_features = vapply(modalities, function(m) {
        mm <- membership[[m]]
        if (!cn %in% colnames(mm)) return(0L)
        as.integer(sum(mm[, cn] > 0))
      }, integer(1)),
      stringsAsFactors = FALSE
    )
  }))
  coverage$expressed <- coverage$n_features >= min_features
  rownames(coverage) <- NULL

  reached <- vapply(split(coverage$expressed, coverage$concept), sum, integer(1))
  summary <- data.frame(
    concept = names(reached),
    n_modalities = as.integer(reached),
    in_all_modalities = as.integer(reached) == length(modalities),
    stringsAsFactors = FALSE
  )
  summary <- summary[order(-summary$n_modalities, summary$concept), ,
                     drop = FALSE]
  rownames(summary) <- NULL

  vocabulary <- summary$concept[summary$n_modalities >= min_modalities]
  # Each modality keeps only the concepts it expresses, so nothing downstream
  # scores a modality on a concept it has no features for.
  membership <- lapply(modalities, function(m) {
    mm <- membership[[m]]
    keep <- intersect(
      colnames(mm),
      coverage$concept[coverage$modality == m & coverage$expressed])
    mm[, intersect(vocabulary, keep), drop = FALSE]
  })
  names(membership) <- modalities

  structure(
    list(
      modalities = modalities,
      feature_space = space,
      vocabulary = vocabulary,
      membership = membership,
      coverage = coverage,
      summary = summary,
      n_supplied = length(sets),
      n_in_all = sum(summary$in_all_modalities),
      n_in_some = sum(summary$n_modalities > 0 & !summary$in_all_modalities),
      n_in_none = sum(summary$n_modalities == 0),
      min_features = as.integer(min_features),
      min_modalities = as.integer(min_modalities)
    ),
    class = "chorale_concepts"
  )
}

#' @export
print.chorale_concepts <- function(x, ...) {
  cat("<chorale_concepts>\n")
  cat("  modalities:", paste(x$modalities, collapse = ", "), "\n")
  cat("  concepts supplied:", x$n_supplied, "\n")
  cat("  expressible in every modality:", x$n_in_all, "\n")
  cat("  expressible in some but not all:", x$n_in_some, "\n")
  cat("  expressible in none:", x$n_in_none, "\n")
  cat("  in the vocabulary (at least", x$min_modalities, "modalities):",
      length(x$vocabulary), "\n")
  cat("  features a modality needs for a concept:", x$min_features, "\n")
  for (m in x$modalities) {
    cat(sprintf("    %-12s %d concepts\n", m, ncol(x$membership[[m]])))
  }
  invisible(x)
}
