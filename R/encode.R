#' Score every modality on the shared concepts, and keep what they miss
#'
#' Each modality is projected onto the vocabulary from [chorale_concepts()],
#' giving one score per concept per sample: the average standardised level of
#' the features that carry that concept. The score is a number about a named
#' piece of biology, and it means the same thing in every modality, which is
#' what lets modalities measured on different individuals be compared at all.
#'
#' What the vocabulary does not explain is not discarded. The concept scores are
#' regressed out of the assay and independent component analysis is run on what
#' remains, so coordinated variation that no concept accounts for stays visible
#' as free dimensions rather than being absorbed into a concept that does not
#' describe it. The free dimensions are orthogonal to the concept scores by
#' construction, which is what stops the two channels reporting the same signal
#' twice.
#'
#' Encoding never sees the phenotype. Nothing in this function reads the design,
#' so the scores cannot have been shaped by the contrast they are later asked
#' about.
#'
#' @param containers A named list of [chorale_load()] containers.
#' @param concepts A `chorale_concepts` object for the same collection.
#' @param n_free Free dimensions per modality: one number, one per modality, a
#'   named vector, or `"auto"` to take the count from parallel analysis on the
#'   residual, which is the largest count the data can support.
#' @param transform Per-modality scale handling, as in [chorale_transform()].
#' @param assay_name Assay to read from each container. Defaults to the first.
#' @param control A [chorale_control()] object.
#' @param seed Integer seed.
#' @param ... Named overrides applied to `control`.
#'
#' @returns An object of class `chorale_encode` with `encodings`, one entry per
#'   modality carrying `concept_scores`, `free_scores` and `free_loadings`;
#'   `variance`, one row per modality giving the share of variance the concepts
#'   carry, the share the free dimensions carry and the share left over; and
#'   `concept_variance`, one row per concept per modality.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60, seed = 1)
#' ids <- sprintf("feature_%05d", seq_len(60))
#' sim$modalities <- lapply(sim$modalities, function(m) {
#'   rownames(m) <- ids
#'   m
#' })
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' sets <- list(one = ids[1:20], two = ids[21:40])
#' cc <- chorale_concepts(containers, sets, min_features = 5)
#' chorale_encode(containers, cc, n_free = 2, n_init = 2)
chorale_encode <- function(containers, concepts, n_free = "auto",
                           transform = "auto", assay_name = NULL,
                           control = chorale_control(), seed = 1L, ...) {
  if (!inherits(concepts, "chorale_concepts")) {
    rlang::abort("`concepts` must be a chorale_concepts object.")
  }
  if (is.null(names(containers))) {
    names(containers) <- paste0("modality_", seq_along(containers))
  }
  modalities <- names(containers)
  if (!identical(sort(modalities), sort(concepts$modalities))) {
    rlang::abort("`concepts` was built on a different set of modalities.")
  }
  control <- chorale_merge_control(control, list(...))
  transform_of <- chorale_transform_spec(transform, modalities)
  free_of <- chorale_free_spec(n_free, modalities)

  encodings <- list()
  designs <- list()
  variance <- list()
  concept_variance <- list()

  for (m in modalities) {
    se <- containers[[m]]
    an <- assay_name %||% SummarizedExperiment::assayNames(se)[1]
    mat <- SummarizedExperiment::assay(se, an)
    tf <- chorale_transform(mat, transform = transform_of[[m]])
    x <- scale(t(tf$matrix))
    x[!is.finite(x)] <- 0

    scored <- chorale_concept_scores(x, concepts$membership[[m]])
    scores <- scored$scores
    total <- sum(x^2)

    fitted <- chorale_concept_reconstruction(x, scores)
    residual <- x - fitted
    concept_share <- if (total > 0) 1 - sum(residual^2) / total else NA_real_

    k <- free_of[[m]]
    ceiling_k <- chorale_n_factors(residual, quantile = control$n_factors_quantile,
                                   max_factors = control$max_factors, seed = seed)
    if (identical(k, "auto")) {
      k <- ceiling_k
    } else {
      k <- as.integer(k)
    }

    free <- NULL
    free_share <- 0
    if (is.finite(k) && k >= 1L) {
      free <- chorale_ica(residual, k, n_init = control$n_init, seed = seed,
                          consensus = control$consensus)
      approx <- free$scores %*% t(free$loadings)
      free_share <- if (total > 0) {
        (sum(residual^2) - sum((residual - approx)^2)) / total
      } else {
        NA_real_
      }
    }

    free_scores <- if (is.null(free)) {
      matrix(numeric(0), nrow = nrow(x), ncol = 0,
             dimnames = list(rownames(x), NULL))
    } else {
      free$scores
    }
    free_loadings <- if (is.null(free)) {
      matrix(numeric(0), nrow = ncol(x), ncol = 0,
             dimnames = list(colnames(x), NULL))
    } else {
      free$loadings
    }
    if (ncol(free_scores) > 0) {
      colnames(free_scores) <- colnames(free_loadings) <-
        paste0("free_", seq_len(ncol(free_scores)))
    }

    # The free dimensions are taken from the residual, so they should carry no
    # concept signal. Reporting the largest correlation that survives is what
    # turns that construction into a checkable claim.
    absorbed <- if (ncol(free_scores) > 0 && ncol(scores) > 0) {
      cc <- abs(suppressWarnings(stats::cor(free_scores, scores)))
      cc[!is.finite(cc)] <- 0
      max(cc)
    } else {
      NA_real_
    }

    encodings[[m]] <- list(
      concept_scores = scores,
      free_scores = free_scores,
      free_loadings = free_loadings,
      concept_weights = scored$weights,
      analysis_matrix = x,
      transform = tf$applied,
      stability = if (is.null(free)) data.frame() else free$stability,
      n_free = ncol(free_scores),
      dropped_concepts = scored$dropped
    )
    designs[[m]] <- as.data.frame(SummarizedExperiment::colData(se))

    variance[[m]] <- data.frame(
      modality = m,
      n_samples = nrow(x),
      n_features = ncol(x),
      n_concepts = ncol(scores),
      n_free = ncol(free_scores),
      free_ceiling = ceiling_k,
      concept_share = round(concept_share, 4),
      free_share = round(free_share, 4),
      residual_share = round(1 - concept_share - free_share, 4),
      max_free_concept_correlation = round(absorbed, 4),
      stringsAsFactors = FALSE
    )
    concept_variance[[m]] <- chorale_concept_variance(x, scores, m)
  }

  structure(
    list(
      modalities = modalities,
      concepts = concepts,
      encodings = encodings,
      designs = designs,
      variance = do.call(rbind, variance),
      concept_variance = do.call(rbind, concept_variance),
      control = control,
      seed = seed
    ),
    class = "chorale_encode"
  )
}

#' Concepts the modalities score in common
#'
#' @param encoding A `chorale_encode` object.
#' @returns A character vector of concept names every modality scores.
#' @keywords internal
#' @noRd
chorale_scored_concepts <- function(encoding) {
  Reduce(intersect, lapply(encoding$encodings,
                           function(e) colnames(e$concept_scores)))
}

#' Average standardised level of the features carrying each concept
#'
#' The score is a weighted mean of standardised feature values, so a large
#' concept and a small one are on one scale and a feature mapping to several
#' genes contributes the fraction [chorale_map()] gave it. Scores are
#' standardised afterwards, so an effect on one concept is comparable with an
#' effect on another.
#'
#' @keywords internal
#' @noRd
chorale_concept_scores <- function(x, membership) {
  empty <- matrix(numeric(0), nrow = nrow(x), ncol = 0,
                  dimnames = list(rownames(x), NULL))
  if (is.null(membership) || ncol(membership) == 0) {
    return(list(scores = empty, weights = membership, dropped = character()))
  }
  shared <- intersect(colnames(x), rownames(membership))
  if (length(shared) == 0) {
    return(list(scores = empty, weights = membership[integer(0), , drop = FALSE],
                dropped = colnames(membership)))
  }
  w <- membership[shared, , drop = FALSE]
  mass <- colSums(w)
  keep <- mass > 0
  dropped <- colnames(w)[!keep]
  w <- w[, keep, drop = FALSE]
  mass <- mass[keep]
  if (ncol(w) == 0) {
    return(list(scores = empty, weights = w, dropped = dropped))
  }
  w <- sweep(w, 2, mass, `/`)
  raw <- x[, shared, drop = FALSE] %*% w

  sdv <- apply(raw, 2, stats::sd)
  constant <- !is.finite(sdv) | sdv == 0
  dropped <- c(dropped, colnames(raw)[constant])
  raw <- raw[, !constant, drop = FALSE]
  w <- w[, !constant, drop = FALSE]
  if (ncol(raw) == 0) {
    return(list(scores = empty, weights = w, dropped = dropped))
  }
  scores <- scale(raw)
  attr(scores, "scaled:center") <- NULL
  attr(scores, "scaled:scale") <- NULL
  dimnames(scores) <- list(rownames(x), colnames(raw))
  list(scores = scores, weights = w, dropped = dropped)
}

#' What the concept scores reconstruct of the assay
#'
#' Curated concepts overlap, so their scores are correlated and the coefficients
#' of a regression on them are unstable: which of two overlapping concepts
#' absorbs a shared signal is close to arbitrary. What the overlap does not make
#' arbitrary is the space the scores span, and it is that space, not the
#' individual coefficients, that decides what is left for the free dimensions.
#' The projection is therefore taken onto the span, through a pseudo-inverse
#' that is defined whether or not the scores are of full rank.
#'
#' The residual is then exactly orthogonal to every concept score, so a free
#' dimension cannot carry concept signal by construction rather than by
#' approximation.
#'
#' @keywords internal
#' @noRd
chorale_concept_reconstruction <- function(x, scores) {
  if (ncol(scores) == 0) return(matrix(0, nrow = nrow(x), ncol = ncol(x),
                                       dimnames = dimnames(x)))
  gram <- crossprod(scores)
  scores %*% (chorale_ginv(gram) %*% crossprod(scores, x))
}

#' Share of a modality's variance each concept carries on its own
#' @keywords internal
#' @noRd
chorale_concept_variance <- function(x, scores, modality) {
  if (ncol(scores) == 0) {
    return(data.frame(modality = character(), concept = character(),
                      variance_share = numeric(), stringsAsFactors = FALSE))
  }
  total <- sum(x^2)
  share <- vapply(seq_len(ncol(scores)), function(j) {
    s <- scores[, j]
    denom <- sum(s^2)
    if (!is.finite(denom) || denom == 0 || total <= 0) return(NA_real_)
    b <- crossprod(s, x) / denom
    sum((s %*% b)^2) / total
  }, numeric(1))
  data.frame(modality = modality, concept = colnames(scores),
             variance_share = round(share, 5), stringsAsFactors = FALSE)
}

#' Resolve the free-dimension count requested for every modality
#' @keywords internal
#' @noRd
chorale_free_spec <- function(n_free, modalities) {
  out <- stats::setNames(rep(list("auto"), length(modalities)), modalities)
  if (is.null(n_free)) return(out)
  if (is.null(names(n_free))) {
    if (length(n_free) == 1) {
      out[] <- list(n_free[[1]])
    } else if (length(n_free) == length(modalities)) {
      out[] <- as.list(n_free)
    } else {
      rlang::abort("`n_free` must be one value, one per modality, or named.")
    }
  } else {
    unknown <- setdiff(names(n_free), modalities)
    if (length(unknown) > 0) {
      rlang::abort(paste0("`n_free` names unknown modalities: ",
                          paste(unknown, collapse = ", "), "."))
    }
    out[names(n_free)] <- as.list(n_free)
  }
  for (m in modalities) {
    v <- out[[m]]
    if (identical(v, "auto")) next
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0) {
      rlang::abort("`n_free` must be \"auto\" or a non-negative count.")
    }
  }
  out
}

#' @export
print.chorale_encode <- function(x, ...) {
  cat("<chorale_encode>\n")
  cat("  concepts scored in every modality:",
      length(chorale_scored_concepts(x)), "\n")
  v <- x$variance
  for (i in seq_len(nrow(v))) {
    cat(sprintf(
      "    %-12s %3d concepts, %2d free dimensions; variance %.2f concepts, %.2f free, %.2f left\n",
      v$modality[i], v$n_concepts[i], v$n_free[i], v$concept_share[i],
      v$free_share[i], v$residual_share[i]))
  }
  cat("  largest free-to-concept correlation:",
      signif(max(v$max_free_concept_correlation, na.rm = TRUE), 3), "\n")
  invisible(x)
}
