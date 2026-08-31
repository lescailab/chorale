#' One latent state estimated from every modality at once
#'
#' [chorale_concept_evidence()] asks each modality the same question separately
#' and combines the answers. That is a meta-analysis: the modality is the unit
#' of estimation, so a modality that cannot estimate a concept's effect
#' contributes nothing to it, and the object being estimated is a scalar each
#' modality could in principle have estimated alone.
#'
#' This function estimates a different object. After encoding, every modality is
#' described on the same vocabulary, and each score is a weighted mean of
#' standardised feature values, so the score matrices share their columns and
#' are on one scale. They can therefore be stacked by row into a single matrix
#' of every sample in the collection against every concept, and decomposed once.
#' The concept loadings that come out are common to the whole collection by
#' construction rather than matched across modalities afterwards.
#'
#' The fitted directions live in concept space, whose dimension is the size of
#' the vocabulary. Stacking increases the number and diversity of observations
#' available to estimate those directions and permits concepts with incomplete
#' modality coverage to contribute wherever they are observed. A direction is
#' not guaranteed to be unique or biologically identified merely because the
#' vocabulary is larger than one modality's sample count; interpretation rests
#' on stability, loadings, phenotype evidence and held-out transfer.
#'
#' @section What is removed before stacking:
#' Two things would otherwise dominate a decomposition of stacked modalities.
#' The first is the difference between modalities themselves, which is removed
#' by centring and scaling each concept within each modality, so only variation
#' between samples of the same modality survives. The second is whatever a
#' modality's own nuisance structure contributes: study, acquisition batch,
#' colonisation, region. Those are named per modality in `nuisance` and
#' regressed out of that modality's scores before stacking. They need not be
#' shared with any other modality, which is what makes them usable here and not
#' in a model the modalities hold in common.
#'
#' @section Concepts a modality cannot express:
#' Coverage differs between modalities, and a concept a modality cannot be
#' scored on is missing rather than zero. Those entries carry no weight in the
#' fit, so a modality contributes to the components only through the concepts it
#' measures. Restricting the vocabulary to concepts every modality expresses
#' would discard most of it and reinstate the narrowest modality as the limit on
#' the whole collection.
#'
#' @param encoding A `chorale_encode` object.
#' @param n_components Number of components. `NULL`, the default, selects by
#'   reproducibility under subsampling with [chorale_select_factors()], held to
#'   the same settings the encoder selects free dimensions under:
#'   `n_select_init`, `n_subsample`, `subsample_fraction` and `reproducibility`
#'   from `control`. The selector needs a complete matrix, so it is run on the
#'   stacked scores with the unobserved entries filled with zero, while the fit
#'   itself uses the observed entries alone. The reproducibility of the selected
#'   rank is therefore a property of the coverage-filled matrix, and a
#'   collection whose coverage is sparse should be read with that in mind or
#'   given `n_components` directly.
#' @param max_components Ceiling on the number selected.
#' @param nuisance Optional named list, one entry per modality, giving design
#'   columns of that modality to regress out of its scores before stacking. A
#'   column naming the phenotype is refused.
#' @param control A [chorale_control()] object.
#' @param seed Integer seed.
#' @param tol,max_iter Convergence tolerance and iteration cap for the weighted
#'   fit.
#' @param ... Named overrides applied to `control`.
#'
#' @returns An object of class `chorale_joint_state` with `scores`, one row per
#'   sample of the whole collection and one column per component; `loadings`,
#'   one row per concept and one column per component; `design`, the stacked
#'   design carrying a `modality` column; `variance`, the share of observed
#'   variance each component carries; and `coverage`, how many concepts each
#'   modality contributed.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 80,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' ids <- sprintf("feature_%05d", seq_len(80))
#' sim$modalities <- lapply(sim$modalities, function(m) {
#'   rownames(m) <- ids
#'   m
#' })
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' cc <- chorale_concepts(containers, list(one = ids[1:20], two = ids[21:45]),
#'                        min_features = 5)
#' enc <- chorale_encode(containers, cc, n_free = 1, n_init = 2)
#' chorale_joint_state(enc, n_components = 1)
chorale_joint_state <- function(encoding, n_components = NULL,
                                max_components = 10L, nuisance = NULL,
                                control = chorale_control(), seed = 1L,
                                tol = 1e-6, max_iter = 200L, ...) {
  if (!inherits(encoding, "chorale_encode")) {
    rlang::abort("`encoding` must be a chorale_encode object.")
  }
  control <- chorale_merge_control(control, list(...))
  modalities <- encoding$modalities
  designs <- encoding$designs

  stacked <- chorale_stack_concept_scores(encoding, nuisance,
                                          control$phenotype_column)
  z <- stacked$scores
  observed <- !is.na(z)
  if (sum(observed) == 0) {
    rlang::abort("No modality carries a concept score to stack.")
  }

  ceiling_k <- min(as.integer(max_components), nrow(z) - 1L, ncol(z))
  if (ceiling_k < 1L) {
    rlang::abort("The stacked collection is too small to carry a component.")
  }
  selection <- NULL
  if (is.null(n_components)) {
    # The selector needs a complete matrix, so a concept a modality cannot
    # express is filled with zero here, which is that concept's within-modality
    # mean. The fit below carries no weight on those entries. Coverage
    # therefore enters the chosen rank and nothing else, and it enters under the
    # settings the encoder selects free dimensions with rather than under a
    # second set of defaults.
    filled <- z
    filled[!observed] <- 0
    selection <- chorale_select_factors(
      filled, max_factors = ceiling_k, n_init = control$n_select_init,
      n_subsample = control$n_subsample,
      subsample_fraction = control$subsample_fraction,
      threshold = control$reproducibility, seed = seed)
    n_components <- attr(selection, "selected") %||% 0L
  }
  n_components <- min(as.integer(n_components), ceiling_k)
  if (n_components < 1L) {
    return(chorale_empty_joint_state(stacked, selection))
  }

  fit <- chorale_weighted_lowrank(z, observed, n_components, seed = seed,
                                  tol = tol, max_iter = max_iter)

  total <- sum(z[observed]^2)
  variance <- vapply(seq_len(n_components), function(r) {
    approximation <- fit$scores[, r, drop = FALSE] %*%
      t(fit$loadings[, r, drop = FALSE])
    if (total > 0) sum(approximation[observed]^2) / total else NA_real_
  }, numeric(1))

  structure(
    list(
      scores = fit$scores,
      loadings = fit$loadings,
      design = stacked$design,
      modality = stacked$modality,
      variance = data.frame(component = colnames(fit$scores),
                            share = round(variance, 4),
                            stringsAsFactors = FALSE),
      coverage = stacked$coverage,
      nuisance = stacked$nuisance,
      n_components = n_components,
      selection = selection,
      iterations = fit$iterations,
      converged = fit$converged
    ),
    class = "chorale_joint_state"
  )
}

#' Stack every modality's concept scores into one matrix
#'
#' @param encoding A `chorale_encode` object.
#' @param nuisance Named list of layer-local design columns to regress out.
#' @param phenotype_column The phenotype column, which nuisance may not name.
#' @returns A list with the stacked score matrix, the stacked design, the
#'   modality label of each row and a coverage table.
#' @keywords internal
#' @noRd
chorale_stack_concept_scores <- function(encoding, nuisance,
                                         phenotype_column) {
  modalities <- encoding$modalities
  designs <- encoding$designs
  nuisance <- nuisance %||% list()
  if (length(nuisance) > 0 && is.null(names(nuisance))) {
    rlang::abort("`nuisance` must be a named list, one entry per modality.")
  }
  unknown <- setdiff(names(nuisance), modalities)
  if (length(unknown) > 0) {
    rlang::abort(paste0("`nuisance` names modalities the encoding does not ",
                        "carry: ", paste(unknown, collapse = ", "), "."))
  }

  vocabulary <- unique(unlist(lapply(encoding$encodings, function(e)
    colnames(e$concept_scores))))
  if (length(vocabulary) == 0) {
    rlang::abort("The encoding carries no concept scores.")
  }

  blocks <- list()
  design_rows <- list()
  labels <- character()
  coverage <- list()
  applied <- list()

  for (m in modalities) {
    scores <- encoding$encodings[[m]]$concept_scores
    if (is.null(scores) || ncol(scores) == 0) next
    design <- designs[[m]]
    design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]

    columns <- nuisance[[m]]
    if (length(columns) > 0) {
      if (phenotype_column %in% columns) {
        rlang::abort(paste0("`nuisance` names the phenotype for modality `", m,
                            "`, which would remove the contrast being tested."),
                     class = "chorale_nuisance_is_phenotype")
      }
      scores <- chorale_residualise_scores(scores, design, columns, m)
    }
    applied[[m]] <- columns %||% character()

    # Centring and scaling within the modality is what makes the rows
    # comparable: what survives is a sample's position among samples measured
    # the same way, not the level a platform happens to read at.
    scores <- scale(scores)
    scores[!is.finite(scores)] <- NA_real_

    full <- matrix(NA_real_, nrow = nrow(scores), ncol = length(vocabulary),
                   dimnames = list(rownames(scores), vocabulary))
    full[, colnames(scores)] <- scores
    blocks[[m]] <- full
    labels <- c(labels, rep(m, nrow(full)))
    design_rows[[m]] <- data.frame(sample_id = rownames(full), modality = m,
                                   stringsAsFactors = FALSE)
    coverage[[m]] <- data.frame(
      modality = m, n_samples = nrow(full),
      n_concepts = sum(colSums(!is.na(full)) > 0),
      stringsAsFactors = FALSE)
  }

  if (length(blocks) == 0) rlang::abort("No modality carries concept scores.")

  z <- do.call(rbind, blocks)
  rownames(z) <- paste(labels, unlist(lapply(blocks, rownames)), sep = ":")

  shared <- Reduce(intersect, lapply(designs, colnames))
  stacked_design <- do.call(rbind, lapply(modalities, function(m) {
    if (is.null(blocks[[m]])) return(NULL)
    d <- designs[[m]]
    d <- d[match(rownames(blocks[[m]]), d$sample_id), shared, drop = FALSE]
    d$modality <- m
    d$sample_id <- paste(m, d$sample_id, sep = ":")
    d
  }))
  rownames(stacked_design) <- NULL

  list(scores = z, design = stacked_design, modality = labels,
       coverage = do.call(rbind, coverage), nuisance = applied,
       vocabulary = vocabulary)
}

#' Regress a modality's own nuisance covariates out of its scores
#'
#' The covariates need not be shared with any other modality, which is the
#' point: study, acquisition batch, colonisation and region are recorded by the
#' modality that has them and by no other, so a model the modalities hold in
#' common cannot adjust for them and a model local to one modality can.
#'
#' @keywords internal
#' @noRd
chorale_residualise_scores <- function(scores, design, columns, modality) {
  missing <- setdiff(columns, colnames(design))
  if (length(missing) > 0) {
    rlang::abort(paste0("Modality `", modality, "` has no design column(s) ",
                        paste(missing, collapse = ", "), "."))
  }
  parts <- lapply(columns, function(cv) {
    v <- design[[cv]]
    if (is.character(v) || is.factor(v)) factor(v) else v
  })
  names(parts) <- columns
  frame <- as.data.frame(parts, stringsAsFactors = FALSE)
  usable <- vapply(frame, function(v) {
    length(unique(v[!is.na(v)])) > 1
  }, logical(1))
  if (!any(usable)) return(scores)
  frame <- frame[, usable, drop = FALSE]

  ok <- stats::complete.cases(frame)
  x <- stats::model.matrix(~., data = frame[ok, , drop = FALSE])
  if (sum(ok) <= ncol(x) || qr(x)$rank < ncol(x)) return(scores)

  out <- scores
  fit <- stats::lm.fit(x, scores[ok, , drop = FALSE])
  # A sample the nuisance model cannot use keeps its own value, so it is
  # neither dropped nor adjusted by a model it did not enter.
  out[ok, ] <- as.matrix(fit$residuals)
  out
}

#' Low-rank fit over the observed entries only
#'
#' Alternating ridge-regularised least squares. Entries a modality cannot
#' supply carry no weight, so a component is built from what each modality
#' actually measures rather than from an imputed value standing in for it.
#'
#' @keywords internal
#' @noRd
chorale_weighted_lowrank <- function(z, observed, k, seed = 1L, tol = 1e-6,
                                     max_iter = 200L, ridge = 1e-6) {
  # The starting point is the decomposition of the matrix with missing entries
  # set to the concept mean, which after within-modality scaling is zero. That
  # is an initialisation and nothing else: every row and column carrying an
  # observation is re-estimated below from its observed entries alone, so no
  # imputed value survives into the fit.
  filled <- z
  filled[!observed] <- 0
  set.seed(seed)
  start <- svd(filled, nu = k, nv = k)
  scores <- start$u %*% diag(start$d[seq_len(k)], nrow = k)
  loadings <- start$v

  rss_of <- function(s, l) {
    approximation <- s %*% t(l)
    sum((z[observed] - approximation[observed])^2)
  }
  previous <- rss_of(scores, loadings)
  converged <- FALSE
  iterations <- 0L

  for (it in seq_len(max_iter)) {
    iterations <- it
    for (i in seq_len(nrow(z))) {
      j <- which(observed[i, ])
      # Only a row with nothing observed is skipped. The ridge makes the normal
      # equations solvable from a single observation, giving the smallest score
      # consistent with what the row actually carries; requiring as many
      # observations as components would instead leave a sparsely covered row
      # holding its initialisation, which is the imputed value this fit exists
      # to avoid.
      if (length(j) == 0L) next
      l <- loadings[j, , drop = FALSE]
      scores[i, ] <- solve(crossprod(l) + ridge * diag(k),
                           crossprod(l, z[i, j]))
    }
    for (j in seq_len(ncol(z))) {
      i <- which(observed[, j])
      if (length(i) == 0L) next
      s <- scores[i, , drop = FALSE]
      loadings[j, ] <- solve(crossprod(s) + ridge * diag(k),
                             crossprod(s, z[i, j]))
    }
    current <- rss_of(scores, loadings)
    if (is.finite(previous) && previous > 0 &&
        abs(previous - current) / previous < tol) {
      converged <- TRUE
      break
    }
    previous <- current
  }

  # The factorisation is identified only up to a rotation of the k components,
  # so they are returned in the orthogonal basis that orders them by the
  # variance they carry, which is what makes one component comparable across
  # runs. The rotation is read off a k-by-k matrix rather than off the
  # reconstruction, which is never formed.
  qs <- qr(scores)
  ql <- qr(loadings)
  middle <- svd(qr.R(qs) %*% t(qr.R(ql)))
  keep <- seq_len(k)
  scores <- qr.Q(qs) %*% middle$u[, keep, drop = FALSE] %*%
    diag(middle$d[keep], nrow = k)
  loadings <- qr.Q(ql) %*% middle$v[, keep, drop = FALSE]

  dimnames(scores) <- list(rownames(z), sprintf("joint_%02d", keep))
  dimnames(loadings) <- list(colnames(z), sprintf("joint_%02d", keep))
  list(scores = scores, loadings = loadings, iterations = iterations,
       converged = converged)
}

#' @keywords internal
#' @noRd
chorale_empty_joint_state <- function(stacked, selection) {
  structure(
    list(scores = matrix(numeric(0), nrow = nrow(stacked$scores), ncol = 0,
                         dimnames = list(rownames(stacked$scores), NULL)),
         loadings = matrix(numeric(0), nrow = ncol(stacked$scores), ncol = 0,
                           dimnames = list(colnames(stacked$scores), NULL)),
         design = stacked$design, modality = stacked$modality,
         variance = data.frame(component = character(), share = numeric()),
         coverage = stacked$coverage, nuisance = stacked$nuisance,
         n_components = 0L, selection = selection,
         iterations = 0L, converged = TRUE),
    class = "chorale_joint_state")
}

#' @export
print.chorale_joint_state <- function(x, ...) {
  cat("<chorale_joint_state>\n")
  cat("  samples:    ", nrow(x$scores), " across ",
      length(unique(x$modality)), " modalities\n", sep = "")
  cat("  concepts:   ", nrow(x$loadings), "\n", sep = "")
  cat("  components: ", x$n_components, "\n", sep = "")
  if (nrow(x$variance) > 0) {
    cat("  variance:   ",
        paste(sprintf("%s %.3f", x$variance$component, x$variance$share),
              collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
