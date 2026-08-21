# What the modalities share of their designs, and what a design term does to a
# score. The terms are resolved once for a collection so that an effect means
# the same thing in every modality it is estimated in.

#' How a covariate can be compared across modalities
#'
#' One decision, used wherever a covariate is admitted and wherever its terms
#' are built, so the two cannot disagree about the same column.
#'
#' A covariate that is numeric in every modality is continuous and contributes
#' one term, whatever how many distinct values any one modality happens to
#' realise. Deciding otherwise would send it down the categorical path, where
#' levels are matched as text: two modalities recording the same age as `6` and
#' `6.0` would then share no level, and a covariate every modality measures on
#' the same scale would be discarded for a difference of formatting.
#'
#' A covariate that is categorical in every modality contributes one term per
#' shared level beyond a reference. A covariate that is numeric in some
#' modalities and categorical in others cannot be compared at all.
#'
#' @param vals A list holding the covariate as each modality records it.
#' @returns `"continuous"`, `"categorical"`, or `NA_character_` where the
#'   modalities disagree about its type.
#' @keywords internal
#' @noRd
chorale_covariate_kind <- function(vals) {
  if (!length(vals)) return(NA_character_)
  numeric_all <- all(vapply(vals, is.numeric, logical(1)))
  categorical_all <- all(vapply(vals, function(v) !is.numeric(v), logical(1)))
  if (numeric_all) return("continuous")
  if (categorical_all) return("categorical")
  NA_character_
}

#' Terms a design profile can carry in every modality
#'
#' A profile entry is only comparable across modalities if it means the same
#' thing in each, so the terms are fixed once for the whole collection rather
#' than per pair. A covariate that is numeric everywhere contributes one term.
#' A categorical covariate contributes one term per level beyond a reference,
#' using only the levels that every design realises, so a level observed in one
#' cohort and not another cannot enter the comparison.
#'
#' @param designs A list of design tables.
#' @param covariates Candidate covariate columns.
#'
#' @returns A named list, one entry per usable covariate: `NA_character_` for a
#'   continuous covariate, otherwise the shared levels with the reference first.
#'   A covariate that cannot be made comparable has no entry.
#' @keywords internal
#' @noRd
chorale_profile_levels <- function(designs, covariates) {
  out <- list()
  for (cv in covariates) {
    if (!all(vapply(designs, function(d) cv %in% colnames(d), logical(1)))) next
    vals <- lapply(designs, function(d) d[[cv]])
    kind <- chorale_covariate_kind(vals)
    if (is.na(kind)) next
    if (identical(kind, "continuous")) {
      out[[cv]] <- NA_character_
      next
    }
    lv <- Reduce(intersect, lapply(vals, function(v) {
      unique(as.character(stats::na.omit(v)))
    }))
    lv <- sort(lv)
    if (length(lv) < 2) next
    out[[cv]] <- lv
  }
  out
}

#' Column names a set of profile terms produces
#' @keywords internal
#' @noRd
chorale_profile_terms <- function(levels) {
  unlist(lapply(names(levels), function(cv) {
    lv <- levels[[cv]]
    if (length(lv) == 1 && is.na(lv)) return(cv)
    paste0(cv, "=", lv[-1])
  }), use.names = FALSE)
}

#' Resolve a common adjusted design signature
#' @keywords internal
#' @noRd
chorale_resolve_signature <- function(designs, phenotype_column = "phenotype",
                                      phenotype_reference = "control",
                                      profile_covariates = NULL) {
  modalities <- names(designs)
  candidates <- profile_covariates %||% chorale_candidate_covariates(designs)
  candidates <- unique(c(phenotype_column, candidates))
  all_columns <- unique(unlist(lapply(designs, colnames)))
  excluded <- list()
  shared <- character()

  for (cv in candidates) {
    absent <- modalities[!vapply(designs, function(d) cv %in% names(d), logical(1))]
    if (length(absent)) {
      excluded[[length(excluded) + 1L]] <- data.frame(
        covariate = cv, reason = "absent from a modality",
        detail = paste(absent, collapse = ", "), stringsAsFactors = FALSE)
      next
    }
    vals <- lapply(designs, `[[`, cv)
    missing_share <- vapply(vals, function(v) mean(is.na(v)), numeric(1))
    if (any(missing_share > 0.5)) {
      excluded[[length(excluded) + 1L]] <- data.frame(
        covariate = cv, reason = "excessive missingness",
        detail = paste0(modalities[missing_share > 0.5], " (",
                        round(100 * missing_share[missing_share > 0.5]), "%)",
                        collapse = ", "), stringsAsFactors = FALSE)
      next
    }
    varying <- vapply(vals, function(v) length(unique(stats::na.omit(v))) >= 2L,
                      logical(1))
    if (!all(varying)) {
      excluded[[length(excluded) + 1L]] <- data.frame(
        covariate = cv, reason = "constant or entirely missing",
        detail = paste(modalities[!varying], collapse = ", "),
        stringsAsFactors = FALSE)
      next
    }
    kind <- chorale_covariate_kind(vals)
    numeric_all <- identical(kind, "continuous")
    categorical_all <- identical(kind, "categorical")
    if (is.na(kind)) {
      excluded[[length(excluded) + 1L]] <- data.frame(
        covariate = cv, reason = "incompatible types", detail = "",
        stringsAsFactors = FALSE)
      next
    }
    if (categorical_all) {
      level_sets <- lapply(vals, function(v) sort(unique(as.character(
        stats::na.omit(v)))))
      common <- Reduce(intersect, level_sets)
      if (length(common) < 2L) {
        excluded[[length(excluded) + 1L]] <- data.frame(
          covariate = cv, reason = "fewer than two shared levels", detail = "",
          stringsAsFactors = FALSE)
        next
      }
      compatible <- all(vapply(level_sets, identical, logical(1), common))
      if (!compatible) {
        if (identical(cv, phenotype_column)) {
          rlang::abort(
            paste0("Phenotype levels must match across modalities; found: ",
                   paste(vapply(level_sets, paste, character(1), collapse = "/"),
                         collapse = "; "), "."),
            class = "chorale_incompatible_phenotype_levels")
        }
        excluded[[length(excluded) + 1L]] <- data.frame(
          covariate = cv, reason = "incompatible levels",
          detail = paste(vapply(level_sets, paste, character(1), collapse = "/"),
                         collapse = "; "), stringsAsFactors = FALSE)
        next
      }
    }
    shared <- c(shared, cv)
  }

  if (!phenotype_column %in% shared) {
    rlang::abort(
      paste0("Every modality must contain an estimable `", phenotype_column,
             "` contrast with compatible levels."),
      class = "chorale_missing_phenotype"
    )
  }

  levels <- chorale_profile_levels(designs, shared)

  # A covariate is only retained where its terms could be built. Admitting one
  # here and failing to resolve its levels would leave a covariate carried into
  # the model matrix with nothing to build a block from.
  unresolved <- setdiff(shared, names(levels))
  if (length(unresolved)) {
    if (phenotype_column %in% unresolved) {
      rlang::abort(
        paste0("`", phenotype_column, "` has no comparable terms across the ",
               "modalities; it must be categorical with at least two shared ",
               "levels."),
        class = "chorale_invalid_phenotype")
    }
    excluded[[length(excluded) + 1L]] <- data.frame(
      covariate = unresolved, reason = "no comparable terms across modalities",
      detail = "", stringsAsFactors = FALSE)
    shared <- setdiff(shared, unresolved)
  }

  phenotype_levels <- levels[[phenotype_column]]
  if (is.null(phenotype_levels) || anyNA(phenotype_levels)) {
    rlang::abort("The mandatory phenotype must be categorical.",
                 class = "chorale_invalid_phenotype")
  }
  ref_at <- match(tolower(phenotype_reference), tolower(phenotype_levels))
  if (is.na(ref_at)) {
    rlang::abort(
      paste0("Phenotype reference `", phenotype_reference,
             "` was not found in the shared levels: ",
             paste(phenotype_levels, collapse = ", "), "."),
      class = "chorale_missing_phenotype_reference"
    )
  }
  levels[[phenotype_column]] <- c(phenotype_levels[ref_at],
                                   phenotype_levels[-ref_at])

  # Add secondary covariates only when their complete block is estimable in
  # every modality after the blocks already retained. This catches duplicate
  # encodings (for example a continuous age and bins derived from it), nested
  # batches and other rank deficiencies without discarding the phenotype model.
  selected <- phenotype_column
  for (cv in setdiff(shared, phenotype_column)) {
    trial <- c(selected, cv)
    trial_spec <- list(covariates = trial, levels = levels[trial])
    full_rank <- all(vapply(designs, function(d) {
      x <- chorale_signature_matrix(d, trial_spec)$x
      ok <- stats::complete.cases(x)
      sum(ok) > ncol(x) && qr(x[ok, , drop = FALSE])$rank == ncol(x)
    }, logical(1)))
    if (full_rank) {
      selected <- trial
    } else {
      excluded[[length(excluded) + 1L]] <- data.frame(
        covariate = cv, reason = "not jointly estimable",
        detail = "rank deficient after shared covariates",
        stringsAsFactors = FALSE)
    }
  }
  shared <- selected
  levels <- levels[shared]

  list(
    phenotype = phenotype_column,
    covariates = shared,
    secondary = setdiff(shared, phenotype_column),
    levels = levels,
    terms = chorale_profile_terms(levels),
    excluded = if (length(excluded)) do.call(rbind, excluded) else
      data.frame(covariate = character(), reason = character(),
                 detail = character()),
    unused_columns = setdiff(all_columns, c(shared, "sample_id"))
  )
}

#' Build a common model matrix without formula-name ambiguity
#' @keywords internal
#' @noRd
chorale_signature_matrix <- function(design, spec, include = spec$covariates) {
  n <- nrow(design)
  missing_levels <- setdiff(intersect(spec$covariates, include), names(spec$levels))
  if (length(missing_levels)) {
    rlang::abort(
      paste0("No terms were resolved for: ", paste(missing_levels, collapse = ", "),
             ". A covariate carried in the signature must have an entry in its levels."),
      class = "chorale_unresolved_covariate")
  }
  blocks <- list()
  term_covariate <- character()
  for (cv in intersect(spec$covariates, include)) {
    lv <- spec$levels[[cv]]
    if (length(lv) == 1L && is.na(lv)) {
      v <- suppressWarnings(as.numeric(design[[cv]]))
      z <- as.numeric(scale(v))
      blocks[[cv]] <- matrix(z, ncol = 1L, dimnames = list(NULL, cv))
      term_covariate <- c(term_covariate, cv)
    } else {
      v <- as.character(design[[cv]])
      block <- vapply(lv[-1L], function(one) as.numeric(v == one), numeric(n))
      if (is.null(dim(block))) block <- matrix(block, ncol = 1L)
      colnames(block) <- paste0(cv, "=", lv[-1L])
      block[is.na(v), ] <- NA_real_
      blocks[[cv]] <- block
      term_covariate <- c(term_covariate, rep(cv, ncol(block)))
    }
  }
  x <- if (length(blocks)) do.call(cbind, blocks) else matrix(numeric(), n, 0L)
  list(x = cbind(`(Intercept)` = 1, x), blocks = blocks,
       term_covariate = stats::setNames(term_covariate, colnames(x)))
}

#' Adjusted effects and their standard errors for every factor
#' @keywords internal
#' @noRd
chorale_adjusted_profile <- function(scores, design, spec) {
  design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]
  mm <- chorale_signature_matrix(design, spec)
  terms <- colnames(mm$x)[-1L]
  effects <- se <- matrix(NA_real_, nrow = ncol(scores), ncol = length(terms),
                          dimnames = list(colnames(scores), terms))
  estimable <- matrix(FALSE, nrow = ncol(scores), ncol = length(terms),
                      dimnames = dimnames(effects))
  covariance <- vector("list", ncol(scores))
  names(covariance) <- colnames(scores)

  common_ok <- stats::complete.cases(mm$x) & stats::complete.cases(scores)
  common_x <- mm$x[common_ok, , drop = FALSE]
  if (sum(common_ok) > ncol(common_x) && qr(common_x)$rank == ncol(common_x)) {
    fit <- stats::lm.fit(common_x, scores[common_ok, , drop = FALSE])
    inv <- chol2inv(qr.R(fit$qr))
    # A single response column comes back as a vector rather than a matrix, so
    # both are restored to matrices before anything indexes them by column.
    residuals <- as.matrix(fit$residuals)
    coefficients <- as.matrix(fit$coefficients)
    sigma2 <- colSums(residuals^2) / fit$df.residual
    effects[,] <- t(coefficients[-1L, , drop = FALSE])
    for (j in seq_len(ncol(scores))) {
      v <- sigma2[j] * inv[-1L, -1L, drop = FALSE]
      covariance[[j]] <- v
      se[j, ] <- sqrt(pmax(diag(v), 0))
      estimable[j, ] <- is.finite(effects[j, ]) & is.finite(se[j, ]) & se[j, ] > 0
    }
  } else for (j in seq_len(ncol(scores))) {
    y <- scores[, j]
    ok <- is.finite(y) & stats::complete.cases(mm$x)
    x <- mm$x[ok, , drop = FALSE]
    yy <- y[ok]
    if (nrow(x) <= ncol(x) || qr(x)$rank < ncol(x)) next
    fit <- stats::lm.fit(x, yy)
    rss <- sum(fit$residuals^2)
    sigma2 <- rss / fit$df.residual
    inv <- chol2inv(qr.R(fit$qr))
    vc <- sigma2 * inv
    b <- fit$coefficients[-1L]
    v <- vc[-1L, -1L, drop = FALSE]
    effects[j, ] <- b
    se[j, ] <- sqrt(pmax(diag(v), 0))
    estimable[j, ] <- is.finite(b) & is.finite(se[j, ]) & se[j, ] > 0
    covariance[[j]] <- v
  }

  z <- effects / se
  z[!is.finite(z)] <- 0
  list(effects = effects, se = se, z = z, estimable = estimable,
       covariance = covariance, blocks = mm$blocks,
       term_covariate = mm$term_covariate)
}
