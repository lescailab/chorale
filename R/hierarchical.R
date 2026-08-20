# Internal machinery for phenotype-led, covariate-refined matching.

#' Deterministic task mapping for repeated resampling
#' @keywords internal
#' @noRd
chorale_deterministic_lapply <- function(x, fun, n_cores = 1L) {
  if (n_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(x, fun, mc.cores = n_cores, mc.preschedule = TRUE,
                       mc.set.seed = FALSE)
  } else {
    lapply(x, fun)
  }
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
    numeric_all <- all(vapply(vals, is.numeric, logical(1)))
    categorical_all <- all(vapply(vals, function(v) !is.numeric(v), logical(1)))
    if (!numeric_all && !categorical_all) {
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
    sigma2 <- colSums(fit$residuals^2) / fit$df.residual
    effects[,] <- t(fit$coefficients[-1L, , drop = FALSE])
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

#' Equal-block secondary signature
#' @keywords internal
#' @noRd
chorale_secondary_signature <- function(profile, spec) {
  z <- profile$z
  if (length(spec$secondary) == 0L) return(z[, integer(), drop = FALSE])
  out <- list()
  for (cv in spec$secondary) {
    at <- which(profile$term_covariate == cv)
    if (!length(at)) next
    block <- z[, at, drop = FALSE] / sqrt(length(at))
    colnames(block) <- paste0("secondary:", colnames(block))
    out[[cv]] <- block
  }
  if (length(out)) do.call(cbind, out) else z[, integer(), drop = FALSE]
}

#' Uncertainty-weighted compatibility for one secondary covariate block
#' @keywords internal
#' @noRd
chorale_secondary_block_affinity <- function(a, b, covariate,
                                             orientation = 1) {
  ia <- which(a$term_covariate == covariate)
  ib <- which(b$term_covariate == covariate)
  if (!length(ia) || length(ia) != length(ib)) {
    return(matrix(0, nrow(a$effects), nrow(b$effects)))
  }
  out <- matrix(0, nrow(a$effects), nrow(b$effects),
                dimnames = list(rownames(a$effects), rownames(b$effects)))
  for (r in seq_len(nrow(out))) for (s in seq_len(ncol(out))) {
    ea <- a$effects[r, ia]
    eb <- b$effects[s, ib] * orientation[r, s]
    va <- a$se[r, ia]^2
    vb <- b$se[s, ib]^2
    ok <- is.finite(ea) & is.finite(eb) & is.finite(va) & is.finite(vb) &
      va > 0 & vb > 0
    if (!any(ok)) next
    q <- mean((ea[ok] - eb[ok])^2 / (va[ok] + vb[ok]))
    strength <- sqrt(mean(ea[ok]^2 / va[ok]) * mean(eb[ok]^2 / vb[ok]))
    out[r, s] <- (1 - exp(-pmax(strength, 0) / 2)) * exp(-q / 20)
  }
  out
}

#' Largest block-weighted secondary similarity in a factor search
#' @keywords internal
#' @noRd
chorale_secondary_max <- function(profiles, spec, covariate) {
  mods <- names(profiles)
  best <- 0
  for (i in seq_along(mods)) for (j in seq_along(mods)) {
    if (j <= i) next
    orient <- matrix(1, nrow(profiles[[mods[i]]]$effects),
                     nrow(profiles[[mods[j]]]$effects))
    sim <- chorale_secondary_block_affinity(profiles[[mods[i]]],
                                             profiles[[mods[j]]], covariate,
                                             orient)
    best <- max(best, sim, na.rm = TRUE)
  }
  best
}

#' Pairwise phenotype affinity with sign and compatibility diagnostics
#' @keywords internal
#' @noRd
chorale_phenotype_affinity <- function(a, b, phenotype_terms) {
  na <- nrow(a$effects)
  nb <- nrow(b$effects)
  affinity <- loss <- matrix(0, na, nb,
                             dimnames = list(rownames(a$effects),
                                             rownames(b$effects)))
  orientation <- matrix(1, na, nb, dimnames = dimnames(affinity))
  signal <- matrix(0, na, nb, dimnames = dimnames(affinity))
  for (i in seq_len(na)) for (j in seq_len(nb)) {
    ea <- a$effects[i, phenotype_terms]
    eb <- b$effects[j, phenotype_terms]
    va <- a$se[i, phenotype_terms]^2
    vb <- b$se[j, phenotype_terms]^2
    ok <- is.finite(ea) & is.finite(eb) & is.finite(va) & is.finite(vb) &
      va > 0 & vb > 0
    if (!any(ok)) next
    denom <- va[ok] + vb[ok]
    same <- mean((ea[ok] - eb[ok])^2 / denom)
    opposite <- mean((ea[ok] + eb[ok])^2 / denom)
    use_opposite <- opposite < same
    q <- min(same, opposite)
    za <- mean((ea[ok]^2) / va[ok])
    zb <- mean((eb[ok]^2) / vb[ok])
    strength <- sqrt(pmax(za, 0) * pmax(zb, 0))
    # Signal strength supplies power against the phenotype-null permutation;
    # the compatibility penalty is deliberately softer than a test of exact
    # coefficient equality because separately standardised modalities need not
    # express the same biological effect at identical magnitude.
    affinity[i, j] <- strength * exp(-q / 20)
    loss[i, j] <- q
    signal[i, j] <- strength
    orientation[i, j] <- if (use_opposite) -1 else 1
  }
  list(affinity = affinity, loss = loss, signal = signal,
       orientation = orientation)
}

#' Build phenotype and secondary pairwise blocks
#' @keywords internal
#' @noRd
chorale_hierarchical_blocks <- function(profiles, spec,
                                        ambiguity_level = 0.95) {
  mods <- names(profiles)
  phenotype_terms <- names(profiles[[1]]$term_covariate)[
    profiles[[1]]$term_covariate == spec$phenotype]
  primary <- secondary <- final <- diagnostics <- list()
  cutoff <- stats::qchisq(ambiguity_level, df = max(1L, length(phenotype_terms)))
  for (i in seq_along(mods)) for (j in seq_along(mods)) {
    if (j <= i) next
    key <- paste(mods[i], mods[j], sep = "|")
    ph <- chorale_phenotype_affinity(profiles[[i]], profiles[[j]], phenotype_terms)
    sec <- matrix(0, nrow(ph$affinity), ncol(ph$affinity),
                  dimnames = dimnames(ph$affinity))
    if (length(spec$secondary)) {
      blocks_by_covariate <- lapply(spec$secondary, function(cv) {
        chorale_secondary_block_affinity(profiles[[i]], profiles[[j]], cv,
                                         ph$orientation)
      })
      sec <- Reduce(`+`, blocks_by_covariate) / length(blocks_by_covariate)
    }
    # Secondary evidence can only order candidates whose phenotype loss is not
    # distinguishable from the best loss for that factor. It cannot promote an
    # incompatible candidate into the phenotype-compatible set.
    candidate <- matrix(FALSE, nrow(sec), ncol(sec), dimnames = dimnames(sec))
    for (r in seq_len(nrow(candidate))) {
      best <- min(ph$loss[r, ], na.rm = TRUE)
      candidate[r, ] <- is.finite(ph$loss[r, ]) & ph$loss[r, ] <= best + cutoff
    }
    sec01 <- matrix((pmax(-1, pmin(1, as.numeric(sec))) + 1) / 2,
                    nrow = nrow(sec), ncol = ncol(sec),
                    dimnames = dimnames(sec))
    combined <- ph$affinity
    if (length(spec$secondary) > 0L && any(abs(sec) > 0, na.rm = TRUE)) {
      for (r in seq_len(nrow(combined))) {
        cand <- which(candidate[r, ])
        if (length(cand)) combined[r, cand] <- 1 + sec01[r, cand]
        non <- which(!candidate[r, ])
        if (length(non)) combined[r, non] <- pmin(ph$affinity[r, non], 0.999999)
      }
    }
    primary[[key]] <- ph$affinity
    secondary[[key]] <- sec
    final[[key]] <- combined
    diagnostics[[key]] <- c(ph, list(candidate = candidate))
  }
  list(primary = primary, secondary = secondary, final = final,
       diagnostics = diagnostics, phenotype_terms = phenotype_terms)
}

#' Replace analytic candidate sets with bootstrap candidate sets
#' @keywords internal
#' @noRd
chorale_bootstrap_candidates <- function(blocks, bootstrap_blocks,
                                         ambiguity_level = 0.95) {
  if (!length(bootstrap_blocks)) return(blocks)
  alpha <- 1 - ambiguity_level
  for (key in names(blocks$diagnostics)) {
    observed <- blocks$diagnostics[[key]]$loss
    draws <- lapply(bootstrap_blocks, function(x) x$diagnostics[[key]]$loss)
    draws <- draws[vapply(draws, function(x) identical(dim(x), dim(observed)),
                          logical(1))]
    if (!length(draws)) next
    candidate <- matrix(FALSE, nrow(observed), ncol(observed),
                        dimnames = dimnames(observed))
    margin_lower <- margin_upper <- matrix(NA_real_, nrow(observed), ncol(observed),
                                            dimnames = dimnames(observed))
    for (i in seq_len(nrow(observed))) {
      for (j in seq_len(ncol(observed))) {
        delta <- vapply(draws, function(x) {
          row <- x[i, ]
          x[i, j] - min(row[is.finite(row)], na.rm = TRUE)
        }, numeric(1))
        delta <- delta[is.finite(delta)]
        if (!length(delta)) next
        qs <- stats::quantile(delta, c(alpha / 2, 1 - alpha / 2),
                              names = FALSE, type = 8)
        # A candidate is removed only when its excess phenotype loss is
        # reproducibly above zero. This keeps ties and uncertain alternatives.
        candidate[i, j] <- qs[1] <= 0
        margin_lower[i, j] <- -qs[2]
        margin_upper[i, j] <- -qs[1]
      }
      if (!any(candidate[i, ])) candidate[i, which.min(observed[i, ])] <- TRUE
    }
    blocks$diagnostics[[key]]$candidate <- candidate
    blocks$diagnostics[[key]]$phenotype_margin_lower <- margin_lower
    blocks$diagnostics[[key]]$phenotype_margin_upper <- margin_upper
    sec <- blocks$secondary[[key]]
    combined <- blocks$primary[[key]]
    sec01 <- matrix((pmax(-1, pmin(1, as.numeric(sec))) + 1) / 2,
                    nrow = nrow(sec), ncol = ncol(sec), dimnames = dimnames(sec))
    if (any(abs(sec) > 0, na.rm = TRUE)) {
      for (i in seq_len(nrow(combined))) {
        yes <- which(candidate[i, ])
        no <- which(!candidate[i, ])
        if (length(yes)) combined[i, yes] <- 1 + sec01[i, yes]
        if (length(no)) combined[i, no] <- pmin(blocks$primary[[key]][i, no],
                                                 0.999999)
      }
    }
    blocks$final[[key]] <- combined
  }
  blocks
}

#' Apply max-statistic eligibility to secondary evidence and rebuild ranking
#' @keywords internal
#' @noRd
chorale_gate_secondary <- function(blocks, profiles, spec, cutoff) {
  mods <- names(profiles)
  for (i in seq_along(mods)) for (j in seq_along(mods)) {
    if (j <= i) next
    key <- paste(mods[i], mods[j], sep = "|")
    dg <- blocks$diagnostics[[key]]
    sec <- matrix(0, nrow(dg$loss), ncol(dg$loss), dimnames = dimnames(dg$loss))
    if (length(spec$secondary)) {
      contributions <- lapply(spec$secondary, function(cv) {
        value <- chorale_secondary_block_affinity(
          profiles[[i]], profiles[[j]], cv, dg$orientation)
        value[value <= cutoff] <- 0
        value
      })
      sec <- Reduce(`+`, contributions) / length(contributions)
    }
    blocks$secondary[[key]] <- sec
    combined <- blocks$primary[[key]]
    candidate <- dg$candidate
    if (any(sec > 0, na.rm = TRUE)) {
      for (r in seq_len(nrow(combined))) {
        yes <- which(candidate[r, ])
        no <- which(!candidate[r, ])
        if (length(yes)) combined[r, yes] <- 1 + sec[r, yes]
        if (length(no)) combined[r, no] <- pmin(combined[r, no], 0.999999)
      }
    }
    blocks$final[[key]] <- combined
  }
  blocks
}

#' Resample rows within optional exchangeability blocks
#' @keywords internal
#' @noRd
chorale_bootstrap_rows <- function(design, exchangeability_blocks = NULL) {
  blocks <- intersect(exchangeability_blocks %||% character(), names(design))
  key <- if (!length(blocks)) rep("all", nrow(design)) else {
    do.call(paste, c(lapply(design[blocks], function(x) {
      out <- as.character(x); out[is.na(out)] <- "missing"; out
    }), sep = "\r"))
  }
  unlist(lapply(split(seq_len(nrow(design)), key), function(i) {
    sample(i, length(i), replace = TRUE)
  }), use.names = FALSE)
}

#' Synchronise precomputed pairwise affinities
#' @keywords internal
#' @noRd
chorale_synchronise_affinity <- function(affinity, factor_names) {
  mods <- names(factor_names)
  k <- vapply(factor_names, length, integer(1))
  if (length(mods) < 2L || any(k < 1L)) {
    return(list(assignment = data.frame(), similarity = affinity,
                n_programmes = 0L))
  }
  offset <- cumsum(c(0L, k))
  total <- sum(k)
  scale_by <- max(abs(unlist(affinity)), na.rm = TRUE)
  if (!is.finite(scale_by) || scale_by == 0) scale_by <- 1
  w <- matrix(0, total, total)
  for (i in seq_along(mods)) {
    ri <- offset[i] + seq_len(k[i])
    w[ri, ri] <- diag(k[i])
    for (j in seq_along(mods)) {
      if (j <= i) next
      rj <- offset[j] + seq_len(k[j])
      key <- paste(mods[i], mods[j], sep = "|")
      block <- affinity[[key]] / scale_by
      block[!is.finite(block)] <- 0
      w[ri, rj] <- block
      w[rj, ri] <- t(block)
    }
  }
  d <- min(k)
  e <- eigen(w, symmetric = TRUE)
  vals <- pmax(e$values[seq_len(d)], 0)
  embedding <- e$vectors[, seq_len(d), drop = FALSE] %*%
    diag(sqrt(vals), d, d)
  blocks <- lapply(seq_along(mods), function(i) {
    embedding[offset[i] + seq_len(k[i]), , drop = FALSE]
  })
  leverage <- vapply(blocks, function(x) sum(x^2), numeric(1))
  anchors <- which(k == min(k))
  anchor <- anchors[which.max(leverage[anchors])]
  reference <- blocks[[anchor]]
  rows <- list()
  for (i in seq_along(mods)) {
    stat <- abs(blocks[[i]] %*% t(reference))
    pick <- chorale_assign(stat)
    for (f in seq_len(k[i])) {
      if (is.na(pick[f])) next
      rows[[length(rows) + 1L]] <- data.frame(
        programme = pick[f], modality = mods[i], factor = factor_names[[i]][f],
        factor_index = f, stringsAsFactors = FALSE)
    }
  }
  assignment <- if (length(rows)) do.call(rbind, rows) else data.frame()
  list(assignment = assignment, similarity = affinity,
       n_programmes = length(unique(assignment$programme)))
}

#' Freedman--Lane phenotype null scores
#' @keywords internal
#' @noRd
chorale_phenotype_null_scores <- function(scores, design, spec,
                                          exchangeability_blocks = NULL) {
  design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]
  reduced <- chorale_signature_matrix(design, spec, include = spec$secondary)$x
  out <- scores
  block_key <- rep("all", nrow(design))
  blocks <- intersect(exchangeability_blocks %||% character(), names(design))
  if (length(blocks)) {
    block_key <- do.call(paste, c(lapply(design[blocks], as.character), sep = "\r"))
  }
  for (j in seq_len(ncol(scores))) {
    y <- scores[, j]
    ok <- is.finite(y) & stats::complete.cases(reduced) & !is.na(block_key)
    if (sum(ok) <= ncol(reduced)) next
    fit <- stats::lm.fit(reduced[ok, , drop = FALSE], y[ok])
    res <- fit$residuals
    key <- block_key[ok]
    perm <- seq_along(res)
    for (g in unique(key)) {
      at <- which(key == g)
      perm[at] <- sample(at, length(at), replace = FALSE)
    }
    out[ok, j] <- fit$fitted.values + res[perm]
  }
  out
}

#' Reduced-model residual permutation for one secondary covariate
#' @keywords internal
#' @noRd
chorale_covariate_null_scores <- function(scores, design, spec, covariate,
                                          exchangeability_blocks = NULL) {
  design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]
  keep <- setdiff(spec$covariates, covariate)
  reduced <- chorale_signature_matrix(design, spec, include = keep)$x
  out <- scores
  block_key <- rep("all", nrow(design))
  blocks <- intersect(exchangeability_blocks %||% character(), names(design))
  if (length(blocks)) {
    block_key <- do.call(paste, c(lapply(design[blocks], as.character), sep = "\r"))
  }
  ok <- stats::complete.cases(reduced) & stats::complete.cases(scores) &
    !is.na(block_key)
  if (sum(ok) <= ncol(reduced) || qr(reduced[ok, , drop = FALSE])$rank < ncol(reduced)) {
    return(out)
  }
  fit <- stats::lm.fit(reduced[ok, , drop = FALSE], scores[ok, , drop = FALSE])
  perm <- seq_len(sum(ok))
  key <- block_key[ok]
  for (g in unique(key)) {
    at <- which(key == g)
    perm[at] <- sample(at, length(at), replace = FALSE)
  }
  out[ok, ] <- fit$fitted.values + fit$residuals[perm, , drop = FALSE]
  out
}
