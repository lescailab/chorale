#' Recover latent components within one modality
#'
#' Runs independent component analysis over several random initialisations and
#' keeps the run reaching the highest average non-Gaussianity, measured as mean
#' absolute excess kurtosis of the recovered sources. ICA is non-convex, so a
#' single run is a draw rather than an estimate; the spread across runs is
#' returned alongside the selected fit and is what [chorale_null()] reports as
#' stability.
#'
#' @param x A samples-by-features numeric matrix, centred and scaled.
#' @param n_factors Integer number of components to recover.
#' @param n_init Integer number of random initialisations.
#' @param seed Integer seed. Initialisation `i` uses `seed + i`.
#'
#' @returns A list with `scores` (samples by factors), `loadings` (features by
#'   factors), and `stability`, a data frame of the objective reached by each
#'   initialisation.
#' @keywords internal
#' @noRd
chorale_ica <- function(x, n_factors, n_init = 20L, seed = 1L) {
  rlang::check_installed("fastICA")
  best <- NULL
  best_obj <- -Inf
  obj <- rep(NA_real_, n_init)

  for (i in seq_len(n_init)) {
    set.seed(seed + i)
    fit <- try(
      fastICA::fastICA(x, n.comp = n_factors, method = "C",
                       maxit = 500, tol = 1e-5, verbose = FALSE),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) next
    s <- scale(fit$S)
    o <- mean(abs(apply(s, 2, chorale_excess_kurtosis)), na.rm = TRUE)
    obj[i] <- o
    if (is.finite(o) && o > best_obj) {
      best_obj <- o
      best <- list(scores = s)
    }
  }

  if (is.null(best)) {
    rlang::abort("Independent component analysis failed at every initialisation.")
  }

  # Loadings in feature space: regress each feature on the recovered sources.
  coefs <- stats::coef(stats::lm(x ~ best$scores))
  loadings <- t(coefs[-1, , drop = FALSE])
  colnames(loadings) <- paste0("factor_", seq_len(n_factors))
  rownames(loadings) <- colnames(x)
  colnames(best$scores) <- colnames(loadings)
  rownames(best$scores) <- rownames(x)

  list(
    scores = best$scores,
    loadings = loadings,
    stability = data.frame(init = seq_len(n_init), objective = obj,
                           stringsAsFactors = FALSE)
  )
}

#' Excess kurtosis of a vector
#' @keywords internal
#' @noRd
chorale_excess_kurtosis <- function(v) {
  v <- v[is.finite(v)]
  n <- length(v)
  if (n < 4) return(NA_real_)
  m <- mean(v)
  s <- stats::sd(v)
  if (s == 0) return(NA_real_)
  sum(((v - m) / s)^4) / n - 3
}

#' Select pure features anchoring each factor
#'
#' A pure feature loads on one factor and negligibly on the others, and is what
#' gives a factor a definition independent of the rest of the loading vector.
#' At least two per factor per modality are required for the shared latent
#' structure to be recoverable.
#'
#' Purity is decided first, on the loadings alone, so the identification
#' argument holds whatever the biology turns out to be. Where more features
#' qualify than are needed, the tie is broken on biology: candidates sharing a
#' curated set with the other candidates for the same factor are preferred. Any
#' qualifying set satisfies the condition equally, so choosing the coherent one
#' costs nothing and yields markers that can be read.
#'
#' @param loadings A features-by-factors numeric matrix.
#' @param purity_ratio A feature is pure for a factor when the largest
#'   competing absolute loading is at most this fraction of its own.
#' @param min_markers Minimum qualifying features per factor.
#' @param max_markers Maximum markers retained per factor.
#' @param prior Optional feature-by-set indicator matrix from
#'   [chorale_geneset_matrix()], used only to break ties.
#'
#' @returns A list with `markers`, the features per factor that clear the
#'   purity threshold; `best_candidates`, the most nearly pure features whether
#'   or not they clear it; `purity_margin`, the competing-to-own loading ratio
#'   of the best candidates; and `pure_feature_condition`, recording which
#'   factors reached `min_markers`.
#' @export
#' @examples
#' set.seed(1)
#' l <- matrix(rnorm(60), nrow = 20, dimnames = list(paste0("f", 1:20), NULL))
#' l[1:2, 2:3] <- 0
#' chorale_markers(l)$markers[[1]]
chorale_markers <- function(loadings, purity_ratio = 0.25, min_markers = 2L,
                            max_markers = 20L, prior = NULL) {
  if (!is.matrix(loadings)) rlang::abort("`loadings` must be a matrix.")
  if (is.null(colnames(loadings))) {
    colnames(loadings) <- paste0("factor_", seq_len(ncol(loadings)))
  }
  if (is.null(rownames(loadings))) {
    rownames(loadings) <- paste0("feature_", seq_len(nrow(loadings)))
  }

  a <- abs(loadings)
  markers <- vector("list", ncol(loadings))
  best_candidates <- vector("list", ncol(loadings))
  purity_margin <- rep(NA_real_, ncol(loadings))
  names(markers) <- names(best_candidates) <- names(purity_margin) <-
    colnames(loadings)

  for (j in seq_len(ncol(loadings))) {
    own <- a[, j]
    other <- if (ncol(a) > 1) {
      apply(a[, -j, drop = FALSE], 1, max)
    } else {
      rep(0, nrow(a))
    }
    qualifies <- own > 0 & other <= purity_ratio * own
    candidates <- rownames(a)[qualifies]

    if (length(candidates) > max_markers) {
      score <- own[candidates]
      score <- score / max(score)
      if (!is.null(prior)) {
        shared <- intersect(candidates, rownames(prior))
        if (length(shared) > 1) {
          sub <- prior[shared, , drop = FALSE] > 0
          # Coherence: how many other candidates each one co-occurs with in
          # at least one curated set.
          co <- tcrossprod(sub) > 0
          coherence <- rowSums(co) - 1
          if (max(coherence) > 0) {
            bonus <- stats::setNames(rep(0, length(candidates)), candidates)
            bonus[shared] <- coherence / max(coherence)
            score <- score + bonus
          }
        }
      }
      candidates <- names(sort(score, decreasing = TRUE))[seq_len(max_markers)]
    }
    markers[[j]] <- candidates

    # The purity margin of the best available features, whether or not they
    # clear the threshold. Where the rotation is imperfectly recovered no
    # feature is strictly pure, and reporting the margin distinguishes that
    # case from a factor that is genuinely diffuse.
    margin <- other / own
    margin[!is.finite(margin)] <- Inf
    best <- names(sort(margin))[seq_len(min(max_markers, length(margin)))]
    best_candidates[[j]] <- best
    purity_margin[j] <- stats::median(margin[best[seq_len(min(min_markers, length(best)))]])
  }

  list(
    markers = markers,
    best_candidates = best_candidates,
    purity_margin = purity_margin,
    purity_ratio = purity_ratio,
    pure_feature_condition = vapply(
      markers, function(m) length(m) >= min_markers, logical(1)
    )
  )
}

#' Shrink non-anchor loadings towards curated sets
#'
#' Expresses the loadings of features that are not markers as a combination of
#' curated gene sets, in the manner of PLIER, so each factor carries a pathway
#' definition from estimation rather than from later annotation. Marker
#' loadings are left untouched, since they carry the purity the identification
#' argument rests on and a dense prior would erode it.
#'
#' @param loadings A features-by-factors numeric matrix.
#' @param prior A feature-by-set indicator matrix from
#'   [chorale_geneset_matrix()], with rows matching `loadings`.
#' @param markers A named list of marker features per factor.
#' @param lambda Ridge penalty on the set coefficients.
#'
#' @returns A list with `loadings`, the constrained matrix, and `set_weights`,
#'   a sets-by-factors matrix giving each factor's composition in set space.
#' @export
#' @examples
#' set.seed(1)
#' l <- matrix(rnorm(40), nrow = 10, dimnames = list(paste0("f", 1:10), NULL))
#' p <- matrix(rbinom(30, 1, 0.4), nrow = 10,
#'             dimnames = list(paste0("f", 1:10), paste0("set", 1:3)))
#' dim(chorale_constrain(l, p, rep(list(character()), 4))$set_weights)
chorale_constrain <- function(loadings, prior, markers, lambda = 1) {
  shared <- intersect(rownames(loadings), rownames(prior))
  if (length(shared) < 2 || ncol(prior) < 1) {
    return(list(
      loadings = loadings,
      set_weights = matrix(numeric(0), nrow = 0, ncol = ncol(loadings))
    ))
  }

  c_mat <- prior[shared, , drop = FALSE]
  out <- loadings
  set_weights <- matrix(
    0, nrow = ncol(c_mat), ncol = ncol(loadings),
    dimnames = list(colnames(c_mat), colnames(loadings))
  )

  ctc <- crossprod(c_mat) + lambda * diag(ncol(c_mat))
  inv <- tryCatch(solve(ctc), error = function(e) chorale_ginv(ctc))

  for (j in seq_len(ncol(loadings))) {
    free <- setdiff(shared, markers[[j]])
    if (length(free) < 2) next
    y <- loadings[shared, j]
    y[!(shared %in% free)] <- 0
    u <- inv %*% crossprod(c_mat, y)
    set_weights[, j] <- as.numeric(u)
    fitted <- as.numeric(c_mat %*% u)
    names(fitted) <- shared
    out[free, j] <- fitted[free]
  }

  list(loadings = out, set_weights = set_weights)
}

#' Moore-Penrose inverse
#' @keywords internal
#' @noRd
chorale_ginv <- function(m, tol = sqrt(.Machine$double.eps)) {
  s <- svd(m)
  keep <- s$d > max(tol * s$d[1], 0)
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

#' Mean factor score per design stratum
#' @keywords internal
#' @noRd
chorale_stratum_means <- function(scores, design, strata_keys) {
  keys <- intersect(strata_keys, colnames(design))
  empty <- matrix(numeric(0), nrow = 0, ncol = ncol(scores))
  if (length(keys) == 0) return(empty)
  design <- design[match(rownames(scores), design$sample_id), , drop = FALSE]
  ok <- stats::complete.cases(design[, keys, drop = FALSE])
  if (!any(ok)) return(empty)
  key <- do.call(paste, c(lapply(keys, function(k) as.character(design[[k]])),
                          sep = "|"))
  agg <- rowsum(scores[ok, , drop = FALSE], key[ok], reorder = TRUE)
  counts <- as.numeric(table(key[ok])[rownames(agg)])
  agg / counts
}

#' Match factors across modalities
#'
#' With disjoint samples the modalities share no animal, so factors cannot be
#' matched by correlating scores. What they do share is the design: a factor
#' measuring the same latent state orders the design strata the same way in
#' every modality that measures it. Agreement between stratum profiles is
#' therefore the matching statistic, and it is calibrated by permuting the
#' stratum labels, so the reported matches are those that beat a null built
#' from the same data.
#'
#' Marginal distributional distance is reported alongside as a shape check. It
#' does not drive the matching: independent components are standardised by
#' construction, so their marginals resemble one another whether or not they
#' measure the same thing.
#'
#' @param fits A named list of per-modality fits, each carrying `scores`.
#' @param designs A named list of per-modality design tables.
#' @param strata_keys Design columns defining a stratum.
#' @param n_perm Number of stratum-label permutations calibrating agreement.
#' @param alpha Retain matches with a permutation p-value below this.
#' @param seed Integer seed for the permutations.
#'
#' @returns A data frame, one row per assigned cross-modality factor pair, with
#'   the anchor agreement, its permutation p-value, the marginal
#'   Kolmogorov-Smirnov distance, and the number of shared strata.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 2, n_per_cell = 2, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' fit$matches
chorale_match <- function(fits, designs,
                          strata_keys = c("phenotype", "age_bin", "sex"),
                          n_perm = 200L, alpha = 0.05, seed = 1L) {
  modalities <- names(fits)
  rows <- list()

  for (i in seq_along(modalities)) {
    for (j in seq_along(modalities)) {
      if (j <= i) next
      a <- modalities[i]
      b <- modalities[j]
      sa <- fits[[a]]$scores
      sb <- fits[[b]]$scores
      anchor_a <- chorale_stratum_means(sa, designs[[a]], strata_keys)
      anchor_b <- chorale_stratum_means(sb, designs[[b]], strata_keys)
      common <- intersect(rownames(anchor_a), rownames(anchor_b))

      if (length(common) < 3) {
        # Without shared strata there is no anchor evidence, and marginal
        # matching alone cannot assign factors to one another.
        next
      }

      pa <- anchor_a[common, , drop = FALSE]
      pb <- anchor_b[common, , drop = FALSE]
      agree <- suppressWarnings(abs(stats::cor(pa, pb, method = "spearman")))
      agree[!is.finite(agree)] <- 0

      # One-to-one assignment, so a factor is claimed by at most one partner.
      # solve_LSAP requires at least as many columns as rows, so the wider
      # orientation is solved and the result mapped back.
      if (nrow(agree) > ncol(agree)) {
        flipped <- clue::solve_LSAP(t(agree), maximum = TRUE)
        assignment <- rep(NA_integer_, nrow(agree))
        assignment[as.integer(flipped)] <- seq_len(ncol(agree))
      } else {
        assignment <- as.integer(clue::solve_LSAP(agree, maximum = TRUE))
      }

      # Permutation null: shuffle the stratum labels of one modality and
      # record the best agreement any assignment could reach.
      set.seed(seed)
      null_best <- replicate(n_perm, {
        perm <- sample(seq_len(nrow(pb)))
        g <- suppressWarnings(abs(stats::cor(pa, pb[perm, , drop = FALSE],
                                             method = "spearman")))
        g[!is.finite(g)] <- 0
        max(g)
      })

      for (ca in seq_len(ncol(sa))) {
        cb <- assignment[ca]
        # Where one modality carries more factors than the other, the surplus
        # factors have no partner and are simply unmatched.
        if (is.na(cb)) next
        stat <- agree[ca, cb]
        p <- (1 + sum(null_best >= stat)) / (1 + n_perm)
        va <- sa[, ca]
        vb <- sb[, cb]
        sign_b <- if (stats::cor(pa[, ca], pb[, cb], method = "spearman") < 0) -1 else 1
        d <- min(
          suppressWarnings(stats::ks.test(va, vb)$statistic),
          suppressWarnings(stats::ks.test(va, -vb)$statistic)
        )
        rows[[length(rows) + 1]] <- data.frame(
          modality_a = a, modality_b = b,
          factor_a = colnames(sa)[ca], factor_b = colnames(sb)[cb],
          sign = sign_b,
          anchor_agreement = as.numeric(stat),
          p_value = p,
          ks_distance = as.numeric(d),
          n_shared_strata = length(common),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out$significant <- out$p_value < alpha
  # With few shared strata the permutation test has a floor below which it
  # cannot reach, whatever the data. Recording it keeps an underpowered
  # comparison distinguishable from a genuine absence of agreement.
  out$p_attainable_floor <- 1 / (1 + n_perm)
  out[order(out$p_value, -out$anchor_agreement), , drop = FALSE]
}

#' Add the age band used for anchoring
#' @keywords internal
#' @noRd
chorale_add_age_bin <- function(design) {
  if ("age_bin" %in% colnames(design)) return(design)
  if (!"age_months" %in% colnames(design)) return(design)
  age <- suppressWarnings(as.numeric(as.character(design$age_months)))
  design$age_bin <- cut(age, breaks = c(0, 4, 9, 18, 100),
                        labels = c("2mo", "6mo", "14mo", "aged"))
  design
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Estimate the shared latent state across modalities
#'
#' The estimator of `AGENT_PLAN.md` Section 8.3, in four steps. Components are
#' recovered per modality by independent component analysis over several
#' initialisations. Pure features are then selected on the loadings alone, so
#' the identification argument holds before any prior is applied, with ties
#' broken towards curated-set coherence. The remaining loadings are shrunk
#' towards those sets, so every factor carries a pathway definition at
#' estimation. Factors are finally matched across modalities by distributional
#' agreement, corroborated by the design strata.
#'
#' @param containers A named list of [SummarizedExperiment::SummarizedExperiment]
#'   objects, one per modality, as returned by [chorale_load()].
#' @param n_factors Integer, or one integer per modality, giving the number of
#'   components to recover, from the detectability gate.
#' @param gene_sets A named list of curated gene sets, as returned by
#'   [chorale_genesets()]. Without it, markers are selected on loadings alone
#'   and factors carry no pathway definition.
#' @param feature_map Optional named list, one entry per modality, each a data
#'   frame from [chorale_map()] harmonising that modality's features to Entrez
#'   identifiers.
#' @param n_init Integer number of random initialisations per modality.
#' @param strata_keys Design columns defining an anchoring stratum.
#' @param assay_name Assay to take from each container.
#' @param seed Integer seed.
#'
#' @returns An object of class `chorale_fit`.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 60,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 2, n_per_cell = 2, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' fit
chorale_fit <- function(containers,
                        n_factors,
                        gene_sets = NULL,
                        feature_map = NULL,
                        n_init = 20L,
                        strata_keys = c("phenotype", "age_bin", "sex"),
                        assay_name = NULL,
                        seed = 1L) {
  if (!is.list(containers) || length(containers) < 2) {
    rlang::abort("`containers` must be a list of at least two modalities.")
  }
  if (is.null(names(containers))) {
    names(containers) <- paste0("modality_", seq_along(containers))
  }
  modalities <- names(containers)

  if (length(n_factors) == 1) n_factors <- rep(n_factors, length(containers))
  names(n_factors) <- modalities

  fits <- list()
  designs <- list()

  for (m in modalities) {
    se <- containers[[m]]
    an <- assay_name %||% SummarizedExperiment::assayNames(se)[1]
    mat <- SummarizedExperiment::assay(se, an)
    design <- as.data.frame(SummarizedExperiment::colData(se))
    design <- chorale_add_age_bin(design)

    x <- scale(t(as.matrix(mat)))
    x[!is.finite(x)] <- 0

    fit <- chorale_ica(x, n_factors[[m]], n_init = n_init, seed = seed)

    prior <- NULL
    if (!is.null(gene_sets)) {
      ids <- rownames(mat)
      weights <- rep(1, length(ids))
      if (!is.null(feature_map) && !is.null(feature_map[[m]])) {
        fm <- feature_map[[m]]
        idx <- match(ids, fm$id)
        mapped <- !is.na(idx)
        ids[mapped] <- fm$ENTREZID[idx[mapped]]
        weights[mapped] <- fm$weight[idx[mapped]]
      }
      prior <- chorale_geneset_matrix(ids, gene_sets, weights = weights)
      if (ncol(prior) > 0) rownames(prior) <- rownames(mat)
    }

    mk <- chorale_markers(fit$loadings, prior = prior)
    if (!is.null(prior) && ncol(prior) > 0) {
      constrained <- chorale_constrain(fit$loadings, prior, mk$markers)
      fit$loadings <- constrained$loadings
      fit$set_weights <- constrained$set_weights
    }
    fit$markers <- mk$markers
    fit$best_candidates <- mk$best_candidates
    fit$purity_margin <- mk$purity_margin
    fit$pure_feature_condition <- mk$pure_feature_condition
    fit$prior <- prior

    fits[[m]] <- fit
    designs[[m]] <- design
  }

  matches <- chorale_match(fits, designs, strata_keys = strata_keys,
                           seed = seed)

  structure(
    list(
      modalities = modalities,
      fits = fits,
      designs = designs,
      matches = matches,
      n_shared = if (nrow(matches) > 0) sum(matches$significant) else 0L,
      n_factors = n_factors,
      strata_keys = strata_keys,
      gene_sets = gene_sets,
      seed = seed
    ),
    class = "chorale_fit"
  )
}

#' @export
print.chorale_fit <- function(x, ...) {
  cat("<chorale_fit>\n")
  cat("  modalities:", paste(x$modalities, collapse = ", "), "\n")
  cat("  factors per modality:", paste(x$n_factors, collapse = ", "), "\n")
  for (m in x$modalities) {
    ok <- sum(x$fits[[m]]$pure_feature_condition)
    cat(sprintf("  %-12s %d of %d factors carry pure features\n",
                m, ok, length(x$fits[[m]]$pure_feature_condition)))
  }
  cat("  assigned factor pairs:",
      if (is.data.frame(x$matches)) nrow(x$matches) else 0,
      "of which significant:", x$n_shared, "\n")
  invisible(x)
}
