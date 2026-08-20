#' Whether the cohorts describe comparable populations
#'
#' A programme recovered across cohorts is a statement about a population, and
#' the population it describes is the one all contributing cohorts represent.
#' Where the cohorts realise the design differently, that population is smaller
#' than any of them, and a claim made outside it rests on cohorts that do not
#' cover it.
#'
#' For each covariate the modalities share, the overlap of their level
#' distributions is reported as the total variation distance and as the common
#' support, the levels every cohort populates. A covariate on which the cohorts
#' barely overlap is one on which their programmes are not comparable, and a
#' level absent from one cohort is a level no claim can be made about.
#'
#' @param fit A `chorale_fit` object.
#' @param min_samples Levels held by fewer than this many samples in a cohort
#'   count as unpopulated there.
#'
#' @returns A data frame with one row per shared covariate, carrying the number
#'   of levels each cohort populates, the levels common to all, the share of
#'   samples those common levels hold, and the largest total variation distance
#'   between any two cohorts.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 80,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2,
#'                    n_ambiguity_boot = 19)
#' chorale_cohort_overlap(fit)
chorale_cohort_overlap <- function(fit, min_samples = 2L) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  designs <- fit$designs
  shared <- chorale_shared_covariates_all(designs, fit$strata_keys)
  if (length(shared) == 0) return(data.frame())

  rows <- lapply(shared, function(cv) {
    levels_by <- lapply(designs, function(d) {
      v <- as.character(d[[cv]])
      tb <- table(v[!is.na(v)])
      names(tb)[tb >= min_samples]
    })
    common <- Reduce(intersect, levels_by)
    covered <- vapply(designs, function(d) {
      v <- as.character(d[[cv]])
      if (length(common) == 0) return(0)
      mean(v %in% common, na.rm = TRUE)
    }, numeric(1))

    # Total variation distance between two cohorts' level distributions: half
    # the sum of absolute differences in proportion, so zero means the cohorts
    # realise the covariate identically and one means they do not overlap.
    all_levels <- sort(unique(unlist(levels_by)))
    props <- lapply(designs, function(d) {
      v <- factor(as.character(d[[cv]]), levels = all_levels)
      p <- table(v) / sum(!is.na(v))
      as.numeric(p)
    })
    tvd <- 0
    ms <- names(designs)
    for (i in seq_along(ms)) {
      for (j in seq_along(ms)) {
        if (j <= i) next
        tvd <- max(tvd, 0.5 * sum(abs(props[[i]] - props[[j]]), na.rm = TRUE))
      }
    }

    data.frame(
      covariate = cv,
      n_levels_each = paste(vapply(levels_by, length, integer(1)), collapse = ", "),
      n_common_levels = length(common),
      common_levels = paste(common, collapse = ", "),
      min_share_covered = round(min(covered), 3),
      max_total_variation = round(tvd, 3),
      comparable = length(common) >= 2 && min(covered) >= 0.5,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' The population a claim from this fit can be made about
#'
#' The common support: the design cells every contributing cohort populates.
#' A programme recovered here describes those animals. Reporting it beside a
#' result is what keeps a claim inside the population the data cover, rather
#' than extending it to cells only one cohort holds.
#'
#' @param fit A `chorale_fit` object.
#' @param min_samples Cells held by fewer than this many samples in a cohort
#'   count as unpopulated there.
#'
#' @returns A list with `cells`, the design cells common to every modality;
#'   `share`, the fraction of each modality's samples those cells hold; and
#'   `restricted`, true where any modality has samples outside the common
#'   support.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 80,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2,
#'                    n_ambiguity_boot = 19)
#' chorale_common_support(fit)$cells
chorale_common_support <- function(fit, min_samples = 2L) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  keys <- intersect(fit$strata_keys, Reduce(intersect, lapply(fit$designs, colnames)))
  if (length(keys) == 0) {
    return(list(cells = character(0), share = numeric(0), restricted = NA))
  }
  cells_by <- lapply(fit$designs, function(d) {
    k <- do.call(paste, c(lapply(keys, function(x) as.character(d[[x]])), sep = "|"))
    k[grepl("NA", k, fixed = TRUE)] <- NA_character_
    tb <- table(k[!is.na(k)])
    names(tb)[tb >= min_samples]
  })
  common <- Reduce(intersect, cells_by)
  share <- vapply(fit$designs, function(d) {
    k <- do.call(paste, c(lapply(keys, function(x) as.character(d[[x]])), sep = "|"))
    if (length(common) == 0) return(0)
    round(mean(k %in% common), 3)
  }, numeric(1))
  list(cells = common, share = share, keys = keys,
       restricted = any(share < 1))
}
