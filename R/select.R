#' Choose how many free dimensions the data support
#'
#' A factor count is usually asserted, or read off a variance curve. Neither
#' says whether the factors at that count are the same factors when the
#' estimator is started again, and independent component analysis is non-convex,
#' so a component recovered at one initialisation and not at another is a draw
#' rather than an estimate.
#'
#' The count is therefore chosen by reproducibility. At each candidate count the
#' factorisation is run repeatedly from different initialisations, and again on
#' subsamples of the samples; every run's components are matched one to one with
#' a reference run's and the matched correlation recorded. The selected count is
#' the largest at which the weakest component still clears a declared threshold,
#' so adding a component that only appears sometimes cannot raise the count.
#'
#' Two kinds of instability are separated because they mean different things.
#' A component that moves between initialisations is a property of the
#' optimiser. A component that moves between subsamples is a property of the
#' data: it rests on which samples were drawn, so it would not survive a
#' replication of the study.
#'
#' Nothing here reads the design. Selection uses the assay alone, so the count
#' cannot have been chosen to suit the contrast the fit is later asked about.
#'
#' @param x A samples-by-features numeric matrix, centred and scaled.
#' @param max_factors Largest count considered. Parallel analysis by
#'   [chorale_n_factors()] is the natural source of this upper bound.
#' @param counts Candidate counts. Defaults to every count from one to
#'   `max_factors`.
#' @param n_init Initialisations compared at each candidate count.
#' @param n_subsample Subsamples compared at each candidate count. Zero skips
#'   the subsample criterion.
#' @param subsample_fraction Share of the samples each subsample draws.
#' @param threshold Matched correlation the weakest component must reach.
#' @param seed Integer seed.
#'
#' @returns A data frame with one row per candidate count, carrying the weakest
#'   and typical matched correlation across initialisations and across
#'   subsamples, and whether the count is admissible. Attribute `selected` holds
#'   the chosen count and `threshold` the value it was held to.
#' @export
#' @examples
#' # Two heavy-tailed sources in noise: beyond them there is nothing to recover
#' # reproducibly, which is where the selection stops.
#' set.seed(1)
#' sources <- matrix(stats::rt(120, df = 3), nrow = 60, ncol = 2)
#' x <- scale(sources %*% matrix(stats::rnorm(80), nrow = 2) +
#'              matrix(stats::rnorm(60 * 40, sd = 2), nrow = 60))
#' sel <- chorale_select_factors(x, max_factors = 3, n_init = 3,
#'                               n_subsample = 3)
#' attr(sel, "selected")
chorale_select_factors <- function(x, max_factors = 10L, counts = NULL,
                                   n_init = 5L, n_subsample = 5L,
                                   subsample_fraction = 0.8,
                                   threshold = 0.75, seed = 1L) {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0
  if (nrow(x) < 4 || ncol(x) < 2) {
    rlang::abort("`x` needs at least four samples and two features.")
  }
  if (threshold <= 0 || threshold > 1) {
    rlang::abort("`threshold` must lie in (0, 1].")
  }
  if (subsample_fraction <= 0 || subsample_fraction >= 1) {
    rlang::abort("`subsample_fraction` must lie strictly between 0 and 1.")
  }
  # A centred matrix has at most `n - 1` independent directions, and at most one
  # per feature, so no candidate count above either is meaningful whatever the
  # caller asked for.
  ceiling_k <- min(as.integer(max_factors), nrow(x) - 1L, ncol(x))
  counts <- counts %||% seq_len(max(ceiling_k, 1L))
  counts <- sort(unique(as.integer(counts)))
  counts <- counts[counts >= 1L & counts <= ceiling_k]
  if (length(counts) == 0) {
    out <- chorale_empty_selection()
    attr(out, "selected") <- 0L
    attr(out, "threshold") <- threshold
    return(out)
  }

  # Four is the smallest subsample a factorisation can be attempted on, and is
  # the floor this function validates its own input against. A subsample below
  # it would fail rather than report a low agreement, and a failure is read here
  # as an unrecoverable count.
  n_keep <- max(4L, floor(nrow(x) * subsample_fraction))
  rows <- lapply(counts, function(k) {
    reference <- try(chorale_single_ica(x, k, seed), silent = TRUE)
    if (inherits(reference, "try-error")) {
      return(data.frame(
        n_factors = k, init_weakest = NA_real_, init_typical = NA_real_,
        subsample_weakest = NA_real_, subsample_typical = NA_real_,
        weakest = NA_real_, admissible = FALSE, stringsAsFactors = FALSE))
    }

    # The strides keep three families of runs on seeds that cannot collide: the
    # reference at `seed`, the initialisations at multiples of 1000, and the
    # subsamples at multiples of 5000 offset by the candidate count. A collision
    # would make two runs identical and report agreement that is really the same
    # run compared with itself.
    init <- chorale_component_agreement(reference, lapply(
      seq_len(n_init), function(i) {
        try(chorale_single_ica(x, k, seed + 1000L * i), silent = TRUE)
      }))

    sub <- if (n_subsample > 0) {
      chorale_component_agreement(reference, lapply(
        seq_len(n_subsample), function(i) {
          set.seed(seed + 5000L * i + k)
          keep <- sample(seq_len(nrow(x)), n_keep)
          try(chorale_single_ica(x[keep, , drop = FALSE], k, seed + i),
              silent = TRUE)
        }))
    } else {
      rep(NA_real_, k)
    }

    weakest <- chorale_finite_min(c(init, sub))
    data.frame(
      n_factors = k,
      init_weakest = round(chorale_finite_min(init), 4),
      init_typical = round(chorale_finite_median(init), 4),
      subsample_weakest = round(chorale_finite_min(sub), 4),
      subsample_typical = round(chorale_finite_median(sub), 4),
      weakest = round(weakest, 4),
      admissible = isTRUE(weakest >= threshold),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  # The largest admissible count, not the largest count reached: a count is
  # admissible only if every count below it is too, since a component that
  # cannot be recovered at a smaller model will not be recovered at a larger
  # one, and an isolated admissible count above an inadmissible one is noise.
  admissible <- out$admissible
  # cumprod over the admissibility flags is zero from the first inadmissible
  # count onwards, so the selected count is the end of the leading run rather
  # than the largest count that happens to be admissible.
  run <- cumprod(as.integer(admissible))
  selected <- if (any(run == 1L)) max(out$n_factors[run == 1L]) else 0L

  attr(out, "selected") <- as.integer(selected)
  attr(out, "threshold") <- threshold
  rownames(out) <- NULL
  out
}

#' Summaries that report absence rather than an infinite bound
#' @keywords internal
#' @noRd
chorale_finite_min <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) NA_real_ else min(v)
}

#' @keywords internal
#' @noRd
chorale_finite_median <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) NA_real_ else stats::median(v)
}

#' An empty selection table
#' @keywords internal
#' @noRd
chorale_empty_selection <- function() {
  data.frame(n_factors = integer(), init_weakest = numeric(),
             init_typical = numeric(), subsample_weakest = numeric(),
             subsample_typical = numeric(), weakest = numeric(),
             admissible = logical(), stringsAsFactors = FALSE)
}

#' One factorisation from one initialisation
#'
#' `chorale_ica()` runs many initialisations and selects among them, which is
#' the wrong unit here: selection compares single runs, so it needs one run at a
#' time.
#'
#' @keywords internal
#' @noRd
chorale_single_ica <- function(x, k, seed) {
  fit <- chorale_ica(x, k, n_init = 1L, seed = seed, consensus = FALSE)
  fit$loadings
}

#' How well each reference component reappears in a set of runs
#'
#' Components are compared in feature space rather than in sample space, because
#' a run on a subsample has no scores for the samples it left out. Sign and
#' order are not identified, so runs are matched one to one on absolute
#' correlation before anything is read off.
#'
#' @keywords internal
#' @noRd
chorale_component_agreement <- function(reference, runs) {
  k <- ncol(reference)
  values <- matrix(NA_real_, nrow = length(runs), ncol = k)
  for (i in seq_along(runs)) {
    s <- runs[[i]]
    if (inherits(s, "try-error") || is.null(s) ||
        !identical(dim(s), dim(reference))) next
    a <- abs(suppressWarnings(stats::cor(reference, s)))
    a[!is.finite(a)] <- 0
    assign <- as.integer(clue::solve_LSAP(a, maximum = TRUE))
    values[i, ] <- a[cbind(seq_len(k), assign)]
  }
  apply(values, 2, function(v) if (all(is.na(v))) NA_real_ else
    mean(v, na.rm = TRUE))
}
