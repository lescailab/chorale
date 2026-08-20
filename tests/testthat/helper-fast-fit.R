# Production uses 999 ambiguity resamples so both tails of a 95% interval are
# represented. Unit tests exercise the same code path with fewer resamples;
# calibration and determinism have dedicated tests.
chorale_fit_production <- chorale_fit
chorale_fit <- function(...) {
  args <- list(...)
  if (is.null(args$n_ambiguity_boot)) args$n_ambiguity_boot <- 19L
  do.call(chorale_fit_production, args)
}
