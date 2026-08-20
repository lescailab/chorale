# Reproducible benchmark for the confirmed repeated-model bottleneck.
# Run from the package root after activating the repository environment.

devtools::load_all(quiet = TRUE)
set.seed(1)
n <- 2000L
k <- 20L
design <- data.frame(
  sample_id = paste0("sample_", seq_len(n)),
  phenotype = rep(c("control", "case"), length.out = n),
  age = stats::rnorm(n),
  sex = rep(c("F", "M"), length.out = n),
  batch = rep(paste0("batch_", 1:4), length.out = n),
  site = rep(paste0("site_", 1:5), length.out = n)
)
scores <- matrix(stats::rnorm(n * k), n, k,
                 dimnames = list(design$sample_id, paste0("factor_", seq_len(k))))
spec <- chorale_resolve_signature(list(a = design, b = design))

legacy <- function(scores, design, spec) {
  mm <- chorale_signature_matrix(design, spec)
  terms <- colnames(mm$x)[-1L]
  effects <- se <- matrix(NA_real_, ncol(scores), length(terms))
  for (j in seq_len(ncol(scores))) {
    ok <- is.finite(scores[, j]) & stats::complete.cases(mm$x)
    x <- mm$x[ok, , drop = FALSE]
    fit <- stats::lm.fit(x, scores[ok, j])
    sigma2 <- sum(fit$residuals^2) / fit$df.residual
    vc <- sigma2 * chol2inv(qr.R(fit$qr))
    effects[j, ] <- fit$coefficients[-1L]
    se[j, ] <- sqrt(pmax(diag(vc[-1L, -1L, drop = FALSE]), 0))
  }
  list(effects = effects, se = se)
}

elapsed <- function(fun, repeats = 30L) {
  unname(system.time(for (i in seq_len(repeats)) fun(scores, design, spec))["elapsed"])
}
old <- elapsed(legacy)
new <- elapsed(chorale_adjusted_profile)
speedup <- old / new
result <- data.frame(samples = n, factors = k, repeats = 30L,
                     legacy_seconds = old, vectorised_seconds = new,
                     speedup = speedup)
print(result, row.names = FALSE)
if (!is.finite(speedup) || speedup < 2) {
  stop("Vectorised adjusted-profile calculation did not reach 2x speedup.")
}
