## Derive the committed test fixtures from <project root>/test_data/.
##
## The fixtures are the smallest subset of the validated production layers
## that exercises every code path: three modalities on disjoint samples, a
## populated anchoring stratum in each layer that resolves covariates, real
## identifiers, and the real missingness pattern of the proteome.
##
## Run from the package root after `analysis/build_test_data.py` has written
## <project root>/test_data/. The output is committed, so continuous
## integration never reads outside the repository.
##
##   Rscript data-raw/fixtures.R

test_data_dir <- normalizePath(file.path("..", "test_data"), mustWork = TRUE)
fixture_dir <- file.path("tests", "testthat", "fixtures")
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)

n_features <- c(RNA = 150L, PROT = 150L, METAB = 150L)
n_per_stratum <- 3L
stratum_keys <- c("phenotype", "age_bin", "sex")

for (layer in c("RNA", "PROT", "METAB")) {
  mat <- nanoparquet::read_parquet(
    file.path(test_data_dir, paste0(layer, "_matrix.parquet"))
  )
  design <- utils::read.delim(
    file.path(test_data_dir, paste0(layer, "_design.tsv")),
    stringsAsFactors = FALSE
  )

  feature_id <- rownames(mat)
  if (is.null(feature_id)) feature_id <- as.character(seq_len(nrow(mat)))

  # Sample selection: up to n_per_stratum per anchoring stratum, so every
  # stratum retained is still able to support stratified permutation. Layers
  # whose covariates do not resolve contribute a plain random subset and
  # participate in marginal matching only.
  resolvable <- stats::complete.cases(design[, stratum_keys, drop = FALSE])
  if (any(resolvable)) {
    d <- design[resolvable, , drop = FALSE]
    key <- interaction(d[, stratum_keys], drop = TRUE)
    picked <- unlist(lapply(split(seq_len(nrow(d)), key), function(idx) {
      utils::head(sample(idx), n_per_stratum)
    }), use.names = FALSE)
    keep_ids <- d$sample_id[sort(picked)]
  } else {
    keep_ids <- sort(sample(design$sample_id, min(24L, nrow(design))))
  }

  keep_ids <- intersect(keep_ids, colnames(mat))
  sub <- mat[, keep_ids, drop = FALSE]

  # Feature selection: most variable by median absolute deviation, computed
  # on the retained samples so the fixture is internally consistent.
  m <- as.matrix(sub)
  mad_stat <- apply(m, 1, function(r) stats::median(abs(r - stats::median(r, na.rm = TRUE)), na.rm = TRUE))
  ord <- order(mad_stat, decreasing = TRUE)
  take <- utils::head(ord, min(n_features[[layer]], nrow(m)))
  m <- m[sort(take), , drop = FALSE]
  keep_feature <- feature_id[sort(take)]

  out <- as.data.frame(m, check.names = FALSE)
  out <- cbind(feature_id = keep_feature, out)

  nanoparquet::write_parquet(
    out, file.path(fixture_dir, paste0(layer, "_matrix.parquet"))
  )
  utils::write.table(
    design[design$sample_id %in% keep_ids, , drop = FALSE],
    file.path(fixture_dir, paste0(layer, "_design.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  message(sprintf(
    "%-6s %4d features x %3d samples", layer, nrow(m), ncol(m)
  ))
}

message("fixtures written to ", fixture_dir)
