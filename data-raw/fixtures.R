## Generate privacy-safe, reproducible package fixtures.
##
## Run from the package root:
##   Rscript data-raw/fixtures.R

set.seed(1042)

layers <- c("RNA", "PROT", "METAB")
destinations <- c(file.path("tests", "testthat", "fixtures"),
                  file.path("inst", "fixtures"))
invisible(lapply(destinations, dir.create, recursive = TRUE,
                 showWarnings = FALSE))

for (layer_index in seq_along(layers)) {
  layer <- layers[layer_index]
  cells <- expand.grid(
    phenotype = c("control", "case"),
    sex = c("F", "M"),
    age = c(1, 2, 3),
    replicate = seq_len(2),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)
  n <- nrow(cells)
  sample_id <- sprintf("%s_sample_%03d", tolower(layer), seq_len(n))
  design <- data.frame(
    sample_id = sample_id,
    modality = layer,
    phenotype = cells$phenotype,
    sex = cells$sex,
    age = cells$age,
    batch = paste0("batch_", 1L + (seq_len(n) %% 2L)),
    stringsAsFactors = FALSE)

  phenotype <- ifelse(design$phenotype == "case", 1, -1)
  sex <- ifelse(design$sex == "M", 1, -1)
  age <- as.numeric(scale(design$age))
  scores <- cbind(
    phenotype + stats::rnorm(n, sd = 0.25),
    0.7 * phenotype + 0.8 * age + stats::rnorm(n, sd = 0.3),
    sex + stats::rnorm(n, sd = 0.35))
  loadings <- matrix(stats::rt(90 * 3, df = 5), nrow = 90)
  assay <- loadings %*% t(scores) + matrix(stats::rnorm(90 * n, sd = 0.7), 90)
  if (layer == "PROT") {
    assay[sample(length(assay), floor(0.03 * length(assay)))] <- NA
  }
  rownames(assay) <- sprintf("%s_feature_%04d", tolower(layer),
                             seq_len(nrow(assay)))
  colnames(assay) <- sample_id

  stored <- data.frame(feature_id = rownames(assay), assay,
                       check.names = FALSE)
  for (destination in destinations) {
    utils::write.table(
      stored, file.path(destination, paste0(layer, "_matrix.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    utils::write.table(
      design, file.path(destination, paste0(layer, "_design.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

message("Synthetic fixtures written to tests/testthat/fixtures and inst/fixtures")
