test_that("the vocabulary reports coverage per concept per modality", {
  containers <- sim_containers()
  ids <- rownames(SummarizedExperiment::assay(containers$A))
  sets <- list(wide = ids[1:20], narrow = ids[1:3],
               absent = paste0("no_such_feature_", 1:20))

  cc <- chorale_concepts(containers, sets, min_features = 5)

  expect_s3_class(cc, "chorale_concepts")
  expect_equal(nrow(cc$coverage), length(sets) * 2L)
  expect_setequal(cc$coverage$concept, names(sets))
  expect_equal(cc$coverage$n_features[cc$coverage$concept == "wide"], c(20L, 20L))
  # A concept below the threshold keeps its true count and is not expressed.
  expect_equal(cc$coverage$n_features[cc$coverage$concept == "narrow"], c(3L, 3L))
  expect_false(any(cc$coverage$expressed[cc$coverage$concept == "narrow"]))
  expect_equal(cc$coverage$n_features[cc$coverage$concept == "absent"], c(0L, 0L))
})

test_that("the count expressible in every modality is stated, not assumed", {
  containers <- sim_containers()
  a_ids <- rownames(SummarizedExperiment::assay(containers$A))
  b_ids <- rownames(SummarizedExperiment::assay(containers$B))
  sets <- list(both = intersect(a_ids, b_ids)[1:20],
               neither = paste0("unmeasured_", 1:20))

  cc <- chorale_concepts(containers, sets, min_features = 5)

  expect_equal(cc$n_supplied, 2L)
  expect_equal(cc$n_in_all, 1L)
  expect_equal(cc$n_in_none, 1L)
  expect_true("both" %in% cc$vocabulary)
  expect_true(cc$summary$in_all_modalities[cc$summary$concept == "both"])
})

test_that("a concept reaching one modality is kept and reported as partial", {
  containers <- sim_containers()
  a_ids <- rownames(SummarizedExperiment::assay(containers$A))
  # Give the second modality features of its own, so a concept can reach one
  # modality and not the other.
  b <- SummarizedExperiment::assay(containers$B)
  rownames(b) <- paste0("other_", seq_len(nrow(b)))
  containers$B <- chorale_load(
    b, as.data.frame(SummarizedExperiment::colData(containers$B)))
  sets <- list(a_only = a_ids[1:20])

  cc <- chorale_concepts(containers, sets, min_features = 5)

  expect_equal(cc$n_in_all, 0L)
  expect_equal(cc$n_in_some, 1L)
  expect_equal(cc$vocabulary, "a_only")
  expect_equal(ncol(cc$membership$A), 1L)
  expect_equal(ncol(cc$membership$B), 0L)

  # Requiring two modalities empties the vocabulary rather than pretending.
  strict <- chorale_concepts(containers, sets, min_features = 5,
                             min_modalities = 2L)
  expect_length(strict$vocabulary, 0L)
})

test_that("membership rows are the modality's features in assay order", {
  containers <- sim_containers()
  ids <- rownames(SummarizedExperiment::assay(containers$A))
  cc <- chorale_concepts(containers, list(one = ids[1:20]), min_features = 5)
  expect_identical(rownames(cc$membership$A), ids)
})

test_that("a modality without feature identifiers is refused by name", {
  containers <- sim_containers()
  a <- SummarizedExperiment::assay(containers$A)
  rownames(a) <- NULL
  containers$A <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = a),
    colData = SummarizedExperiment::colData(containers$A))
  expect_error(chorale_concepts(containers, list(one = letters)),
               "Modality 'A'")
})

test_that("a lipid modality reaches the same vocabulary through its classes", {
  containers <- sim_containers(n_features = 40)
  names(containers) <- c("RNA", "LIPID")
  lipids <- c(paste0("SM ", 30:44, ":1;O2"), paste0("Cer ", 30:44, ":1;O2"),
              paste0("PC ", 30:39, ":1"))
  b <- SummarizedExperiment::assay(containers$LIPID)
  b <- b[seq_along(lipids), , drop = FALSE]
  rownames(b) <- lipids
  containers$LIPID <- chorale_load(
    b, as.data.frame(SummarizedExperiment::colData(containers$LIPID)))
  rna_ids <- rownames(SummarizedExperiment::assay(containers$RNA))
  sets <- stats::setNames(
    list(rna_ids[1:20], rna_ids[5:25]),
    c("REACTOME_SPHINGOLIPID_METABOLISM", "REACTOME_CARNITINE_SHUTTLE"))

  cc <- chorale_concepts(containers, sets,
                         feature_space = c(LIPID = "lipid"),
                         min_features = 5)

  expect_equal(unname(cc$feature_space[["LIPID"]]), "lipid")
  expect_true("REACTOME_SPHINGOLIPID_METABOLISM" %in%
                colnames(cc$membership$LIPID))
  expect_true("REACTOME_SPHINGOLIPID_METABOLISM" %in%
                colnames(cc$membership$RNA))
})

test_that("printing states what the vocabulary covers", {
  containers <- sim_containers()
  ids <- rownames(SummarizedExperiment::assay(containers$A))
  cc <- chorale_concepts(containers, list(one = ids[1:20]), min_features = 5)
  expect_output(print(cc), "expressible in every modality: 1")
})
