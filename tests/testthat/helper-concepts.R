# The simulator names each modality's features after that modality, since a
# disjoint collection generally shares none. A shared vocabulary is what these
# tests are about, so the two modalities are given one feature space here.
sim_containers <- function(n_features = 60, n_modalities = 2, seed = 1,
                           shared_features = TRUE) {
  sim <- chorale_simulate(n_modalities = n_modalities, n_features = n_features,
                          seed = seed)
  if (shared_features) {
    ids <- sprintf("feature_%05d", seq_len(n_features))
    sim$modalities <- lapply(sim$modalities, function(m) {
      rownames(m) <- ids
      m
    })
  }
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  names(containers) <- LETTERS[seq_along(containers)]
  containers
}

encode_fixture <- function(n_features = 80, n_modalities = 2, seed = 1) {
  sim <- chorale_simulate(n_modalities = n_modalities, n_features = n_features,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 4, n_per_cell = 3, effect_size = 3,
                          seed = seed)
  ids <- sprintf("feature_%05d", seq_len(n_features))
  sim$modalities <- lapply(sim$modalities, function(m) {
    rownames(m) <- ids
    m
  })
  containers <- Map(chorale_load, sim$modalities, sim$col_data)
  names(containers) <- LETTERS[seq_along(containers)]
  sets <- list(one = ids[1:20], two = ids[21:45], three = ids[40:70])
  list(sim = sim, containers = containers, sets = sets, ids = ids,
       concepts = chorale_concepts(containers, sets, min_features = 5))
}

# A collection in which a continuous covariate is associated with the phenotype
# and drives a concept, while no concept has an effect of its own once that
# covariate is adjusted for. The adjusted phenotype coefficient is null by
# construction, so a calibrated null must reject at the nominal rate.
confounded_collection <- function(n_samples = 60L, n_features = 90L,
                                  slope = 1.5, seed = 1L) {
  set.seed(seed)
  ids <- sprintf("feature_%04d", seq_len(n_features))
  k <- floor((n_features - 15L) / 3L)
  sets <- list(driven = ids[seq_len(k)],
               quiet = ids[k + 5L + seq_len(k)],
               other = ids[2L * k + 10L + seq_len(k)])
  containers <- list()
  for (i in seq_len(2L)) {
    sample_id <- sprintf("m%d_s%03d", i, seq_len(n_samples))
    age <- stats::runif(n_samples, 2, 20)
    # The phenotype follows the covariate, so the two are collinear and a
    # permutation that ignored the covariate would break their relation.
    phenotype <- ifelse(age + stats::rnorm(n_samples, sd = 3) >
                          stats::median(age), "case", "control")
    x <- matrix(stats::rnorm(n_features * n_samples), nrow = n_features,
                dimnames = list(ids, sample_id))
    # The concept responds to the covariate alone.
    x[sets$driven, ] <- x[sets$driven, ] +
      slope * matrix(rep(scale(age)[, 1], each = length(sets$driven)),
                     nrow = length(sets$driven))
    containers[[i]] <- chorale_load(
      x, data.frame(sample_id = sample_id, phenotype = phenotype, age = age,
                    stringsAsFactors = FALSE))
  }
  names(containers) <- c("A", "B")
  list(containers = containers, sets = sets)
}
