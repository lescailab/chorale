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
