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
