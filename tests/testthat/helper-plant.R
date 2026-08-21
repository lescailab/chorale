toy_profiles <- function(n_features = 60, seed = 1) {
  sim <- chorale_simulate(n_modalities = 2, n_features = n_features,
                          n_shared_factors = 2, n_private_factors = 1,
                          n_strains = 3, n_per_cell = 3, seed = seed)
  out <- lapply(seq_len(2), function(i) {
    chorale_data_profile(sim$modalities[[i]], sim$col_data[[i]],
                         covariates = c("phenotype", "sex"),
                         layer = paste0("layer_", i))
  })
  names(out) <- c("A", "B")
  out
}

toy_membership <- function(profiles, sets) {
  lapply(profiles, function(p) {
    ids <- p$feature$feature_id
    m <- vapply(sets, function(s) as.numeric(seq_along(ids) %in% s),
                numeric(length(ids)))
    rownames(m) <- ids
    m
  })
}
