# Two modalities measured on the same individuals, where a named concept is
# driven by one latent state per individual in both of them. The pairing is what the
# benchmark withholds, so it has to exist before it can be destroyed.
paired_concept_data <- function(n_samples = 60L, n_features = 90L,
                                concept_effect = 1, phenotype_effect = 0.5,
                                seed = 1L) {
  set.seed(seed)
  ids <- sprintf("feature_%04d", seq_len(n_features))
  sets <- list(planted = ids[1:25], quiet = ids[31:55], other = ids[61:85])
  sample_id <- sprintf("s%03d", seq_len(n_samples))
  phenotype <- rep(c("control", "case"), length.out = n_samples)
  sex <- rep(c("F", "F", "M", "M"), length.out = n_samples)

  planted_state <- stats::rnorm(n_samples)
  quiet_state <- stats::rnorm(n_samples)
  build <- function() {
    x <- matrix(stats::rnorm(n_features * n_samples), nrow = n_features,
                dimnames = list(ids, sample_id))
    add <- function(members, state) {
      x[members, ] <<- x[members, ] +
        concept_effect * matrix(rep(state, each = length(members)),
                                nrow = length(members))
    }
    add(sets$planted, planted_state)
    add(sets$quiet, quiet_state)
    x[, phenotype == "case"] <- x[, phenotype == "case"] + phenotype_effect
    x
  }

  list(a = build(), b = build(), sets = sets,
       design = data.frame(sample_id = sample_id, phenotype = phenotype,
                           sex = sex, stringsAsFactors = FALSE))
}
