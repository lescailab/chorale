#' Plant named pathways as the shared programmes of a simulation
#'
#' Simulated data in which the shared programmes are curated pathways rather
#' than anonymous factors. Each planted programme occupies the features that
#' belong to one set in every modality that can express it: the genes of the
#' set in a transcriptome or proteome, and the lipid classes the same set acts
#' on in a lipidome. Recovery can then be scored on biology, which no real
#' cohort supports, since no real cohort records which features belong to a
#' programme.
#'
#' A concept is planted and recovered by name, so the planting and scoring
#' vocabularies are the same one. The question is not whether a concept\'s
#' composition can be rediscovered but whether the design shows it separating
#' cases from controls in every modality that expresses it, and whether the
#' concepts that were not planted stay quiet.
#'
#' Two dials govern how far a planted signal departs from a clean set
#' indicator. `member_fraction` is the share of a set's features that carry the
#' programme, and `leak_fraction` is the share of the programme's loading mass
#' placed on features outside the set. A planted truth is therefore recorded as
#' the features that actually received loading, not as the nominal pathway, and
#' identifier mapping losses are absorbed the same way.
#'
#' @param profiles A named list of [chorale_data_profile()] objects, one per
#'   modality, carrying the real feature identifiers a pathway is defined on.
#' @param membership A named list matching `profiles`, each a feature-by-set
#'   matrix placing that modality's features in the planting vocabulary, from
#'   [chorale_geneset_matrix()] for genes or [chorale_metabolite_matrix()] for
#'   lipids. Row names must be the profile's feature identifiers.
#' @param plant_sets The planting collection, a named list of member
#'   identifiers, as returned by [chorale_genesets()].
#' @param score_sets The collection recovery is scored in, which is the planting
#'   collection itself. The realised overlap between planted sets travels with
#'   the result.
#' @param n_programmes Number of pathways to plant.
#' @param n_private_factors Modality-private factors carrying no pathway and no
#'   cross-modality partner.
#' @param min_modalities Modalities a set must reach to be plantable.
#' @param min_features Features a set must match in a modality to reach it.
#' @param member_fraction Share of a set's features in a modality that carry the
#'   programme.
#' @param leak_fraction Share of the programme's loading mass placed on features
#'   outside the set.
#' @param n_markers Members per programme made pure, so the pure features of a
#'   recovered factor can be checked against the pathway.
#' @param background_sd Standard deviation of the loadings of features carrying
#'   no programme.
#' @param seed Integer seed.
#' @param ... Passed to [chorale_simulate()], for example `effect_size` or
#'   `noise_sd`.
#'
#' @returns A list with `sim`, the [chorale_simulate()] output; `plant`, one row
#'   per programme and modality recording the set, its size, the features that
#'   received loading and the markers among them; `sets`, the planted sets in
#'   their planting vocabulary; and `concepts`, their names.
#' @export
#' @examplesIf FALSE
#' plant <- chorale_plant(profiles, membership, kegg, reactome)
chorale_plant <- function(profiles, membership, plant_sets, score_sets,
                          n_programmes = 3L,
                          n_private_factors = 2L,
                          min_modalities = 2L,
                          min_features = 10L,
                          member_fraction = 0.6,
                          leak_fraction = 0.2,
                          n_markers = 5L,
                          background_sd = 0.2,
                          seed = 1L,
                          ...) {
  if (!is.list(profiles) || length(profiles) < 2) {
    rlang::abort("`profiles` must be a list of at least two modalities.")
  }
  if (is.null(names(profiles))) {
    names(profiles) <- paste0("modality_", seq_along(profiles))
  }
  if (!identical(sort(names(membership)), sort(names(profiles)))) {
    rlang::abort("`membership` must have one entry per modality of `profiles`.")
  }
  if (member_fraction <= 0 || member_fraction > 1) {
    rlang::abort("`member_fraction` must lie in (0, 1].")
  }
  if (leak_fraction < 0 || leak_fraction >= 1) {
    rlang::abort("`leak_fraction` must lie in [0, 1).")
  }
  membership <- membership[names(profiles)]

  for (m in names(profiles)) {
    rn <- rownames(membership[[m]])
    if (is.null(rn) || !identical(rn, profiles[[m]]$feature$feature_id)) {
      rlang::abort(paste0(
        "Row names of `membership$", m,
        "` must be the feature identifiers of the matching profile, in order."
      ))
    }
  }

  candidates <- chorale_plantable_sets(membership, min_features, min_modalities)
  if (nrow(candidates) == 0) {
    rlang::abort(paste0(
      "No set reaches ", min_modalities, " modalities with at least ",
      min_features, " features."
    ))
  }
  overlap <- chorale_set_overlap(plant_sets[candidates$set], score_sets)
  candidates$max_jaccard <- overlap$max_jaccard[match(candidates$set,
                                                      overlap$set)]

  # The widest programmes first, since a programme reaching every modality is
  # what the integration is being asked to recover, and among those the ones
  # matching most features, which carry the clearest planted truth.
  candidates <- candidates[order(-candidates$n_modalities,
                                 -candidates$min_features), , drop = FALSE]
  chosen <- utils::head(candidates, as.integer(n_programmes))
  if (nrow(chosen) < n_programmes) {
    rlang::warn(paste0("Only ", nrow(chosen), " sets are plantable."))
  }
  n_programmes <- nrow(chosen)

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  })
  set.seed(seed)

  built <- lapply(names(profiles), function(m) {
    chorale_plant_loadings(
      profile = profiles[[m]], membership = membership[[m]],
      sets = chosen$set, n_private_factors = n_private_factors,
      member_fraction = member_fraction, leak_fraction = leak_fraction,
      n_markers = n_markers, background_sd = background_sd
    )
  })
  names(built) <- names(profiles)

  # Two programmes with the same design response are one programme as far as
  # the design channel can tell, so the planted signatures are spread over the
  # terms the designs share rather than each taking one term in turn.
  dots <- list(...)
  if (is.null(dots$signature)) {
    dots$signature <- chorale_spread_signature(n_programmes,
                                               chorale_signature_terms(unname(profiles)))
  }
  if (is.null(dots$effect_size)) dots$effect_size <- 3

  sim <- do.call(chorale_simulate, c(list(
    n_modalities = length(profiles),
    n_features = NULL,
    n_shared_factors = n_programmes,
    n_private_factors = n_private_factors,
    profile = unname(profiles),
    loadings = lapply(built, function(b) b$loadings),
    seed = seed
  ), dots))
  # The modalities keep the names of the profiles they were built from, so a
  # planted pathway can be traced to the layer it was planted in.
  names(sim$modalities) <- names(profiles)
  names(sim$col_data) <- names(profiles)
  names(sim$truth$loadings) <- names(profiles)
  names(sim$truth$scores) <- names(profiles)
  names(sim$truth$markers) <- names(profiles)
  for (i in seq_along(sim$col_data)) {
    sim$col_data[[i]]$modality <- names(profiles)[i]
  }

  rows <- list()
  for (m in names(profiles)) {
    b <- built[[m]]
    for (k in seq_len(n_programmes)) {
      rows[[length(rows) + 1L]] <- data.frame(
        programme = paste0("shared_", k),
        set = chosen$set[k],
        modality = m,
        n_members = if (chosen$set[k] %in% colnames(membership[[m]])) {
          sum(membership[[m]][, chosen$set[k]] > 0)
        } else {
          0L
        },
        n_planted = length(b$planted[[k]]),
        n_leaked = length(b$leaked[[k]]),
        max_jaccard = chosen$max_jaccard[k],
        stringsAsFactors = FALSE
      )
    }
  }
  plant <- do.call(rbind, rows)
  plant$planted <- unlist(lapply(built, function(b) b$planted),
                          recursive = FALSE, use.names = FALSE)
  plant$markers <- unlist(lapply(built, function(b) b$markers),
                          recursive = FALSE, use.names = FALSE)

  list(sim = sim, plant = plant,
       sets = plant_sets[chosen$set],
       concepts = chosen$set,
       programmes = stats::setNames(chosen$set, paste0("shared_", seq_len(n_programmes))))
}

#' Sets reaching enough modalities to be planted
#' @keywords internal
#' @noRd
chorale_plantable_sets <- function(membership, min_features, min_modalities) {
  counts <- lapply(membership, function(mm) {
    if (is.null(mm) || ncol(mm) == 0) return(numeric(0))
    cs <- colSums(mm > 0)
    cs[cs >= min_features]
  })
  all_sets <- unique(unlist(lapply(counts, names), use.names = FALSE))
  if (length(all_sets) == 0) {
    return(data.frame(set = character(), n_modalities = integer(),
                      min_features = numeric(), stringsAsFactors = FALSE))
  }
  reached <- vapply(all_sets, function(s) {
    sum(vapply(counts, function(cs) s %in% names(cs), logical(1)))
  }, integer(1))
  smallest <- vapply(all_sets, function(s) {
    v <- vapply(counts, function(cs) if (s %in% names(cs)) cs[[s]] else NA_real_,
                numeric(1))
    min(v, na.rm = TRUE)
  }, numeric(1))
  out <- data.frame(set = all_sets, n_modalities = reached,
                    min_features = smallest, stringsAsFactors = FALSE)
  out[out$n_modalities >= min_modalities, , drop = FALSE]
}

#' Largest overlap of each planting set with any scoring set
#'
#' The Jaccard index of two sets is the size of their intersection over the size
#' of their union. Planting from a set that nearly coincides with a scoring set
#' would let the pathway channel recover it without doing any work, so the
#' largest such overlap is what admits or refuses a planting set.
#'
#' @param plant_sets A named list of member identifiers to be planted.
#' @param score_sets A named list of member identifiers the recovery is scored
#'   in.
#'
#' @returns A data frame with one row per planting set, carrying the scoring set
#'   it overlaps most and the Jaccard index of that overlap.
#' @export
#' @examples
#' a <- list(one = c("1", "2", "3"), two = c("7", "8"))
#' b <- list(x = c("1", "2", "3", "4"), y = c("9"))
#' chorale_set_overlap(a, b)
chorale_set_overlap <- function(plant_sets, score_sets) {
  if (length(plant_sets) == 0) {
    return(data.frame(set = character(), nearest = character(),
                      max_jaccard = numeric(), stringsAsFactors = FALSE))
  }
  score_sets <- lapply(score_sets, function(s) unique(as.character(s)))
  rows <- lapply(names(plant_sets), function(s) {
    a <- unique(as.character(plant_sets[[s]]))
    if (length(a) == 0 || length(score_sets) == 0) {
      return(data.frame(set = s, nearest = NA_character_, max_jaccard = 0,
                        stringsAsFactors = FALSE))
    }
    j <- vapply(score_sets, function(b) {
      inter <- length(intersect(a, b))
      if (inter == 0) return(0)
      inter / length(union(a, b))
    }, numeric(1))
    best <- which.max(j)
    data.frame(set = s, nearest = names(score_sets)[best],
               max_jaccard = unname(j[best]), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Loadings placing each planted set on the features that belong to it
#' @keywords internal
#' @noRd
chorale_plant_loadings <- function(profile, membership, sets,
                                   n_private_factors, member_fraction,
                                   leak_fraction, n_markers, background_sd) {
  ids <- profile$feature$feature_id
  p <- length(ids)
  n_shared <- length(sets)
  k_total <- n_shared + n_private_factors

  l <- matrix(stats::rnorm(p * k_total, sd = background_sd), nrow = p,
              dimnames = list(ids, NULL))

  planted <- vector("list", n_shared)
  leaked <- vector("list", n_shared)
  markers <- vector("list", n_shared)
  names(planted) <- names(leaked) <- names(markers) <-
    paste0("shared_", seq_len(n_shared))

  for (k in seq_len(n_shared)) {
    s <- sets[k]
    members <- if (s %in% colnames(membership)) {
      which(membership[, s] > 0)
    } else {
      integer(0)
    }
    if (length(members) == 0) {
      planted[[k]] <- character(0)
      leaked[[k]] <- character(0)
      markers[[k]] <- character(0)
      next
    }
    n_carry <- max(1L, round(member_fraction * length(members)))
    carriers <- sample(members, n_carry)
    l[carriers, k] <- l[carriers, k] + stats::runif(n_carry, 0.5, 1.5)

    # Loading mass outside the set, so a recovered factor is not a clean set
    # indicator and the pathway channel has to find the set among features that
    # do not belong to it.
    n_leak <- round(leak_fraction / (1 - leak_fraction) * n_carry)
    outside <- setdiff(seq_len(p), members)
    n_leak <- min(n_leak, length(outside))
    leak_idx <- if (n_leak > 0) sample(outside, n_leak) else integer(0)
    if (n_leak > 0) {
      l[leak_idx, k] <- l[leak_idx, k] + stats::runif(n_leak, 0.5, 1.5)
    }

    # Pure features are drawn from the set's own members, so the markers a
    # recovered factor reports can be checked against the pathway that was
    # planted.
    take <- utils::head(carriers[order(-abs(l[carriers, k]))],
                        min(as.integer(n_markers), n_carry))
    l[take, -k] <- 0

    planted[[k]] <- ids[carriers]
    leaked[[k]] <- ids[leak_idx]
    markers[[k]] <- ids[take]
  }

  # A private factor occupies a block of features of its own, so it is
  # recoverable without corresponding to any pathway or to any other modality.
  if (n_private_factors > 0) {
    block <- max(10L, round(0.02 * p))
    for (j in seq_len(n_private_factors)) {
      idx <- sample(seq_len(p), min(block, p))
      l[idx, n_shared + j] <- l[idx, n_shared + j] + stats::runif(length(idx), 0.5, 1.5)
    }
  }

  list(loadings = l, planted = planted, leaked = leaked, markers = markers)
}


#' Probability that a member of the truth set outranks a non-member
#' @keywords internal
#' @noRd
chorale_rank_auc <- function(score, is_truth) {
  ok <- is.finite(score)
  score <- score[ok]
  is_truth <- is_truth[ok]
  n1 <- sum(is_truth)
  n0 <- sum(!is_truth)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[is_truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Design signatures spread over the terms the designs share
#'
#' Programmes are told apart on the design by the direction of their response,
#' so planting several of them needs directions that differ. Where the shared
#' terms are few, giving each programme a term of its own would repeat a
#' direction; spreading the programmes evenly over the space the terms span
#' keeps every pair distinguishable however few terms there are.
#'
#' @keywords internal
#' @noRd
chorale_spread_signature <- function(n_programmes, terms) {
  n <- length(terms)
  sig <- matrix(0, nrow = n_programmes, ncol = n,
                dimnames = list(NULL, terms))
  if (n == 1) {
    sig[, 1] <- 1
    return(sig)
  }
  angles <- seq(0, pi, length.out = n_programmes + 1L)[seq_len(n_programmes)]
  for (k in seq_len(n_programmes)) {
    partner <- ((k - 1L) %% (n - 1L)) + 2L
    sig[k, 1] <- cos(angles[k])
    sig[k, partner] <- sin(angles[k])
  }
  sig
}
