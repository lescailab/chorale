#' Terms a design profile can carry in every modality
#'
#' A profile entry is only comparable across modalities if it means the same
#' thing in each, so the terms are fixed once for the whole collection rather
#' than per pair. A covariate that is numeric everywhere contributes one term.
#' A categorical covariate contributes one term per level beyond a reference,
#' using only the levels that every design realises, so a level observed in one
#' cohort and not another cannot enter the comparison.
#'
#' @param designs A list of design tables.
#' @param covariates Candidate covariate columns.
#'
#' @returns A named list, one entry per usable covariate: `NA_character_` for a
#'   continuous covariate, otherwise the shared levels with the reference first.
#'   A covariate that cannot be made comparable has no entry.
#' @keywords internal
#' @noRd
chorale_profile_levels <- function(designs, covariates) {
  out <- list()
  for (cv in covariates) {
    if (!all(vapply(designs, function(d) cv %in% colnames(d), logical(1)))) next
    vals <- lapply(designs, function(d) d[[cv]])
    kind <- chorale_covariate_kind(vals)
    if (is.na(kind)) next
    if (identical(kind, "continuous")) {
      out[[cv]] <- NA_character_
      next
    }
    lv <- Reduce(intersect, lapply(vals, function(v) {
      unique(as.character(stats::na.omit(v)))
    }))
    lv <- sort(lv)
    if (length(lv) < 2) next
    out[[cv]] <- lv
  }
  out
}


#' How a covariate can be compared across modalities
#'
#' One decision, used wherever a covariate is admitted and wherever its terms
#' are built, so the two cannot disagree about the same column.
#'
#' A covariate that is numeric in every modality is continuous and contributes
#' one term, whatever how many distinct values any one modality happens to
#' realise. Deciding otherwise would send it down the categorical path, where
#' levels are matched as text: two modalities recording the same age as `6` and
#' `6.0` would then share no level, and a covariate every modality measures on
#' the same scale would be discarded for a difference of formatting.
#'
#' A covariate that is categorical in every modality contributes one term per
#' shared level beyond a reference. A covariate that is numeric in some
#' modalities and categorical in others cannot be compared at all.
#'
#' @param vals A list holding the covariate as each modality records it.
#' @returns `"continuous"`, `"categorical"`, or `NA_character_` where the
#'   modalities disagree about its type.
#' @keywords internal
#' @noRd
chorale_covariate_kind <- function(vals) {
  if (!length(vals)) return(NA_character_)
  numeric_all <- all(vapply(vals, is.numeric, logical(1)))
  categorical_all <- all(vapply(vals, function(v) !is.numeric(v), logical(1)))
  if (numeric_all) return("continuous")
  if (categorical_all) return("categorical")
  NA_character_
}

#' Column names a set of profile terms produces
#' @keywords internal
#' @noRd
chorale_profile_terms <- function(levels) {
  unlist(lapply(names(levels), function(cv) {
    lv <- levels[[cv]]
    if (length(lv) == 1 && is.na(lv)) return(cv)
    paste0(cv, "=", lv[-1])
  }), use.names = FALSE)
}

#' Signed standardised effect of a contrast on a factor
#' @keywords internal
#' @noRd
chorale_contrast_effect <- function(y, v, reference, level) {
  a <- y[!is.na(v) & v == reference]
  b <- y[!is.na(v) & v == level]
  if (length(a) < 2 || length(b) < 2) return(0)
  pooled <- sqrt((stats::var(a) + stats::var(b)) / 2)
  if (!is.finite(pooled) || pooled == 0) return(0)
  d <- (mean(b) - mean(a)) / pooled
  if (is.finite(d)) d else 0
}

#' Solve the assignment across all modalities at once
#'
#' Comparing modalities a pair at a time and reconciling the results afterwards
#' admits assignments that contradict one another: a transcriptome factor may be
#' assigned to one proteome factor directly and to a different one by way of the
#' metabolome. Merging the two afterwards propagates a single wrong pairing
#' through the whole programme.
#'
#' The assignment is therefore solved once over the collection. Every pairwise
#' cosine between design profiles is placed in one block matrix, and its leading
#' eigenvectors give each factor a position in a common space of latent
#' programmes; every pairwise comparison contributes to every position, so
#' correct evidence reinforces across the collection while inconsistent evidence
#' cancels. Rounding each modality's block back to an assignment against one
#' frame yields correspondences that agree around every cycle by construction,
#' which is what removes the failure mode rather than repairing it.
#'
#' This is permutation synchronisation as formulated by Pachauri, Kondor and
#' Singh (NeurIPS 2013); Chen, Guibas and Huang (ICML 2014) establish that the
#' joint solution recovers the true assignment even when a large fraction of the
#' pairwise inputs behave as outliers, which no independently solved pairwise
#' procedure can do.
#'
#' @param profiles A named list of factors-by-terms profile matrices, one per
#'   modality, with the same terms in the same order.
#'
#' @returns A list with `assignment`, a data frame of one row per
#'   (programme, modality, factor); `similarity`, the pairwise cosine blocks;
#'   and `n_programmes`.
#' @keywords internal
#' @noRd
chorale_synchronise <- function(profiles) {
  mods <- names(profiles)
  k <- vapply(profiles, nrow, integer(1))
  if (length(mods) < 2 || any(k < 1)) {
    return(list(assignment = data.frame(), similarity = list(),
                n_programmes = 0L))
  }
  offset <- cumsum(c(0L, k))
  names(offset) <- c(mods, "")
  total <- sum(k)

  # Profile entries are standardised effects, so their inner product is already
  # comparable across modalities and carries magnitude as well as direction: a
  # factor that barely moves with the design contributes little, which is the
  # behaviour a matching statistic should have.
  similarity <- list()
  for (i in seq_along(mods)) {
    for (j in seq_along(mods)) {
      if (j <= i) next
      block <- profiles[[i]] %*% t(profiles[[j]])
      block[!is.finite(block)] <- 0
      similarity[[paste(mods[i], mods[j], sep = "|")]] <- block
    }
  }

  # The block matrix is put on a common footing with the identity blocks on its
  # diagonal by dividing through by the largest agreement observed; relative
  # magnitudes, which carry the evidence, are untouched.
  scale_by <- max(abs(unlist(similarity)), na.rm = TRUE)
  if (!is.finite(scale_by) || scale_by == 0) scale_by <- 1
  w <- matrix(0, total, total)
  for (i in seq_along(mods)) {
    ri <- offset[i] + seq_len(k[i])
    w[ri, ri] <- diag(k[i])
    for (j in seq_along(mods)) {
      if (j <= i) next
      rj <- offset[j] + seq_len(k[j])
      block <- abs(similarity[[paste(mods[i], mods[j], sep = "|")]]) / scale_by
      w[ri, rj] <- block
      w[rj, ri] <- t(block)
    }
  }

  d <- min(k)
  e <- eigen(w, symmetric = TRUE)
  vals <- pmax(e$values[seq_len(d)], 0)
  embedding <- e$vectors[, seq_len(d), drop = FALSE] %*% diag(sqrt(vals), d, d)

  blocks <- lapply(seq_along(mods), function(i) {
    embedding[offset[i] + seq_len(k[i]), , drop = FALSE]
  })
  names(blocks) <- mods

  # The frame only names the programmes; the evidence placing every factor in
  # the common space came from the whole collection. The modality carrying the
  # fewest factors defines it, so every programme slot is realisable, ties
  # broken towards the block with the most mass in the leading subspace.
  leverage <- vapply(blocks, function(b) sum(b^2), numeric(1))
  candidates <- which(k == min(k))
  anchor <- candidates[which.max(leverage[candidates])]

  reference <- blocks[[anchor]]
  rows <- list()
  for (i in seq_along(mods)) {
    stat <- abs(blocks[[i]] %*% t(reference))
    stat[!is.finite(stat)] <- 0
    pick <- chorale_assign(stat)
    for (f in seq_len(k[i])) {
      slot <- pick[f]
      if (is.na(slot)) next
      rows[[length(rows) + 1]] <- data.frame(
        programme = slot,
        modality = mods[i],
        factor = rownames(profiles[[i]])[f],
        factor_index = f,
        stringsAsFactors = FALSE
      )
    }
  }

  assignment <- do.call(rbind, rows)
  list(assignment = assignment, similarity = similarity,
       n_programmes = if (is.null(assignment)) 0L else
         length(unique(assignment$programme)))
}

#' Agreement between two members of a candidate programme
#' @keywords internal
#' @noRd
chorale_pair_agreement <- function(similarity, ma, mb, fa, fb) {
  key <- paste(ma, mb, sep = "|")
  if (!is.null(similarity[[key]])) return(similarity[[key]][fa, fb])
  key <- paste(mb, ma, sep = "|")
  if (!is.null(similarity[[key]])) return(similarity[[key]][fb, fa])
  NA_real_
}

#' Joint statistic of one candidate programme
#'
#' The mean agreement over every pair of members, evaluated as one quantity, so
#' a programme is judged by what all its modalities do together rather than by
#' the strongest pair it happens to contain.
#'
#' @keywords internal
#' @noRd
chorale_programme_statistic <- function(similarity, members) {
  n <- nrow(members)
  if (n < 2) return(NA_real_)
  vals <- numeric(0)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (j <= i) next
      vals <- c(vals, abs(chorale_pair_agreement(
        similarity, members$modality[i], members$modality[j],
        members$factor_index[i], members$factor_index[j]
      )))
    }
  }
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
}

#' Every subset of a candidate programme carrying at least two modalities
#' @keywords internal
#' @noRd
chorale_member_subsets <- function(n) {
  if (n > 6L) {
    rlang::abort(
      paste0("Enumerating modality subsets is only supported up to six ",
             "modalities; ", n, " were supplied."),
      class = "chorale_too_many_modalities"
    )
  }
  grid <- expand.grid(rep(list(c(FALSE, TRUE)), n))
  keep <- rowSums(grid) >= 2
  lapply(which(keep), function(r) which(as.logical(grid[r, ])))
}

#' Best joint statistic any programme and any modality subset reaches
#'
#' The null has to see the same freedom the estimator had. The synchronisation
#' is solved again on permuted labels, every candidate programme it produces is
#' scored, and every modality subset of every candidate is scored, so the
#' calibrating value is the best a complete run of the procedure could reach by
#' chance rather than the best a single pre-chosen tuple could reach.
#'
#' @keywords internal
#' @noRd
chorale_best_statistic <- function(sync) {
  if (sync$n_programmes == 0) return(0)
  best <- 0
  for (pr in unique(sync$assignment$programme)) {
    members <- sync$assignment[sync$assignment$programme == pr, , drop = FALSE]
    if (nrow(members) < 2) next
    for (idx in chorale_member_subsets(nrow(members))) {
      s <- chorale_programme_statistic(sync$similarity,
                                       members[idx, , drop = FALSE])
      if (is.finite(s) && s > best) best <- s
    }
  }
  best
}
