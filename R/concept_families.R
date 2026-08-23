#' Groups of concepts that describe the same biology
#'
#' Curated vocabularies name the same process at several resolutions, so a
#' change spread across a dozen overlapping concepts can leave every one of them
#' below threshold while the group as a whole is unambiguous. Testing concepts
#' one at a time cannot see that, because each is compared with the largest
#' statistic anywhere in the vocabulary and none of them is large.
#'
#' Families are formed from the vocabulary itself, by how far the concepts share
#' members, and never from the data. A family is therefore fixed before anything
#' is tested, exactly as the vocabulary is.
#'
#' @param concepts A `chorale_concepts` object.
#' @param min_overlap Jaccard overlap at which two concepts join the same
#'   family under complete linkage, so every pair within a family overlaps at
#'   least this much.
#' @param min_size Smallest family reported. A concept in no family of this size
#'   is reported as its own family of one, so nothing is dropped.
#'
#' @returns A data frame with one row per concept, carrying `concept`, `family`
#'   and `family_size`.
#' @export
#' @examples
#' sets <- list(a = letters[1:10], b = letters[2:11], c = LETTERS[1:10])
#' ids <- c(letters, LETTERS)
#' m <- matrix(rnorm(length(ids) * 8), nrow = length(ids),
#'             dimnames = list(ids, paste0("s", 1:8)))
#' design <- data.frame(sample_id = paste0("s", 1:8),
#'                      phenotype = rep(c("case", "control"), 4))
#' cc <- chorale_concepts(list(chorale_load(m, design)), sets, min_features = 3)
#' chorale_concept_families(cc, min_overlap = 0.5)
chorale_concept_families <- function(concepts, min_overlap = 0.25,
                                     min_size = 2L) {
  if (!inherits(concepts, "chorale_concepts")) {
    rlang::abort("`concepts` must be a chorale_concepts object.")
  }
  if (min_overlap <= 0 || min_overlap > 1) {
    rlang::abort("`min_overlap` must lie in (0, 1].")
  }
  sets <- concepts$sets
  names_kept <- names(sets)
  n <- length(sets)
  if (n == 0) {
    return(data.frame(concept = character(), family = character(),
                      family_size = integer(), stringsAsFactors = FALSE))
  }
  if (n == 1) {
    return(data.frame(concept = names_kept, family = names_kept,
                      family_size = 1L, stringsAsFactors = FALSE))
  }

  members <- lapply(sets, unique)
  sizes <- lengths(members)
  distance <- matrix(1, n, n, dimnames = list(names_kept, names_kept))
  for (i in seq_len(n - 1L)) {
    for (j in seq(i + 1L, n)) {
      shared <- length(intersect(members[[i]], members[[j]]))
      union_size <- sizes[i] + sizes[j] - shared
      jaccard <- if (union_size > 0) shared / union_size else 0
      distance[i, j] <- distance[j, i] <- 1 - jaccard
    }
  }
  diag(distance) <- 0

  tree <- stats::hclust(stats::as.dist(distance), method = "complete")
  cut <- stats::cutree(tree, h = 1 - min_overlap)
  size <- table(cut)
  family <- ifelse(size[as.character(cut)] >= min_size,
                   paste0("family_", sprintf("%04d", cut)), names_kept)
  out <- data.frame(concept = names_kept, family = unname(family),
                    stringsAsFactors = FALSE)
  counts <- table(out$family)
  out$family_size <- as.integer(counts[out$family])
  out[order(-out$family_size, out$family, out$concept), , drop = FALSE]
}

#' Whether a family of concepts moves with the phenotype as a group
#'
#' The statistic is the sum of the combined per-concept statistics across the
#' family, standardised by the family's size. It is directional, so concepts
#' pointing the same way accumulate and concepts pointing opposite ways cancel:
#' a family is supported when its members agree, not when any one of them is
#' large.
#'
#' Calibration reuses the permutations [chorale_concept_evidence()] already
#' drew. That matters more here than anywhere else: concepts within a family
#' overlap heavily and are therefore strongly dependent, and only a null that
#' carries that dependence gives the sum a correct distribution. Drawing fresh
#' permutations per family would treat dependent concepts as independent and
#' report families as supported when they are not.
#'
#' @param evidence A `chorale_concept_evidence` object.
#' @param families A data frame from [chorale_concept_families()].
#' @param control A [chorale_control()] object.
#' @param ... Named overrides applied to `control`.
#'
#' @returns An object of class `chorale_family_evidence` with `families`, one
#'   row per family per phenotype term carrying the family statistic, the
#'   permutation p-values, the error control and the largest single concept
#'   statistic inside the family; and `membership`, the source concept-to-family
#'   mapping retained for reporting and audit.
#' @export
chorale_family_evidence <- function(evidence, families,
                                    control = chorale_control(), ...) {
  if (!inherits(evidence, "chorale_concept_evidence")) {
    rlang::abort("`evidence` must be a chorale_concept_evidence object.")
  }
  control <- chorale_merge_control(control, list(...))
  directional_null <- evidence$null_signed_by_concept
  if (is.null(directional_null)) {
    rlang::abort(paste0("`evidence` carries no signed per-concept null. Refit ",
                        "with chorale_concept_evidence() from this version."),
                 class = "chorale_missing_concept_null")
  }

  joint <- evidence$joint
  signed <- joint$joint_z
  keys <- joint$key

  family_of <- stats::setNames(families$family, families$concept)
  rows <- list()
  observed_by_term <- split(seq_along(keys), joint$term)

  for (term in names(observed_by_term)) {
    at <- observed_by_term[[term]]
    fam <- family_of[joint$concept[at]]
    fam[is.na(fam)] <- joint$concept[at][is.na(fam)]
    columns <- match(keys[at], colnames(directional_null))

    groups <- split(seq_along(at), fam)
    statistic <- vapply(groups, function(g) {
      abs(sum(signed[at][g], na.rm = TRUE)) / sqrt(length(g))
    }, numeric(1))
    null_statistic <- vapply(groups, function(g) {
      cols <- columns[g]
      cols <- cols[!is.na(cols)]
      if (length(cols) == 0) return(rep(NA_real_, nrow(directional_null)))
      abs(rowSums(directional_null[, cols, drop = FALSE])) / sqrt(length(cols))
    }, numeric(nrow(directional_null)))
    if (is.null(dim(null_statistic))) {
      null_statistic <- matrix(null_statistic, nrow = nrow(directional_null))
    }
    null_max <- apply(null_statistic, 1, max, na.rm = TRUE)

    rows[[term]] <- data.frame(
      family = names(groups),
      term = term,
      n_concepts = lengths(groups),
      statistic = round(statistic, 4),
      mean_z = round(vapply(groups, function(g)
        mean(signed[at][g], na.rm = TRUE), numeric(1)), 4),
      best_single_z = round(vapply(groups, function(g)
        max(abs(signed[at][g]), na.rm = TRUE), numeric(1)), 4),
      p_value = vapply(seq_along(groups), function(i) {
        (1 + sum(null_statistic[, i] >= statistic[i], na.rm = TRUE)) /
          (1 + nrow(null_statistic))
      }, numeric(1)),
      p_family = vapply(statistic, function(s) {
        (1 + sum(null_max >= s, na.rm = TRUE)) / (1 + length(null_max))
      }, numeric(1)),
      q_value = chorale_permutation_fdr(statistic, null_statistic),
      members = vapply(groups, function(g)
        paste(joint$concept[at][g], collapse = "; "), character(1)),
      stringsAsFactors = FALSE)
  }

  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(family = character(), term = character(),
               n_concepts = integer(), statistic = numeric(),
               mean_z = numeric(), best_single_z = numeric(),
               p_value = numeric(), p_family = numeric(), q_value = numeric(),
               members = character(), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out$significant <- out$q_value <= control$alpha
  out$family_significant <- out$p_family <= control$alpha
  out <- out[order(-out$statistic), , drop = FALSE]
  rownames(out) <- NULL

  structure(list(families = out, membership = families,
                 alpha = control$alpha,
                 n_permutations = evidence$n_permutations),
            class = "chorale_family_evidence")
}

#' @export
print.chorale_family_evidence <- function(x, ...) {
  cat("<chorale_family_evidence>\n")
  if (nrow(x$families) == 0) {
    cat("  no family to test\n")
    return(invisible(x))
  }
  show <- x$families[, setdiff(colnames(x$families), "members"), drop = FALSE]
  print(utils::head(show, 10))
  invisible(x)
}
