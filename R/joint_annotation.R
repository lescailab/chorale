#' What each component of the joint state is made of
#'
#' A component of [chorale_joint_state()] is a direction in concept space. Its
#' label carries no meaning: `joint_05` is the fifth column of a decomposition,
#' and a table reporting that it moves with the phenotype says nothing about
#' which biology moved. The direction is a weighted combination of named
#' concepts, so the biology is in the loadings, and this reports it.
#'
#' Concepts are named in two groups, those the component ranks high and those it
#' ranks low, because a direction is a contrast rather than a list. Where the
#' component's phenotype effect is supplied the direction is oriented by it, so
#' the concepts named as high are the ones higher in cases; where it is not, the
#' orientation a decomposition returns is arbitrary and the table says so.
#'
#' Where a family assignment is supplied, the share of the component's loading
#' mass each family holds is computed and the largest is named with its members.
#' A component whose mass concentrates in one family describes that family; one
#' whose mass spreads over many is a combination the vocabulary has no single
#' name for, and the share reported makes the difference visible.
#'
#' @param state A `chorale_joint_state` object.
#' @param evidence Optional `chorale_joint_evidence` object. Where supplied,
#'   each component is oriented so its effect on the named term is positive, and
#'   the term is recorded.
#' @param families Optional data frame with `concept` and `family` columns, as
#'   returned by [chorale_concept_families()].
#' @param n_top Concepts named in each direction.
#' @param n_members Members named for the leading family.
#' @param term Which term of `evidence` orients the components. `NULL` takes the
#'   first term it carries.
#'
#' @returns A data frame with one row per component carrying the share of
#'   variance it holds, how many concepts load on it, how it was oriented, the
#'   concepts it ranks high and low with their loadings, and, where a family
#'   assignment was supplied, the family holding most of its loading mass.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 80,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' ids <- sprintf("feature_%05d", seq_len(80))
#' sim$modalities <- lapply(sim$modalities, function(m) {
#'   rownames(m) <- ids
#'   m
#' })
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' cc <- chorale_concepts(containers, list(one = ids[1:20], two = ids[21:45]),
#'                        min_features = 5)
#' enc <- chorale_encode(containers, cc, n_free = 1, n_init = 2)
#' state <- chorale_joint_state(enc, n_components = 1)
#' chorale_joint_concepts(state, n_top = 2)
chorale_joint_concepts <- function(state, evidence = NULL, families = NULL,
                                   n_top = 10L, n_members = 5L, term = NULL) {
  if (!inherits(state, "chorale_joint_state")) {
    rlang::abort("`state` must be a chorale_joint_state object.")
  }
  if (!is.null(evidence) && !inherits(evidence, "chorale_joint_evidence")) {
    rlang::abort("`evidence` must be a chorale_joint_evidence object.")
  }
  loadings <- state$loadings
  if (is.null(loadings) || ncol(loadings) == 0 || nrow(loadings) == 0) {
    return(chorale_empty_joint_concepts())
  }
  families <- chorale_check_family_table(families, rownames(loadings))

  orientation <- chorale_joint_orientation(state, evidence, term)
  share <- state$variance$share[match(colnames(loadings),
                                      state$variance$component)]

  rows <- lapply(seq_len(ncol(loadings)), function(k) {
    component <- colnames(loadings)[k]
    v <- loadings[, k]
    v <- v[is.finite(v)]
    sign <- orientation$sign[[component]] %||% 1
    v <- v * sign
    leading <- chorale_leading_family(v, families, n_members)
    data.frame(
      component = component,
      share = share[k],
      n_concepts = sum(v != 0),
      orientation = orientation$label,
      concepts_high = chorale_name_loadings(v, n_top, high = TRUE),
      concepts_low = chorale_name_loadings(v, n_top, high = FALSE),
      leading_family = leading$family,
      leading_family_share = leading$share,
      leading_family_members = leading$members,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(families)) {
    out$leading_family <- NULL
    out$leading_family_share <- NULL
    out$leading_family_members <- NULL
  }
  rownames(out) <- NULL
  out
}

#' The empty annotation, so a state with no component returns the same columns
#' @keywords internal
#' @noRd
chorale_empty_joint_concepts <- function() {
  data.frame(component = character(), share = numeric(),
             n_concepts = integer(), orientation = character(),
             concepts_high = character(), concepts_low = character(),
             stringsAsFactors = FALSE)
}

#' Check a family assignment against the vocabulary the loadings carry
#' @keywords internal
#' @noRd
chorale_check_family_table <- function(families, vocabulary) {
  if (is.null(families)) return(NULL)
  if (!is.data.frame(families) ||
      !all(c("concept", "family") %in% colnames(families))) {
    rlang::abort("`families` must carry `concept` and `family` columns.")
  }
  families <- families[families$concept %in% vocabulary, , drop = FALSE]
  if (nrow(families) == 0) {
    rlang::abort("`families` names no concept the joint state carries.")
  }
  families
}

#' How each component is oriented, and what that orientation rests on
#'
#' A decomposition fixes neither the sign nor the order of its components, so
#' the concepts a component ranks high are the ones it ranks high only up to a
#' sign. Where the phenotype effect is known the sign is fixed by it and the
#' two groups of concepts mean something; where it is not, the table says the
#' orientation is arbitrary rather than letting it be read as a direction.
#'
#' @keywords internal
#' @noRd
chorale_joint_orientation <- function(state, evidence, term) {
  arbitrary <- list(sign = list(), label = "arbitrary")
  if (is.null(evidence) || nrow(evidence$components) == 0) return(arbitrary)
  components <- evidence$components
  term <- term %||% components$term[1]
  components <- components[components$term == term, , drop = FALSE]
  if (nrow(components) == 0) {
    rlang::abort(paste0("`evidence` carries no term `", term, "`."))
  }
  sign <- as.list(ifelse(components$effect < 0, -1, 1))
  names(sign) <- components$component
  list(sign = sign, label = paste0("high = raised in ", term))
}

#' Name the concepts at one end of a component, with their loadings
#' @keywords internal
#' @noRd
chorale_name_loadings <- function(v, n_top, high) {
  v <- if (high) v[v > 0] else v[v < 0]
  if (length(v) == 0) return(NA_character_)
  v <- sort(v, decreasing = high)
  v <- v[seq_len(min(as.integer(n_top), length(v)))]
  paste0(names(v), " (", format(round(v, 3), trim = TRUE), ")",
         collapse = "; ")
}

#' The family holding most of a component's loading mass
#' @keywords internal
#' @noRd
chorale_leading_family <- function(v, families, n_members) {
  absent <- list(family = NA_character_, share = NA_real_,
                 members = NA_character_)
  if (is.null(families)) return(absent)
  # Loading mass is the squared loading, so it is the concept's contribution to
  # the component's norm and a share of it sums to one across the vocabulary.
  # Using the absolute loading instead would give a quantity whose shares depend
  # on how many concepts load weakly.
  weight <- v^2
  total <- sum(weight)
  if (!is.finite(total) || total <= 0) return(absent)
  family <- families$family[match(names(v), families$concept)]
  keep <- !is.na(family)
  if (!any(keep)) return(absent)
  by_family <- tapply(weight[keep], family[keep], sum)
  leading <- names(by_family)[which.max(by_family)]
  members <- families$concept[families$family == leading]
  members <- members[order(-weight[match(members, names(v))])]
  named <- members[seq_len(min(as.integer(n_members), length(members)))]
  if (length(members) > length(named)) {
    named <- c(named, paste0("and ", length(members) - length(named), " more"))
  }
  list(family = leading,
       share = round(unname(max(by_family)) / total, 4),
       members = paste(named, collapse = "; "))
}
