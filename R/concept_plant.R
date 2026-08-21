#' Score a fit against the concepts that were planted in it
#'
#' On simulated data the answer is known, which no real cohort supplies: a
#' recorded cohort does not say which named concepts take part in the disease.
#' Two questions decide whether the estimator recovered the planting. Were the
#' planted concepts called? And did the rest of the vocabulary stay quiet? A
#' method that calls everything answers the first and fails the second, so both
#' are reported and neither on its own is the result.
#'
#' Ranking is reported beside calling, because a threshold is a choice and a
#' ranking is not. The probability that a planted concept outranks one that was
#' not planted says how far the two are separated whatever threshold is applied.
#'
#' @param fit A `chorale_concept_fit`.
#' @param planted The concepts that were planted: a character vector of names,
#'   or the result of [chorale_plant()], from which the planted set names are
#'   taken.
#' @param alpha Threshold a concept is called at. Defaults to the fit's own.
#' @param criterion Which control decides a call: `"fdr"` holds the false
#'   discovery rate over the vocabulary at `alpha`, `"family_wise"` holds the
#'   probability of any false call at `alpha`. The two answer different
#'   questions and a recovery figure means little without saying which was
#'   applied.
#'
#' @returns A list with `per_concept`, one row per concept in the vocabulary
#'   carrying whether it was planted, its statistic, its rank and whether it was
#'   called; and `summary`, a one-row data frame carrying recall among the
#'   planted concepts, the false positive rate among the rest, and the
#'   probability that a planted concept outranks one that was not.
#' @export
#' @examples
#' fx <- chorale_concept_example(seed = 1)
#' fit <- chorale_concept_fit(fx$containers, fx$sets, n_free = 0,
#'                            n_permutations = 99)
#' chorale_score_concepts(fit, fx$planted)$summary
chorale_score_concepts <- function(fit, planted, alpha = NULL,
                                   criterion = c("fdr", "family_wise")) {
  criterion <- match.arg(criterion)
  if (!inherits(fit, "chorale_concept_fit")) {
    rlang::abort("`fit` must be a chorale_concept_fit object.")
  }
  planted <- chorale_planted_names(planted)
  alpha <- alpha %||% fit$evidence$alpha
  j <- fit$evidence$joint
  if (nrow(j) == 0) {
    return(list(per_concept = data.frame(),
                summary = data.frame(n_planted = length(planted),
                                     recall = NA_real_,
                                     false_positive_rate = NA_real_,
                                     rank_auc = NA_real_,
                                     stringsAsFactors = FALSE)))
  }

  statistic <- abs(j$joint_z)
  per_concept <- data.frame(
    concept = j$concept,
    planted = j$concept %in% planted,
    n_modalities = j$n_modalities,
    statistic = round(statistic, 4),
    rank = rank(-statistic, ties.method = "min"),
    q_value = j$q_value,
    p_family = j$p_family,
    criterion = criterion,
    stringsAsFactors = FALSE
  )
  per_concept$called <- if (identical(criterion, "fdr")) {
    per_concept$q_value <= alpha
  } else {
    per_concept$p_family <= alpha
  }
  per_concept <- per_concept[order(per_concept$rank), , drop = FALSE]
  rownames(per_concept) <- NULL

  is_planted <- per_concept$planted
  summary <- data.frame(
    n_concepts = nrow(per_concept),
    n_planted = sum(is_planted),
    recall = if (any(is_planted)) {
      round(mean(per_concept$called[is_planted]), 3)
    } else {
      NA_real_
    },
    false_positive_rate = if (any(!is_planted)) {
      round(mean(per_concept$called[!is_planted]), 3)
    } else {
      NA_real_
    },
    rank_auc = round(chorale_rank_auc(per_concept$statistic, is_planted), 3),
    criterion = criterion,
    stringsAsFactors = FALSE
  )
  list(per_concept = per_concept, summary = summary)
}

#' The names of the planted concepts, however the planting is supplied
#' @keywords internal
#' @noRd
chorale_planted_names <- function(planted) {
  if (is.character(planted)) return(planted)
  if (is.list(planted)) {
    if (!is.null(planted$concepts)) return(as.character(planted$concepts))
    if (!is.null(planted$planted)) return(unname(as.character(planted$planted)))
  }
  rlang::abort("`planted` must be concept names or the result of chorale_plant().")
}

#' Recovery of planted concepts across a grid of regimes
#'
#' A single simulation says whether the estimator worked once. A curve says
#' where it stops working. Each row of the grid is a regime, several collections
#' are drawn from it, and the recovery of the planted concepts is averaged over
#' them.
#'
#' Three dials matter and each asks a different question. **Sample size** asks
#' whether bulk cohorts are large enough. **Effect size** asks how strong a
#' concept's separation of cases from controls has to be. **Vocabulary
#' coverage** asks how far the vocabulary has to reach into a modality's
#' measured features, which is the question a modality with an incomplete
#' annotation raises.
#'
#' @param grid A data frame of regimes. Any of `n_samples`, `n_features`,
#'   `effect`, `coverage`, `n_modalities`, `n_concepts` and `n_planted` may be
#'   columns; those absent take their default. A `label` column is carried
#'   through.
#' @param n_rep Collections drawn per regime.
#' @param n_permutations Permutations calibrating each fit.
#' @param n_free Free dimensions per modality.
#' @param alpha Threshold a concept is called at.
#' @param criterion Which control decides a call, as in
#'   [chorale_score_concepts()].
#' @param seed Integer seed; replicate `r` of regime `i` uses
#'   `seed + 100 * i + r`.
#'
#' @returns `grid` with the mean recovery over the replicates joined on, the
#'   number of replicates that could be evaluated, and the reason where any
#'   could not.
#' @export
#' @examples
#' grid <- data.frame(label = c("strong", "absent"), effect = c(1, 0))
#' chorale_validate_concepts(grid, n_rep = 1, n_permutations = 99)
chorale_validate_concepts <- function(grid, n_rep = 3L, n_permutations = 999L,
                                      n_free = 0L, alpha = 0.05,
                                      criterion = c("fdr", "family_wise"),
                                      seed = 1L) {
  criterion <- match.arg(criterion)
  if (!is.data.frame(grid) || nrow(grid) == 0) {
    rlang::abort("`grid` must be a data frame with at least one regime.")
  }
  default <- list(n_samples = 60L, n_features = 90L, effect = 1,
                  coverage = 1, n_modalities = 2L, n_concepts = 3L,
                  n_planted = 1L)
  get_col <- function(row, name) {
    if (name %in% names(grid)) grid[[name]][row] else default[[name]]
  }

  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    metrics <- vector("list", n_rep)
    failures <- character()
    for (r in seq_len(n_rep)) {
      s <- seed + 100L * i + r
      fx <- chorale_concept_example(
        n_samples = as.integer(get_col(i, "n_samples")),
        n_features = as.integer(get_col(i, "n_features")),
        effect = get_col(i, "effect"),
        n_modalities = as.integer(get_col(i, "n_modalities")),
        n_concepts = as.integer(get_col(i, "n_concepts")),
        n_planted = as.integer(get_col(i, "n_planted")),
        coverage = get_col(i, "coverage"),
        seed = s)
      fit <- try(chorale_concept_fit(
        fx$containers, fx$sets, n_free = n_free,
        n_permutations = n_permutations, alpha = alpha, seed = s),
        silent = TRUE)
      if (inherits(fit, "try-error")) {
        # A regime the estimator cannot be run in is a result about the regime,
        # so the reason travels with the row rather than leaving an empty cell.
        failures <- c(failures, conditionMessage(attr(fit, "condition")))
        next
      }
      metrics[[r]] <- chorale_score_concepts(fit, fx$planted, alpha,
                                             criterion = criterion)$summary
    }
    metrics <- do.call(rbind, metrics)
    agg <- if (is.null(metrics)) {
      data.frame(recall = NA_real_, false_positive_rate = NA_real_,
                 rank_auc = NA_real_, n_evaluated = 0L,
                 reason = failures[1] %||% NA_character_,
                 stringsAsFactors = FALSE)
    } else {
      data.frame(
        recall = round(mean(metrics$recall, na.rm = TRUE), 3),
        false_positive_rate = round(mean(metrics$false_positive_rate,
                                         na.rm = TRUE), 3),
        rank_auc = round(mean(metrics$rank_auc, na.rm = TRUE), 3),
        n_evaluated = nrow(metrics),
        reason = if (length(failures)) failures[1] else NA_character_,
        stringsAsFactors = FALSE)
    }
    out[[i]] <- cbind(grid[i, , drop = FALSE], agg, n_rep = n_rep)
  }
  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}
