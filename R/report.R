#' Emit interpretable outputs from a fit
#'
#' Writes the outputs of `AGENT_PLAN.md` Section 9. Every factor reported
#' carries a pathway definition and its marker features; a factor that cannot
#' be named is written to the same tables marked unresolved, rather than
#' dropped, since a factor that is identified but unnameable is a finding about
#' the gene sets rather than an absence.
#'
#' Cross-modality disagreement is reported alongside agreement. Proteome and
#' transcriptome inconsistency in Alzheimer's brain is enriched in the amyloid
#' plaque microenvironment and reflects amyloid-delayed protein turnover
#' (Yarbro et al., Nat Commun 16:1533, 2025), so factors moving in opposite
#' directions across modalities are ranked and kept rather than discarded.
#'
#' @param fit A `chorale_fit` object, as returned by [chorale_fit()].
#' @param bound A `chorale_bound` object, as returned by [chorale_bound()].
#' @param null A `chorale_null` object, as returned by [chorale_null()].
#' @param path Directory to write into. Created if absent.
#' @param n_top_sets Number of curated sets naming each factor.
#'
#' @returns Invisibly, a character vector of the files written.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2)
#' out <- chorale_report(fit, chorale_bound(fit), NULL,
#'                       path = withr::local_tempdir())
#' basename(out)
chorale_report <- function(fit, bound = NULL, null = NULL, path,
                           n_top_sets = 5L) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  written <- character()

  factors <- chorale_factor_table(fit, n_top_sets = n_top_sets)
  written <- c(written, chorale_write(factors, path, "factors.tsv"))

  markers <- chorale_marker_table(fit)
  written <- c(written, chorale_write(markers, path, "markers.tsv"))

  for (m in fit$modalities) {
    written <- c(written, chorale_write(
      chorale_loading_table(fit, m), path, paste0("loadings_", m, ".tsv")
    ))
    written <- c(written, chorale_write(
      chorale_score_table(fit, m), path, paste0("scores_", m, ".tsv")
    ))
  }

  associations <- chorale_association_table(fit)
  written <- c(written, chorale_write(associations, path, "associations.tsv"))

  concordance <- chorale_concordance_table(fit)
  written <- c(written, chorale_write(concordance, path, "concordance.tsv"))

  bounds_tbl <- if (!is.null(bound)) bound$bounds else data.frame()
  written <- c(written, chorale_write(bounds_tbl, path, "bounds.tsv"))

  controls <- chorale_control_table(fit, null)
  written <- c(written, chorale_write(controls, path, "controls.tsv"))

  written <- c(written, chorale_write_gmt(fit, factors, path))
  written <- c(written, chorale_write_mae(fit, path))
  written <- c(written, chorale_write_html(fit, factors, markers, associations,
                                           concordance, bounds_tbl, controls,
                                           path))
  invisible(written)
}

#' @keywords internal
#' @noRd
chorale_write <- function(x, path, file) {
  f <- file.path(path, file)
  utils::write.table(x, f, sep = "\t", row.names = FALSE, quote = FALSE)
  f
}

#' Pathway definition and resolution status per factor
#' @keywords internal
#' @noRd
chorale_factor_table <- function(fit, n_top_sets = 5L) {
  rows <- list()
  for (m in fit$modalities) {
    f <- fit$fits[[m]]
    sw <- f$set_weights
    shared_factors <- unique(c(
      fit$matches$factor_a[fit$matches$modality_a == m],
      fit$matches$factor_b[fit$matches$modality_b == m]
    ))
    for (j in colnames(f$loadings)) {
      sets <- NA_character_
      scores <- NA_character_
      if (!is.null(sw) && nrow(sw) > 0 && j %in% colnames(sw)) {
        w <- sw[, j]
        w <- w[order(abs(w), decreasing = TRUE)]
        w <- w[abs(w) > 0][seq_len(min(n_top_sets, sum(abs(w) > 0)))]
        if (length(w) > 0) {
          sets <- paste(names(w), collapse = ";")
          scores <- paste(signif(as.numeric(w), 3), collapse = ";")
        }
      }
      has_markers <- length(f$markers[[j]]) >= 2
      named <- !is.na(sets)
      rows[[length(rows) + 1]] <- data.frame(
        modality = m,
        factor = j,
        shared = j %in% shared_factors,
        pathway_definition = sets,
        pathway_weight = scores,
        n_markers = length(f$markers[[j]]),
        purity_margin = signif(f$purity_margin[[j]], 3),
        pure_feature_condition = has_markers,
        status = if (named && has_markers) {
          "resolved"
        } else if (!has_markers) {
          "unresolved: no pure features"
        } else {
          "unresolved: no pathway definition"
        },
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Marker features per factor per modality
#' @keywords internal
#' @noRd
chorale_marker_table <- function(fit) {
  rows <- list()
  for (m in fit$modalities) {
    f <- fit$fits[[m]]
    for (j in names(f$markers)) {
      ms <- f$markers[[j]]
      strict <- length(ms) > 0
      if (!strict) ms <- f$best_candidates[[j]][seq_len(min(2, length(f$best_candidates[[j]])))]
      if (length(ms) == 0) next
      rows[[length(rows) + 1]] <- data.frame(
        modality = m, factor = j, feature = ms,
        loading = signif(f$loadings[ms, j], 4),
        clears_purity_threshold = strict,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame() else out
}

#' @keywords internal
#' @noRd
chorale_loading_table <- function(fit, m) {
  l <- fit$fits[[m]]$loadings
  data.frame(feature = rownames(l), signif(as.data.frame(l), 5),
             check.names = FALSE, stringsAsFactors = FALSE)
}

#' @keywords internal
#' @noRd
chorale_score_table <- function(fit, m) {
  s <- fit$fits[[m]]$scores
  d <- fit$designs[[m]]
  d <- d[match(rownames(s), d$sample_id), , drop = FALSE]
  data.frame(d, signif(as.data.frame(s), 5), check.names = FALSE,
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Association of each factor with the design covariates
#' @keywords internal
#' @noRd
chorale_association_table <- function(fit, n_perm = 999L) {
  rows <- list()
  for (m in fit$modalities) {
    s <- fit$fits[[m]]$scores
    d <- fit$designs[[m]]
    d <- d[match(rownames(s), d$sample_id), , drop = FALSE]
    for (cov in intersect(c("phenotype", "age_months", "sex", "strain"),
                          colnames(d))) {
      v <- d[[cov]]
      if (length(unique(stats::na.omit(v))) < 2) next
      for (j in colnames(s)) {
        y <- s[, j]
        ok <- is.finite(y) & !is.na(v)
        if (sum(ok) < 6) next
        stat <- chorale_effect(y[ok], v[ok])
        set.seed(1)
        null <- replicate(n_perm, chorale_effect(y[ok], sample(v[ok]))$stat)
        p <- (1 + sum(null >= stat$stat)) / (1 + n_perm)
        rows[[length(rows) + 1]] <- data.frame(
          modality = m, factor = j, covariate = cov,
          effect_size = signif(stat$effect, 4),
          statistic = signif(stat$stat, 4),
          p_permutation = p,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out$p_adjusted <- stats::p.adjust(out$p_permutation, method = "BH")
  out[order(out$p_permutation), , drop = FALSE]
}

#' Effect of a covariate on a factor score
#' @keywords internal
#' @noRd
chorale_effect <- function(y, v) {
  if (is.numeric(v) && length(unique(v)) > 2) {
    r <- suppressWarnings(stats::cor(y, v, method = "spearman"))
    if (!is.finite(r)) r <- 0
    return(list(effect = r, stat = abs(r)))
  }
  v <- as.factor(as.character(v))
  if (nlevels(v) == 2) {
    g <- split(y, v)
    # Standardised mean difference, so the effect is comparable across factors.
    pooled <- sqrt(mean(vapply(g, stats::var, numeric(1)), na.rm = TRUE))
    if (!is.finite(pooled) || pooled == 0) pooled <- 1
    eff <- (mean(g[[2]]) - mean(g[[1]])) / pooled
    return(list(effect = eff, stat = abs(eff)))
  }
  fitted <- stats::lm(y ~ v)
  eta <- summary(fitted)$r.squared
  list(effect = eta, stat = eta)
}

#' Cross-modality agreement and the direction of disagreement
#' @keywords internal
#' @noRd
chorale_concordance_table <- function(fit) {
  if (nrow(fit$matches) == 0) return(data.frame())
  rows <- list()
  for (i in seq_len(nrow(fit$matches))) {
    r <- fit$matches[i, ]
    a <- chorale_stratum_means(fit$fits[[r$modality_a]]$scores,
                               fit$designs[[r$modality_a]], fit$strata_keys)
    b <- chorale_stratum_means(fit$fits[[r$modality_b]]$scores,
                               fit$designs[[r$modality_b]], fit$strata_keys)
    common <- intersect(rownames(a), rownames(b))
    if (length(common) < 3) next
    pa <- a[common, r$factor_a]
    pb <- b[common, r$factor_b] * r$sign
    agreement <- suppressWarnings(stats::cor(pa, pb, method = "spearman"))
    rows[[length(rows) + 1]] <- data.frame(
      modality_a = r$modality_a, modality_b = r$modality_b,
      factor_a = r$factor_a, factor_b = r$factor_b,
      agreement = signif(agreement, 4),
      direction = if (r$sign < 0) "opposed" else "aligned",
      n_strata = length(common),
      discordant_strata = sum(sign(pa) != sign(pb)),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out[order(out$direction, -out$agreement), , drop = FALSE]
}

#' Controls accompanying the fit
#' @keywords internal
#' @noRd
chorale_control_table <- function(fit, null) {
  rows <- list()
  if (!is.null(null)) {
    rows[[length(rows) + 1]] <- data.frame(
      control = "phenotype permutation",
      value = signif(null$p_phenotype, 4),
      detail = sprintf("%d permutations, strongest observed agreement %.3f",
                       null$n_permutations, null$observed_agreement),
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      control = "modality shuffle",
      value = if (isTRUE(null$modality_null$applicable)) {
        signif(null$modality_null$agreement, 4)
      } else {
        NA_real_
      },
      detail = if (isTRUE(null$modality_null$applicable)) {
        "agreement recovered after reassigning samples across modalities"
      } else {
        paste("not applicable:", null$modality_null$reason)
      },
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(null$stability))) {
      rows[[length(rows) + 1]] <- data.frame(
        control = paste0("initialisation stability: ", null$stability$modality[i]),
        value = signif(null$stability$objective_cv[i], 4),
        detail = sprintf("%d initialisations, %d failed",
                         null$stability$n_init[i], null$stability$n_failed[i]),
        stringsAsFactors = FALSE
      )
    }
  }
  rows[[length(rows) + 1]] <- data.frame(
    control = "factors clearing the pure-feature condition",
    value = sum(vapply(fit$modalities, function(m) {
      sum(fit$fits[[m]]$pure_feature_condition)
    }, numeric(1))),
    detail = paste("of", sum(fit$n_factors), "factors across all modalities"),
    stringsAsFactors = FALSE
  )
  do.call(rbind, rows)
}

#' Factor definitions in GMT form
#' @keywords internal
#' @noRd
chorale_write_gmt <- function(fit, factors, path) {
  f <- file.path(path, "factors.gmt")
  lines <- character()
  for (m in fit$modalities) {
    for (j in names(fit$fits[[m]]$markers)) {
      ms <- fit$fits[[m]]$markers[[j]]
      if (length(ms) == 0) next
      row <- factors[factors$modality == m & factors$factor == j, , drop = FALSE]
      desc <- if (nrow(row) > 0 && !is.na(row$pathway_definition[1])) {
        row$pathway_definition[1]
      } else {
        "unresolved"
      }
      lines <- c(lines, paste(c(paste0(m, ":", j), desc, ms), collapse = "\t"))
    }
  }
  writeLines(lines, f)
  f
}

#' The full model as a MultiAssayExperiment
#'
#' Stores the per-modality assays with their sample metadata and carries the
#' fitted scores, loadings, markers and matches as metadata, so a reader can
#' reopen the model without rerunning it.
#'
#' @keywords internal
#' @noRd
chorale_write_mae <- function(fit, path) {
  f <- file.path(path, "chorale_mae.rds")
  if (!rlang::is_installed("MultiAssayExperiment")) {
    saveRDS(list(fits = fit$fits, designs = fit$designs,
                 matches = fit$matches), f)
    return(f)
  }
  experiments <- lapply(fit$modalities, function(m) {
    l <- fit$fits[[m]]$loadings
    SummarizedExperiment::SummarizedExperiment(
      assays = list(loadings = l)
    )
  })
  names(experiments) <- fit$modalities
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experiments,
    metadata = list(
      scores = lapply(fit$fits, `[[`, "scores"),
      markers = lapply(fit$fits, `[[`, "markers"),
      designs = fit$designs,
      matches = fit$matches,
      n_factors = fit$n_factors,
      strata_keys = fit$strata_keys
    )
  )
  saveRDS(mae, f)
  f
}

#' Self-contained report
#' @keywords internal
#' @noRd
chorale_write_html <- function(fit, factors, markers, associations,
                               concordance, bounds, controls, path) {
  f <- file.path(path, "report.html")
  esc <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  tbl <- function(d, caption) {
    if (is.null(d) || nrow(d) == 0) {
      return(paste0("<h2>", esc(caption), "</h2><p>No rows.</p>"))
    }
    head_row <- paste0("<tr>", paste0("<th>", esc(colnames(d)), "</th>",
                                      collapse = ""), "</tr>")
    body <- apply(d, 1, function(r) {
      paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>")
    })
    paste0("<h2>", esc(caption), "</h2><table>", head_row,
           paste(body, collapse = ""), "</table>")
  }

  cards <- ""
  for (i in seq_len(nrow(factors))) {
    r <- factors[i, ]
    ms <- markers[markers$modality == r$modality & markers$factor == r$factor, ,
                  drop = FALSE]
    cards <- paste0(
      cards,
      "<div class='card'><h3>", esc(r$modality), " / ", esc(r$factor), "</h3>",
      "<p class='status ", if (r$status == "resolved") "ok" else "warn", "'>",
      esc(r$status), if (isTRUE(r$shared)) " &middot; shared" else "", "</p>",
      "<p><strong>Pathway definition:</strong> ",
      esc(ifelse(is.na(r$pathway_definition), "none", r$pathway_definition)),
      "</p><p><strong>Markers:</strong> ",
      esc(if (nrow(ms) > 0) paste(ms$feature, collapse = ", ") else "none"),
      "</p><p><strong>Purity margin:</strong> ", esc(r$purity_margin),
      "</p></div>"
    )
  }

  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<title>chorale report</title><style>",
    "body{font-family:system-ui,sans-serif;margin:2rem;max-width:1100px;line-height:1.5}",
    "table{border-collapse:collapse;margin:1rem 0;font-size:.85rem;width:100%;display:block;overflow-x:auto}",
    "th,td{border:1px solid #ccc;padding:.3rem .5rem;text-align:left}",
    "th{background:#f2f2f2}",
    ".card{border:1px solid #ddd;border-radius:6px;padding:.75rem 1rem;margin:.5rem 0}",
    ".status{font-size:.8rem;padding:.1rem .4rem;border-radius:3px;display:inline-block}",
    ".ok{background:#e6f4ea}.warn{background:#fdecea}",
    "</style></head><body>",
    "<h1>chorale report</h1>",
    "<p>Modalities: ", esc(paste(fit$modalities, collapse = ", ")),
    ". Factors per modality: ", esc(paste(fit$n_factors, collapse = ", ")),
    ". Matched shared pairs: ", nrow(fit$matches), ".</p>",
    "<h2>Factor cards</h2>", cards,
    tbl(controls, "Controls"),
    tbl(concordance, "Cross-modality concordance"),
    tbl(utils::head(associations, 50), "Design associations (top 50)"),
    tbl(bounds, "Identified set per coupling"),
    tbl(factors, "Factors"),
    "</body></html>"
  )
  writeLines(html, f)
  f
}
