#' Emit interpretable outputs from a fit
#'
#' Writes the complete fit, controls and interpretation tables. Every factor
#' carries a pathway definition and its marker features; a factor that cannot
#' be named is written to the same tables marked unresolved, rather than
#' dropped, since an unnamed factor is a finding about
#' the gene sets rather than an absence.
#'
#' Cross-modality disagreement is reported alongside agreement. Factors moving
#' in opposite directions are aligned for comparison and their orientation is
#' retained rather than silently discarded.
#'
#' @param fit A fit to write out: a `chorale_concept_fit` or a `chorale_fit`.
#' @param bound A `chorale_bound` object, as returned by [chorale_bound()].
#'   Applies to a `chorale_fit` only.
#' @param null A `chorale_null` object, as returned by [chorale_null()].
#' @param path Directory to write into. Created if absent.
#' @param n_top_sets Number of curated sets naming each factor.
#' @param ... Passed to the method.
#'
#' @returns Invisibly, a character vector of the files written.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2,
#'                    n_ambiguity_boot = 19)
#' out <- chorale_report(fit, chorale_bound(fit), NULL,
#'                       path = withr::local_tempdir())
#' basename(out)
chorale_report <- function(fit, ...) {
  UseMethod("chorale_report")
}

#' @export
chorale_report.default <- function(fit, ...) {
  rlang::abort("`fit` must be a chorale_concept_fit or a chorale_fit object.")
}

#' @rdname chorale_report
#' @export
chorale_report.chorale_fit <- function(fit, bound = NULL, null = NULL, path,
                                       n_top_sets = 5L, ...) {
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

  written <- c(written, chorale_write(fit$matches, path, "matches.tsv"))
  written <- c(written, chorale_write(
    fit$excluded_covariates %||% data.frame(), path, "excluded_covariates.tsv"))

  written <- c(written, chorale_write(
    chorale_integrated_table(fit, factors, markers, associations, concordance),
    path, "programmes.tsv"))
  # Each of these summarises the whole fit, so they are computed once here and
  # handed to the writer rather than recomputed for the rendered report.
  derived <- list(
    leave_one_out = chorale_leave_one_out(fit),
    pathway = chorale_pathway_table(fit),
    fdr = chorale_fdr(fit, associations = associations),
    added_value = chorale_added_value(fit, associations = associations),
    cohorts = chorale_cohort_overlap(fit)
  )
  written <- c(written, chorale_write(
    derived$leave_one_out, path, "leave_one_out.tsv"))
  written <- c(written, chorale_write(
    derived$pathway, path, "pathway_evidence.tsv"))
  written <- c(written, chorale_write(
    derived$fdr, path, "false_discovery.tsv"))
  written <- c(written, chorale_write(
    derived$added_value, path, "added_value.tsv"))
  written <- c(written, chorale_write(
    derived$cohorts, path, "cohort_overlap.tsv"))
  written <- c(written, chorale_write_gmt(fit, factors, path))
  written <- c(written, chorale_write_mae(fit, path))
  written <- c(written, chorale_write_html(fit, factors, markers, associations,
                                           concordance, bounds_tbl, controls,
                                           path, derived = derived))
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
    # A factor counts as shared only where the programme carrying it survived
    # its null. An assignment on its own is not evidence of sharing.
    supported <- if (is.data.frame(fit$matches) && nrow(fit$matches) > 0) {
      fit$matches[fit$matches$significant, , drop = FALSE]
    } else {
      fit$matches
    }
    shared_factors <- if (is.data.frame(supported) && nrow(supported) > 0) {
      unique(c(
        supported$factor_a[supported$modality_a == m],
        supported$factor_b[supported$modality_b == m]
      ))
    } else {
      character(0)
    }
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
        modality = m, factor = j,
        feature = chorale_readable_features(ms),
        feature_id = ms,
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
    profile <- fit$design_profiles[[m]]
    if (is.null(profile) || !ncol(profile$effects)) next
    for (j in rownames(profile$effects)) {
      for (term in colnames(profile$effects)) {
        effect <- profile$effects[j, term]
        se <- profile$se[j, term]
        z <- effect / se
        p <- if (is.finite(z)) 2 * stats::pnorm(-abs(z)) else NA_real_
        rows[[length(rows) + 1]] <- data.frame(
          modality = m, factor = j,
          covariate = unname(profile$term_covariate[term]),
          contrast = term,
          effect_size = signif(effect, 4),
          standard_error = signif(se, 4),
          statistic = signif(abs(z), 4),
          p_permutation = p,
          inference = "adjusted Wald diagnostic; programme support uses FWER permutations",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  # Retained as a compatibility column. Search-adjusted inference lives in the
  # phenotype-led programme test rather than a second, mismatched BH layer.
  out$p_adjusted <- out$p_permutation
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
      detail = sprintf("%d permutations, strongest joint evidence %.3f",
                       null$n_permutations, null$observed_agreement),
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      control = "modality shuffle",
      value = if (isTRUE(null$modality_null$applicable)) {
        signif(null$modality_null$p_value, 4)
      } else {
        NA_real_
      },
      detail = if (isTRUE(null$modality_null$applicable)) {
        sprintf(paste("p-value of the observed evidence against %d shuffles",
                      "reassigning samples across modalities"),
                null$modality_null$n_shuffles)
      } else {
        paste("not applicable:", null$modality_null$reason)
      },
      stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(null$stability))) {
      rows[[length(rows) + 1]] <- data.frame(
        control = paste0("factor stability: ", null$stability$modality[i]),
        value = signif(null$stability$subspace_agreement[i], 4),
        detail = sprintf(paste("mean agreement of the recovered factors across",
                               "%d initialisations, weakest %.3f; %d failed"),
                         null$stability$n_init[i], null$stability$subspace_min[i],
                         null$stability$n_failed[i]),
        stringsAsFactors = FALSE
      )
    }
  } else {
    rows[[length(rows) + 1]] <- data.frame(
      control = "external phenotype-refit control", value = NA_real_,
      detail = "not run; phenotype support inside the fit used its fixed-decomposition null",
      stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      control = "modality shuffle", value = NA_real_, detail = "not run",
      stringsAsFactors = FALSE)
  }
  for (m in fit$modalities) {
    rec <- fit$fits[[m]]$reconstruction
    if (is.null(rec)) next
    rows[[length(rows) + 1]] <- data.frame(
      control = paste0("curated projection: ", m),
      value = signif(rec$projected, 4),
      detail = sprintf(
        paste("variance the curated sets explain, against %.3f for the fitted",
              "loadings; the sets describe the factors, they do not replace",
              "them"),
        rec$fitted),
      stringsAsFactors = FALSE
    )
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

#' Integrated view: what the shared factors are, across modalities
#'
#' One row per shared factor pair, carrying the biology, the markers each
#' modality contributes, and the phenotype effect each modality reports. This is
#' the result integration produces and a single modality cannot.
#'
#' @keywords internal
#' @noRd
chorale_integrated_table <- function(fit, factors, markers, associations,
                                    concordance) {
  # Only programmes whose joint evidence beat its null are reported here. An
  # assignment the null did not reject is not a result, and substituting one
  # would render a negative analysis as a positive integration.
  cols <- c("programme", "n_modalities", "modalities", "programme_pathway",
            "modality", "factor", "markers", "phenotype_effect", "phenotype_p",
            "joint_statistic", "joint_p", "pathway_statistic", "pathway_p",
            "n_shared_sets", "evidence", "resolution_status",
            "secondary_evidence", "phenotype_column", "phenotype_reference")
  pg <- chorale_programmes(fit, significant_only = TRUE)
  if (nrow(pg) == 0) {
    # An empty result is still a result, so it is written with its columns
    # rather than as a blank file.
    return(stats::setNames(
      as.data.frame(rep(list(character(0)), length(cols)),
                    stringsAsFactors = FALSE),
      cols
    ))
  }

  eff <- function(mod, fac) {
    if (nrow(associations) == 0) return(c(NA_real_, NA_real_))
    r <- associations[associations$modality == mod &
                        associations$factor == fac &
                        associations$covariate == "phenotype", , drop = FALSE]
    if (nrow(r) == 0) return(c(NA_real_, NA_real_))
    c(as.numeric(r$effect_size[1]), as.numeric(r$p_permutation[1]))
  }
  pathway <- function(mod, fac) {
    r <- factors[factors$modality == mod & factors$factor == fac, , drop = FALSE]
    if (nrow(r) == 0 || is.na(r$pathway_definition[1])) return(NA_character_)
    strsplit(r$pathway_definition[1], ";")[[1]][1]
  }
  mk <- function(mod, fac, n = 4L) {
    if (nrow(markers) == 0) return("")
    r <- markers[markers$modality == mod & markers$factor == fac, , drop = FALSE]
    if (nrow(r) == 0) return("")
    paste(chorale_readable_features(utils::head(r$feature, n)), collapse = ", ")
  }

  # The pathway channel is separate but shares fitted factors with the design
  # channel, so every programme says which evidence is available.
  pg <- chorale_evidence_label(pg, fit$pathway_evidence)

  pg$pathway <- vapply(seq_len(nrow(pg)), function(i)
    pathway(pg$modality[i], pg$factor[i]) %||% NA_character_, character(1))
  pg$markers <- vapply(seq_len(nrow(pg)), function(i)
    mk(pg$modality[i], pg$factor[i]), character(1))
  e <- t(vapply(seq_len(nrow(pg)), function(i)
    eff(pg$modality[i], pg$factor[i]), numeric(2)))
  pg$phenotype_effect <- round(e[, 1], 3)
  pg$phenotype_p <- signif(e[, 2], 3)

  # One pathway name per programme: the first its members resolve.
  by_prog <- split(pg, pg$programme)
  pg$programme_pathway <- unlist(lapply(by_prog, function(d) {
    nm <- d$pathway[!is.na(d$pathway)]
    rep(if (length(nm)) nm[1] else "unresolved", nrow(d))
  }), use.names = FALSE)[order(order(pg$programme))]

  pg[, intersect(cols, colnames(pg))]
}

#' Show a feature by its gene symbol where the identifier is opaque
#'
#' A marker reported as `ENSMUSG00000030789` carries no meaning to a reader.
#' Where the identifier resolves to a symbol it is shown as one, with the
#' original retained when it does not.
#'
#' @keywords internal
#' @noRd
chorale_readable_features <- function(x) {
  x <- as.character(x)
  ens <- grepl("^ENSMUS[GT]", x)
  if (!any(ens) || !rlang::is_installed("org.Mm.eg.db")) return(x)
  mapped <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = unique(sub("\\..*$", "", x[ens])),
      keytype = "ENSEMBL", columns = "SYMBOL"
    )),
    error = function(e) NULL
  )
  if (is.null(mapped) || nrow(mapped) == 0) return(x)
  mapped <- mapped[!is.na(mapped$SYMBOL), , drop = FALSE]
  mapped <- mapped[!duplicated(mapped$ENSEMBL), , drop = FALSE]
  key <- stats::setNames(mapped$SYMBOL, mapped$ENSEMBL)
  hit <- unname(key[sub("\\..*$", "", x)])
  unname(ifelse(is.na(hit), x, hit))
}

#' Escape text for HTML
#' @keywords internal
#' @noRd
chorale_esc <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' A table with its legend
#' @keywords internal
#' @noRd
chorale_html_table <- function(d, caption, legend) {
  head_block <- paste0("<h3>", chorale_esc(caption), "</h3>",
                       "<p class='legend'>", legend, "</p>")
  if (is.null(d) || nrow(d) == 0) {
    return(paste0(head_block, "<p class='empty'>No rows.</p>"))
  }
  hdr <- paste0("<tr>", paste0("<th>", chorale_esc(colnames(d)), "</th>",
                               collapse = ""), "</tr>")
  body <- apply(d, 1, function(r) {
    paste0("<tr>", paste0("<td>", chorale_esc(r), "</td>", collapse = ""), "</tr>")
  })
  paste0(head_block, "<div class='scroll'><table>", hdr,
         paste(body, collapse = ""), "</table></div>")
}

#' Dot plot of the phenotype effect per shared programme, per modality
#'
#' The series are the modalities, so colour carries identity and is taken in
#' fixed slot order. Every point is directly labelled with its modality, so
#' identity never rests on colour alone.
#'
#' @keywords internal
#' @noRd
chorale_svg_effects <- function(integrated) {
  if (nrow(integrated) == 0) return("")
  long <- integrated[, c("programme", "modality", "phenotype_effect")]
  names(long)[3] <- "effect"
  long <- long[is.finite(long$effect), , drop = FALSE]
  if (nrow(long) == 0) return("")

  mods <- unique(long$modality)
  slot <- stats::setNames(seq_along(mods), mods)

  progs <- unique(long$programme)
  row_h <- 46
  pad_l <- 92; pad_r <- 200; pad_t <- 18; pad_b <- 46
  w <- 760
  h <- pad_t + length(progs) * row_h + pad_b
  plot_w <- w - pad_l - pad_r

  lim <- max(1e-6, max(abs(long$effect)) * 1.15)
  xs <- function(v) pad_l + (v + lim) / (2 * lim) * plot_w

  parts <- c()
  # zero reference, recessive
  z <- xs(0)
  parts <- c(parts, sprintf(
    "<line x1='%.1f' y1='%d' x2='%.1f' y2='%d' class='axis-strong'/>",
    z, pad_t, z, pad_t + length(progs) * row_h))
  for (tick in c(-lim, -lim / 2, 0, lim / 2, lim)) {
    parts <- c(parts, sprintf(
      "<text x='%.1f' y='%d' class='tick' text-anchor='middle'>%.1f</text>",
      xs(tick), h - pad_b + 26, tick))
  }
  parts <- c(parts, sprintf(
    "<text x='%.1f' y='%d' class='axis-title' text-anchor='middle'>standardised phenotype effect (case minus control)</text>",
    pad_l + plot_w / 2, h - 8))

  for (i in seq_along(progs)) {
    pr <- progs[i]
    y <- pad_t + (i - 0.5) * row_h
    sub <- long[long$programme == pr, , drop = FALSE]
    parts <- c(parts, sprintf(
      "<text x='%d' y='%.1f' class='rowlab' text-anchor='end'>%s</text>",
      pad_l - 14, y + 4, chorale_esc(pr)))
    # The two modalities are dodged vertically, so a programme whose effects
    # nearly coincide still shows both marks rather than one hiding the other.
    n_sub <- nrow(sub)
    dodge <- if (n_sub > 1) seq(-5.5, 5.5, length.out = n_sub) else 0
    if (n_sub > 1 && all(is.finite(sub$effect))) {
      ordx <- order(sub$effect)
      parts <- c(parts, sprintf(
        "<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' class='connect'/>",
        xs(sub$effect[ordx[1]]), y + dodge[1],
        xs(sub$effect[ordx[n_sub]]), y + dodge[n_sub]))
    }
    for (j in seq_len(nrow(sub))) {
      cx <- xs(sub$effect[j]); s <- slot[[sub$modality[j]]]
      parts <- c(parts, sprintf(
        "<circle cx='%.1f' cy='%.1f' r='6' class='pt s%d'><title>%s, %s: %.3f</title></circle>",
        cx, y + dodge[j], s, chorale_esc(pr), chorale_esc(sub$modality[j]),
        sub$effect[j]))
    }
    lab <- paste(sprintf("%s %.2f", sub$modality, sub$effect), collapse = "  ")
    parts <- c(parts, sprintf(
      "<text x='%d' y='%.1f' class='ptlab'>%s</text>",
      w - pad_r + 10, y + 4, chorale_esc(lab)))
  }

  legend <- paste0(
    "<div class='legend-row'>",
    paste(vapply(mods, function(m) sprintf(
      "<span class='key'><i class='sw s%d'></i>%s</span>", slot[[m]],
      chorale_esc(m)), character(1)), collapse = ""),
    "</div>")

  paste0(legend, "<svg viewBox='0 0 ", w, " ", h,
         "' width='100%' role='img' aria-label='Phenotype effect per shared programme in each modality'>",
         paste(parts, collapse = ""), "</svg>")
}

#' Dumbbell of the identified set, without and with anchors
#'
#' The same quantity under two conditionings, so one hue in two steps rather
#' than two identities.
#'
#' @keywords internal
#' @noRd
chorale_svg_bounds <- function(bounds, integrated = NULL) {
  if (is.null(bounds) || nrow(bounds) == 0) return("")
  # Rows carry the programme identifier used in the table above, so the two
  # figures can be read against one another.
  bounds$label <- paste0(bounds$modality_a, "/", bounds$modality_b)
  if (!is.null(integrated) && nrow(integrated) > 0) {
    member <- stats::setNames(integrated$programme,
                              paste(integrated$modality, integrated$factor))
    hit <- member[paste(bounds$modality_a, bounds$factor_a)]
    hit2 <- member[paste(bounds$modality_b, bounds$factor_b)]
    prog <- ifelse(is.na(hit), hit2, hit)
    bounds$label <- ifelse(is.na(prog), bounds$label,
                           paste0(prog, "  ", bounds$label))
  }
  bounds <- utils::head(bounds[order(bounds$width_anchored), , drop = FALSE], 12)
  row_h <- 34
  pad_l <- 158; pad_r <- 40; pad_t <- 14; pad_b <- 46
  w <- 760
  h <- pad_t + nrow(bounds) * row_h + pad_b
  plot_w <- w - pad_l - pad_r
  xs <- function(v) pad_l + (v + 1) / 2 * plot_w

  parts <- c(sprintf(
    "<line x1='%.1f' y1='%d' x2='%.1f' y2='%d' class='axis-strong'/>",
    xs(0), pad_t, xs(0), pad_t + nrow(bounds) * row_h))
  for (tick in c(-1, -0.5, 0, 0.5, 1)) {
    parts <- c(parts, sprintf(
      "<text x='%.1f' y='%d' class='tick' text-anchor='middle'>%.1f</text>",
      xs(tick), h - pad_b + 26, tick))
  }
  parts <- c(parts, sprintf(
    "<text x='%.1f' y='%d' class='axis-title' text-anchor='middle'>correlation the data cannot exclude</text>",
    pad_l + plot_w / 2, h - 8))

  for (i in seq_len(nrow(bounds))) {
    b <- bounds[i, ]
    y <- pad_t + (i - 0.5) * row_h
    lab <- b$label
    parts <- c(parts, sprintf(
      "<text x='%d' y='%.1f' class='rowlab' text-anchor='end'>%s</text>",
      pad_l - 12, y + 4, chorale_esc(lab)))
    parts <- c(parts, sprintf(
      "<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' class='band-wide'><title>%s without anchors: %.2f to %.2f</title></line>",
      xs(b$lower_no_anchor), y, xs(b$upper_no_anchor), y,
      chorale_esc(lab), b$lower_no_anchor, b$upper_no_anchor))
    parts <- c(parts, sprintf(
      "<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' class='band-narrow'><title>%s with anchors: %.2f to %.2f</title></line>",
      xs(b$lower_anchored), y, xs(b$upper_anchored), y,
      chorale_esc(lab), b$lower_anchored, b$upper_anchored))
  }

  legend <- paste0(
    "<div class='legend-row'>",
    "<span class='key'><i class='sw wide'></i>without anchors</span>",
    "<span class='key'><i class='sw narrow'></i>with anchors</span></div>")

  paste0(legend, "<svg viewBox='0 0 ", w, " ", h,
         "' width='100%' role='img' aria-label='Identified set per coupling, without and with anchors'>",
         paste(parts, collapse = ""), "</svg>")
}

#' Pathway corroboration per programme, ready to render
#' @keywords internal
#' @noRd
chorale_pathway_table <- function(fit) {
  pe <- fit$pathway_evidence
  if (is.null(pe) || nrow(pe) == 0) return(data.frame())
  pg <- chorale_programmes(fit, significant_only = TRUE)
  if (nrow(pg) == 0) return(data.frame())
  pe <- pe[pe$programme %in% pg$programme, , drop = FALSE]
  if (nrow(pe) == 0) return(data.frame())
  span <- pg[!duplicated(pg$programme), c("programme", "modalities")]
  merge(span, pe, by = "programme", sort = FALSE)
}

#' Self-contained report
#' @keywords internal
#' @noRd
chorale_write_html <- function(fit, factors, markers, associations,
                               concordance, bounds, controls, path,
                               derived = NULL) {
  if (is.null(derived)) {
    derived <- list(
      leave_one_out = chorale_leave_one_out(fit),
      pathway = chorale_pathway_table(fit),
      fdr = chorale_fdr(fit),
      added_value = chorale_added_value(fit),
      cohorts = chorale_cohort_overlap(fit)
    )
  }
  f <- file.path(path, "report.html")
  integrated <- chorale_integrated_table(fit, factors, markers, associations,
                                         concordance)
  n_prog <- if (nrow(integrated) > 0) length(unique(integrated$programme)) else 0
  n_multi <- if (nrow(integrated) > 0) {
    length(unique(integrated$programme[integrated$n_modalities >= 3]))
  } else 0

  css <- "
:root{--surface-1:#fcfcfb;--surface-2:#f4f3f0;--text-primary:#0b0b0b;--text-secondary:#52514e;--text-muted:#767570;--rule:#dcdbd6;--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--seq-light:#86b6ef;--seq-dark:#1c5cab;--good:#0ca30c;--warn:#fab219}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme=light])){--surface-1:#1a1a19;--surface-2:#242423;--text-primary:#fff;--text-secondary:#c3c2b7;--text-muted:#9a998f;--rule:#3a3a37;--s1:#3987e5;--s2:#d95926;--s3:#199e70;--seq-light:#3987e5;--seq-dark:#86b6ef}}
:root[data-theme=dark]{--surface-1:#1a1a19;--surface-2:#242423;--text-primary:#fff;--text-secondary:#c3c2b7;--text-muted:#9a998f;--rule:#3a3a37;--s1:#3987e5;--s2:#d95926;--s3:#199e70;--seq-light:#3987e5;--seq-dark:#86b6ef}
*{box-sizing:border-box}
body{font-family:system-ui,-apple-system,sans-serif;margin:0;padding:2rem 1.25rem 4rem;background:var(--surface-1);color:var(--text-primary);line-height:1.55}
main{max-width:informalmax;max-width:70rem;margin:0 auto}
h1{font-size:1.7rem;margin:0 0 .25rem}
h2{font-size:1.25rem;margin:2.5rem 0 .5rem;padding-bottom:.3rem;border-bottom:1px solid var(--rule)}
h3{font-size:1rem;margin:1.5rem 0 .3rem}
p.lede{color:var(--text-secondary);margin:.2rem 0 1.2rem}
p.legend{color:var(--text-secondary);font-size:.86rem;margin:.1rem 0 .6rem;max-width:60rem}
p.empty{color:var(--text-muted);font-style:italic}
.hero{display:flex;gap:2.5rem;flex-wrap:wrap;background:var(--surface-2);border:1px solid var(--rule);border-radius:10px;padding:1.1rem 1.4rem;margin:1rem 0 .4rem}
.hero .n{font-size:2.6rem;font-weight:650;line-height:1.1}
.hero .cap{color:var(--text-secondary);font-size:.85rem;max-width:16rem}
.scroll{overflow-x:auto;-webkit-overflow-scrolling:touch}
table{border-collapse:collapse;font-size:.82rem;min-width:100%}
th,td{border-bottom:1px solid var(--rule);padding:.36rem .6rem;text-align:left;white-space:nowrap}
th{background:var(--surface-2);font-weight:600;position:sticky;top:0}
.legend-row{display:flex;gap:1.2rem;flex-wrap:wrap;margin:.5rem 0 .2rem;font-size:.85rem;color:var(--text-secondary)}
.key{display:inline-flex;align-items:center;gap:.4rem}
.sw{width:14px;height:14px;border-radius:3px;display:inline-block}
.sw.s1{background:var(--s1)}.sw.s2{background:var(--s2)}.sw.s3{background:var(--s3)}
.sw.wide{background:var(--seq-light)}.sw.narrow{background:var(--seq-dark)}
svg{display:block;margin:.2rem 0 .5rem;overflow:visible}
.pt{stroke:var(--surface-1);stroke-width:2}
.pt.s1{fill:var(--s1)}.pt.s2{fill:var(--s2)}.pt.s3{fill:var(--s3)}
.connect{stroke:var(--rule);stroke-width:2}
.axis-strong{stroke:var(--text-muted);stroke-width:1.5}
.tick,.ptlab{font-size:11px;fill:var(--text-secondary)}
.rowlab{font-size:12px;fill:var(--text-primary);font-weight:550}
.axis-title{font-size:11.5px;fill:var(--text-secondary)}
.band-wide{stroke:var(--seq-light);stroke-width:9;stroke-linecap:round}
.band-narrow{stroke:var(--seq-dark);stroke-width:9;stroke-linecap:round}
.card{border:1px solid var(--rule);border-radius:8px;padding:.7rem .95rem;margin:.5rem 0;background:var(--surface-2)}
.card h4{margin:0 0 .3rem;font-size:.95rem}
details.how{border:1px solid var(--rule);border-radius:8px;padding:.6rem .95rem;margin:.6rem 0 0;background:var(--surface-2)}
details.how summary{cursor:pointer;font-weight:600;font-size:.92rem}
details.how p{color:var(--text-secondary);font-size:.88rem;max-width:60rem}
.tag{font-size:.72rem;padding:.1rem .45rem;border-radius:99px;border:1px solid var(--rule);color:var(--text-secondary);margin-left:.4rem}
.tag.ok{color:var(--good);border-color:var(--good)}
.tag.warn{color:var(--warn);border-color:var(--warn)}
"

  cards <- ""
  if (nrow(integrated) > 0) {
    for (pr in unique(integrated$programme)) {
      d <- integrated[integrated$programme == pr, , drop = FALSE]
      members <- paste0(vapply(seq_len(nrow(d)), function(i) paste0(
        "<strong>", chorale_esc(d$modality[i]), "</strong> ",
        chorale_esc(d$factor[i]),
        ", markers ", chorale_esc(d$markers[i]),
        " &middot; phenotype effect ", chorale_esc(d$phenotype_effect[i]),
        " (p ", chorale_esc(d$phenotype_p[i]), ")"), character(1)),
        collapse = "<br>")
      cards <- paste0(cards,
        "<div class='card'><h4>", chorale_esc(pr), " &middot; ",
        chorale_esc(d$programme_pathway[1]),
        "<span class='tag ", if (d$n_modalities[1] >= 3) "ok" else "",
        "'>", d$n_modalities[1], " modalities</span>",
        if ("joint_p" %in% colnames(d) && !is.na(d$joint_p[1]))
          paste0("<span class='tag'>joint p ", signif(d$joint_p[1], 3), "</span>")
        else "",
        "<span class='tag'>", chorale_esc(d$modalities[1]), "</span></h4>",
        "<p class='legend'>", members, "</p></div>")
    }
  }

  headline <- if (nrow(integrated) == 0) {
    paste0(
      "<p class='lede'>No programme survived its null in this run, so there is no integrated ",
      "result to report. The assignment step always returns a correspondence between factors, ",
      "because it is an assignment; what it does not do is establish that the correspondence ",
      "means anything. Here it did not, and the supporting detail below describes each modality ",
      "on its own.</p>")
  } else {
    paste0("<div class='hero'>",
      "<div><div class='n'>", n_prog, "</div>",
      "<div class='cap'>latent programmes seen in more than one modality</div></div>",
      "<div><div class='n'>", n_multi, "</div>",
      "<div class='cap'>of them carried by all three modalities at once</div></div>",
      "<div><div class='n'>", length(fit$modalities), "</div>",
      "<div class='cap'>modalities integrated, on disjoint animals</div></div>",
      "</div>")
  }

  integrated_display <- integrated

  html <- paste0(
    "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
    "<meta name='viewport' content='width=device-width,initial-scale=1'>",
    "<title>CHORALE integrated report</title><style>", css, "</style></head><body><main>",

    "<h1>What integration found</h1>",
    "<p class='lede'>Each modality was measured on different animals, so nothing here comes from ",
    "matching individuals. A programme below is a latent axis recovered separately in two ",
    "modalities whose behaviour across the case/control design agrees, which is what makes it ",
    "visible only on integration.</p>",
    headline,
    "<details class='how'><summary>How two datasets with no animal in common are compared</summary>",
    "<p>The modalities share no animal, so nothing can be matched individual to individual. What ",
    "they do share is the experimental design: every cohort has cases and controls, and usually ",
    "ages and sexes too.</p>",
    "<p>A factor is a latent axis, estimated inside one modality from its own animals. Of each ",
    "factor one question is asked, and it can be answered without leaving that modality: ",
    "<em>how much higher is this factor in cases than in controls?</em> The answer is one number ",
    "per factor per modality. Age and sex, where recorded, add further numbers of the same kind, ",
    "and together they form that factor's design profile.</p>",
    "<p>Two factors in two modalities are candidates for reflecting the same biology when their ",
    "design profiles agree. That comparison is what <strong>anchoring on the design</strong> ",
    "means: the design is the common reference frame, doing the job matched individuals would do ",
    "if they existed. Agreement alone proves nothing, since unrelated factors agree by chance, so ",
    "the case/control labels are shuffled among animals many times and the observed agreement is ",
    "judged against what shuffling produces.</p>",
    "<p>A programme is not a pair, and it is not assembled from pairs. All the modalities are ",
    "compared at once: every pairwise agreement is placed in one matrix, and the assignment of ",
    "factors to programmes is solved over the whole collection in a single step. Correspondences ",
    "found this way agree around every cycle, so a programme cannot be created by chaining one ",
    "modality to the next through a single mistaken link, and how many modalities a programme ",
    "spans is decided by the same evidence rather than afterwards.</p>",
    "<p>Knowing two factors behave alike across the design still does not say how they are coupled ",
    "within one animal, and no amount of data of this kind will, because each animal was measured ",
    "once. What the data support is a range of couplings, which the second figure reports.</p>",
    "</details>",

    "<h2>Shared biological programmes</h2>",
    "<p class='legend'>Positive effect means higher in cases than in controls, in ",
    "units of the pooled standard deviation. Points joined by a line are the same programme seen ",
    "in two modalities. Matching tests whether the design profiles point the same way, in either ",
    "direction, so a programme whose modalities move oppositely is matched and its opposition ",
    "reported. Hover any point for its value.</p>",
    chorale_svg_effects(integrated),
    chorale_html_table(integrated_display, "Programmes, their biology and their markers",
      paste0("One row per modality that carries a programme, so a programme measured in three ",
             "modalities occupies three rows sharing a <em>programme</em> identifier. ",
             "<em>n_modalities</em> and <em>modalities</em> record how far it reaches. ",
             "<em>programme_pathway</em> is the curated set that best describes the ",
             "factor's loadings. <em>markers</em> are features loading on that programme and ",
             "almost nothing else, which is what makes the axis interpretable. ",
             "<em>phenotype_effect</em> and <em>phenotype_p</em> are the case/control effect and ",
             "its permutation p-value within that modality; comparing the sign across rows of one ",
             "programme shows whether the modalities move together or oppositely, and opposed ",
             "programmes are kept, since transcript and protein disagreement is itself reported. ",
             "<em>joint_statistic</em> and <em>joint_p</em> judge the programme as one object: ",
             "the agreement is averaged over every pair inside it, and the null permutes ",
             "covariate labels in every contributing modality at once, keeping the best value any ",
             "combination of factors could reach. Requiring three modalities to agree at once is ",
             "far harder by chance than requiring two, so this is where measuring more modalities ",
             "pays. The null reruns the whole procedure, assignment included, so the p-value ",
             "accounts for the freedom the estimator had rather than conditioning on the ",
             "programme it chose. Only programmes whose joint evidence beat that null appear ",
             "here. <em>evidence</em> says whether the biology corroborates ",
             "the design, which the section below sets out.")),

    "<h2>Does the biology agree as well as the design?</h2>",
    "<p class='legend'>A programme is recovered by asking whether factors in ",
    "different modalities respond to the case/control design the same way. ",
    "Whether they also implicate the same biology is a separate question, and ",
    "it is asked separately here. The factors were fitted without the curated ",
    "sets, so agreement on biology cannot have been built into them, and a ",
    "programme carrying both kinds of evidence has corroboration from a ",
    "separate channel. The channels are not statistically independent because ",
    "both use the fitted factors. The null permutes score residuals relative ",
    "to the unchanged assay and recalculates loadings, preserving feature ",
    "correlation and annotations.</p>",
    chorale_html_table(derived$pathway,
      "Biological corroboration, by programme",
      paste0("<em>pathway_statistic</em> is the agreement of the members' ",
             "pathway profiles, averaged over every pair inside the programme, ",
             "and <em>pathway_p</em> calibrates it against the ",
             "score-residual null. <em>n_shared_sets</em> is how much ",
             "vocabulary the modalities have in common: where it is small the ",
             "channel has little to say, and an empty result means the ",
             "question could not be asked rather than that the answer was no. ",
             "<em>evidence</em> records which channels the programme rests ",
             "on.")),

    "<h2>Multiplicity and added value</h2>",
    "<p class='legend'>A run tests every factor, every programme and every ",
    "programme's biology, so the surviving results carry a false-discovery rate ",
    "controlled within each of those levels. Beside it, each programme is set ",
    "against the strongest single modality among its members: a programme that ",
    "does not beat the best single modality was visible in one layer and ",
    "accompanied by the others, whatever its joint p-value.</p>",
    chorale_html_table(derived$fdr, "False discovery across the search",
      paste0("<em>level</em> is what was tested: a factor against the design, a ",
             "programme against its joint null, or a programme's biology. ",
             "<em>q_value</em> is the Benjamini-Hochberg rate within that ",
             "level; levels are corrected separately, since pooling them would ",
             "penalise a programme for how many factors happened to be fitted.")),
    chorale_html_table(derived$added_value,
      "Whether a programme needed more than one modality",
      paste0("<em>best_single_p</em> is the strongest phenotype association any ",
             "one member factor reaches alone, and <em>margin</em> is the gap ",
             "from the joint evidence on the negative log scale. ",
             "<em>worst_leave_one_out_delta</em> is the largest change in the ",
             "statistic when a modality is removed, read from the statistic ",
             "rather than the p-value because a permutation p-value saturates ",
             "at its floor. <em>needs_multiple</em> requires both: the ",
             "programme beats every single modality, and removing any modality ",
             "weakens it.")),
    chorale_html_table(derived$cohorts,
      "The population the cohorts share",
      paste0("A programme describes the population every contributing cohort ",
             "represents. <em>common_levels</em> are the levels all cohorts ",
             "populate and <em>min_share_covered</em> the smallest fraction of ",
             "a cohort those levels hold. <em>max_total_variation</em> is the ",
             "largest distance between two cohorts' distributions of the ",
             "covariate, zero meaning they realise it identically. A covariate ",
             "marked not <em>comparable</em> is one no claim should cross.")),

    "<h2>What each modality contributes</h2>",
    "<p class='legend'>A programme carried by three modalities is only a result of integration ",
    "if it needs all three. Each row removes one modality and rescores the programme without ",
    "it: a small <em>delta</em> means the programme was resting on the modalities that remain, ",
    "and a large negative one means the removed modality was carrying the evidence.</p>",
    chorale_html_table(derived$leave_one_out,
      "Joint evidence with one modality removed",
      paste0("<em>dropped</em> is the modality left out, <em>joint_statistic</em> and ",
             "<em>joint_p</em> the evidence for the remainder, and <em>delta</em> the change ",
             "from the full programme. Empty where no programme spans three modalities.")),

    "<h2>How far the data constrain the coupling</h2>",
    "<p class='legend'>Each animal was measured once, so how strongly two programmes move ",
    "together within an animal cannot be recovered. What the data support is a range. Each row ",
    "is the range of correlations they cannot exclude. Taking the design into account narrows ",
    "it: the average of a programme within each design group is observable in every modality ",
    "carrying it, so the part of the relationship running between groups is already fixed, and ",
    "only the variation inside a group stays free. Where a bar remains wide, that width is the ",
    "result rather than a missing number. A correlation concerns two quantities, so a programme ",
    "spanning three modalities contributes one row per pair, all describing that one programme.</p>",
    chorale_svg_bounds(bounds, integrated),
    chorale_html_table(bounds, "Identified set per coupling",
      paste0("<em>lower/upper_no_anchor</em> is the range compatible with the marginal ",
             "scores alone; <em>lower/upper_anchored</em> also conditions on the design. ",
             "For example, [0.35, 0.70] allows any value in that positive range, while ",
             "[-0.95, 0.96] leaves direction undetermined. Width reflects missing pairing; ",
             "a bootstrap sensitivity envelope reflects resampling variation.")),

    "<h2>Controls</h2>",
    "<p class='legend'>Every result above is shown with the check that would have caught it if it ",
    "were an artefact rather than biology. A check that cannot be computed on these data says so ",
    "rather than passing quietly.</p>",
    chorale_html_table(controls, "Controls accompanying this run",
      paste0("<em>phenotype permutation</em> refits after shuffling the case/control label within ",
             "stratum; a small value means the observed agreement exceeds what shuffled labels ",
             "produce. <em>modality shuffle</em> reassigns samples across modalities many times and ",
             "reports the p-value of the observed evidence against that null; it is defined only ",
             "where the modalities share a feature space. <em>factor stability</em> is the mean ",
             "agreement of the recovered factors across restarts, matched one to one; independent ",
             "component analysis is non-convex, so a low value means the recovery is a draw rather ",
             "than an estimate.")),
    chorale_html_table(concordance, "Cross-modality concordance",
      paste0("<em>agreement</em> is the rank correlation of the two modalities' stratum profiles. ",
             "<em>discordant_strata</em> counts design strata where the two disagree in sign.")),

    "<h2>Supporting detail</h2>",
    "<p class='legend'>The tables below describe each modality on its own. They support the ",
    "cross-modal correspondences above rather than establishing one by themselves.</p>",
    chorale_html_table(factors, "Every factor recovered, by modality",
      paste0("<em>shared</em> marks factors entering a programme above. <em>n_markers</em> counts ",
             "features clearing the purity threshold and <em>purity_margin</em> reports how ",
             "cleanly the best ones separate, so a factor with weak markers is visible as such. ",
             "<em>status</em> records why a factor is unresolved: a factor that cannot be named ",
             "is reported, never dropped.")),
    chorale_html_table(utils::head(associations, 40),
      "Design associations, strongest 40",
      paste0("Association of each factor with each design covariate, within modality. ",
             "<em>effect_size</em> is a standardised mean difference for a two-level covariate ",
             "and a rank correlation otherwise. <em>p_adjusted</em> controls the false discovery ",
             "rate across all rows.")),
    chorale_html_table(utils::head(markers, 60), "Marker features, first 60",
      paste0("Features loading on one factor and negligibly on the others. ",
             "<em>clears_purity_threshold</em> distinguishes a marker meeting the condition from ",
             "the nearest available feature where none does.")),

    "<p class='legend' style='margin-top:2rem'>Per-sample factor scores and full feature loadings ",
    "are written alongside this report as <code>scores_&lt;modality&gt;.tsv</code> and ",
    "<code>loadings_&lt;modality&gt;.tsv</code>. The complete model is in ",
    "<code>chorale_mae.rds</code>; factor definitions in GMT form are in <code>factors.gmt</code>.</p>",
    "</main></body></html>"
  )
  writeLines(html, f)
  f
}
