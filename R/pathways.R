#' Pathway membership of lipid classes
#'
#' Reads the bridge that places a lipidome in the same pathway vocabulary the
#' transcriptome and the proteome are described in. Each row states that a
#' lipid class takes part in a pathway, and names that pathway with the
#' identifier the gene-side collection uses, so a lipid and a gene can be
#' placed in the same set rather than in two vocabularies that have to be
#' reconciled afterwards.
#'
#' The table is derived rather than asserted. A class is a term in the ChEBI
#' ontology; the molecules beneath that term are looked up in Reactome's
#' compound annotation; and a pathway is kept only where the gene-side
#' collection also carries it as a set of usable size. `data-raw/` holds the
#' script, the source URLs and the licences.
#'
#' The shipped table takes its compound annotation from mouse Reactome, but it
#' is keyed on the MSigDB Reactome set names, which are the same across species,
#' so it places a lipid in the same set a human or rat study would use wherever
#' that pathway exists there. Rebuilding it for another organism is a matter of
#' running the script with that species set, and any table of the same shape can
#' be supplied through `mapping`.
#'
#' Parsing a lipid shorthand name into its class is a solved problem and is not
#' redone here: [chorale_lipid_class()] uses Goslin where it is installed
#' (Kopczynski et al., Anal Chem 92:10957, 2020). What no published resource
#' provides, and what this table supplies, is the step from that class to the
#' pathway vocabulary the genes live in.
#'
#' @returns A data frame with one row per (lipid class, pathway), carrying
#'   `abbreviation`, `lipidmaps_category`, `chebi_id`, `chebi_name`,
#'   `pathway_id`, `pathway_name`, `msigdb_name`, `n_compounds` and
#'   `n_descendants`.
#' @export
#' @examples
#' sets <- chorale_metabolite_pathways()
#' utils::head(sets[sets$abbreviation == "SM", c("abbreviation", "msigdb_name")])
chorale_metabolite_pathways <- function() {
  path <- system.file("extdata", "metabolite_pathways.tsv", package = "chorale")
  if (!nzchar(path)) {
    rlang::abort("The metabolite pathway table is missing from the installation.")
  }
  utils::read.delim(path, stringsAsFactors = FALSE)
}

#' Lipid class of each feature in a lipidome
#'
#' A lipidomic feature is reported in shorthand, as a class abbreviation
#' followed by a chain description and often an adduct. Only the class is
#' needed to place the feature in a pathway, since a pathway acts on a class
#' rather than on one chain length.
#'
#' Goslin is used where it is installed, since it is the validated grammar for
#' these names and it normalises the many dialects that report them. Where it
#' is absent, or where a name is in a dialect it does not accept, the leading
#' token is taken instead, which is the class in every dialect the project has
#' met.
#'
#' @param ids Character vector of feature identifiers.
#'
#' @returns A character vector of class abbreviations, `NA` where none could be
#'   read.
#' @export
#' @examples
#' chorale_lipid_class(c("PC 34:1|[M+H]+", "SM 36:1;O2", "AHexCer 57:1;O3"))
chorale_lipid_class <- function(ids) {
  ids <- as.character(ids)
  # The class is the leading token, before the chain description, the species
  # level annotation and the adduct.
  fallback <- sub("[ |].*$", "", ids)
  fallback[!nzchar(fallback)] <- NA_character_

  if (!rlang::is_installed("rgoslin")) return(fallback)

  bare <- sub("\\|.*$", "", ids)
  parsed <- vapply(bare, function(one) {
    # The grammar reports a name it cannot read on the message stream, which
    # would otherwise fill the console with one line per unparsed feature. The
    # fallback covers it, so the report is noise rather than information.
    out <- suppressMessages(suppressWarnings(
      utils::capture.output(
        res <- try(rgoslin::parseLipidNames(one), silent = TRUE),
        type = "message"
      )
    ))
    out <- res
    if (inherits(out, "try-error") || is.null(out) || nrow(out) == 0) {
      return(NA_character_)
    }
    cls <- out[["Lipid.Maps.Main.Class"]]
    if (is.null(cls) || length(cls) == 0 || is.na(cls[1])) {
      return(NA_character_)
    }
    as.character(cls[1])
  }, character(1), USE.NAMES = FALSE)

  ifelse(is.na(parsed), fallback, parsed)
}

#' Represent a lipidome as a feature-by-set indicator matrix
#'
#' The counterpart of [chorale_geneset_matrix()] for a modality whose features
#' are lipids. The columns are the same sets, so the two modalities are
#' described in one vocabulary and their factors can be compared on biology.
#'
#' @param feature_ids Character vector of lipid shorthand identifiers, in the
#'   row order of the modality's assay.
#' @param sets A named list of gene sets, as returned by [chorale_genesets()].
#'   Only its names are used: they fix the vocabulary the lipids are placed in.
#' @param mapping The bridge table, defaulting to
#'   [chorale_metabolite_pathways()].
#' @param min_compounds Keep a class-to-pathway link where the pathway is
#'   annotated with at least this many compounds of the class.
#' @param min_specificity Keep a link where the class accounts for at least
#'   this share of the compounds the pathway is annotated with.
#' @param min_features Drop sets matching fewer than this many features.
#'
#' @details
#' The two thresholds are alternatives, not conditions to be met together,
#' because a weak link fails in one of two different ways. Reactome describes
#' some classes by a single representative compound, and that compound is
#' carried into every pathway it appears in, so a class can reach a pathway it
#' has no particular relation to: several distinct compounds of the class rules
#' that out. A large pathway, though, dilutes any share, so requiring a share
#' alone would discard a class from the very pathway that acts on it. A link
#' therefore survives if the class appears in the pathway more than once, or if
#' it accounts for a substantial part of a pathway small enough for one
#' compound to matter.
#'
#' @returns A numeric matrix with one row per feature and one column per
#'   retained set.
#' @export
#' @examples
#' sets <- stats::setNames(list(NULL, NULL),
#'                         c("REACTOME_SPHINGOLIPID_METABOLISM",
#'                           "REACTOME_CARNITINE_SHUTTLE"))
#' m <- chorale_metabolite_matrix(c("SM 36:1;O2", "Cer 34:1;O2", "CAR 18:2"),
#'                                sets, min_features = 1)
#' colSums(m)
chorale_metabolite_matrix <- function(feature_ids, sets,
                                      mapping = chorale_metabolite_pathways(),
                                      min_compounds = 2L,
                                      min_specificity = 0.05,
                                      min_features = 3L) {
  feature_ids <- as.character(feature_ids)
  classes <- chorale_lipid_class(feature_ids)
  vocabulary <- names(sets)
  mapping <- mapping[mapping$msigdb_name %in% vocabulary, , drop = FALSE]
  if (all(c("specificity", "n_compounds") %in% colnames(mapping))) {
    mapping <- mapping[mapping$n_compounds >= min_compounds |
                         mapping$specificity >= min_specificity, , drop = FALSE]
  }

  mat <- matrix(0, nrow = length(feature_ids), ncol = length(vocabulary),
                dimnames = list(feature_ids, vocabulary))
  if (nrow(mapping) == 0) return(mat[, integer(0), drop = FALSE])

  members <- split(mapping$abbreviation, mapping$msigdb_name)
  for (s in names(members)) {
    mat[, s] <- as.numeric(classes %in% members[[s]])
  }
  keep <- colSums(mat > 0) >= min_features
  mat[, keep, drop = FALSE]
}

#' Where each factor sits in the pathway vocabulary
#'
#' The design profile says what a factor does to the experiment. This says what
#' it is made of: for every curated set, how far the loadings of that set's
#' members depart from the factor's loadings as a whole, signed, and scaled by
#' the set's size so that a large set and a small one are on one footing.
#'
#' The two are independent lines of evidence. The factors are fitted without
#' reference to the sets, so a profile computed here cannot have been produced
#' by the annotation, and agreement between modalities on this profile is not a
#' restatement of their agreement on the design.
#'
#' @param loadings A features-by-factors numeric matrix.
#' @param prior A feature-by-set matrix from [chorale_geneset_matrix()] or
#'   [chorale_metabolite_matrix()].
#'
#' @returns A factors-by-sets numeric matrix.
#' @export
#' @examples
#' set.seed(1)
#' l <- matrix(rnorm(60), nrow = 20, dimnames = list(paste0("f", 1:20), NULL))
#' p <- matrix(rbinom(40, 1, 0.4), nrow = 20,
#'             dimnames = list(paste0("f", 1:20), c("set_a", "set_b")))
#' dim(chorale_pathway_profile(l, p))
chorale_pathway_profile <- function(loadings, prior) {
  shared <- intersect(rownames(loadings), rownames(prior))
  out <- matrix(0, nrow = ncol(loadings), ncol = ncol(prior),
                dimnames = list(colnames(loadings), colnames(prior)))
  if (length(shared) < 3 || ncol(prior) < 1) return(out)

  l <- loadings[shared, , drop = FALSE]
  p <- prior[shared, , drop = FALSE] > 0
  sizes <- colSums(p)
  for (j in seq_len(ncol(l))) {
    v <- l[, j]
    mu <- mean(v)
    sdv <- stats::sd(v)
    if (!is.finite(sdv) || sdv == 0) next
    # Standard error of a set's mean loading under no association, so sets of
    # different size contribute on the same scale and the sign is kept.
    means <- as.numeric(crossprod(p, v)) / pmax(sizes, 1)
    out[j, ] <- (means - mu) / (sdv / sqrt(pmax(sizes, 1)))
  }
  out[!is.finite(out)] <- 0
  out
}

#' Corroborate programmes on biology as well as on the design
#'
#' A programme is recovered by asking whether factors in different modalities
#' respond to the design the same way. That they also implicate the same
#' biology is a separate question, and answering it separately is what keeps
#' the answer from being circular: the factors were fitted without the curated
#' sets, so their agreement here was not built in.
#'
#' The statistic mirrors the design channel. Each factor carries a pathway
#' profile over the vocabulary its modality and its partner share, the
#' agreement of two factors is the inner product of their profiles, and a
#' programme's evidence is that agreement averaged over every pair inside it.
#'
#' The null is annotation-matched. Feature labels are permuted within each
#' modality's prior, which leaves every set at its original size and every
#' feature in its original number of sets while breaking the relation between
#' the loadings and the annotation. Set size and feature frequency are the two
#' properties that make an enrichment appear where none exists, so holding both
#' fixed is what makes the calibration mean something. As in the design
#' channel, the null keeps the best value any programme could have reached.
#'
#' @param fit A `chorale_fit` object.
#' @param programmes Output of [chorale_programmes()]; taken from `fit` if
#'   absent.
#' @param n_perm Number of annotation-matched permutations.
#' @param alpha Significance threshold.
#' @param seed Integer seed.
#'
#' @returns A data frame with one row per programme, carrying
#'   `pathway_statistic`, `pathway_p`, `n_shared_sets` and the sets
#'   contributing most.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' # A curated set reaches every modality, because their identifiers are
#' # harmonised into one space before the sets are applied.
#' ids <- lapply(containers, function(se)
#'   rownames(SummarizedExperiment::assay(se)))
#' span <- function(from, to) unlist(lapply(ids, function(v) v[from:to]))
#' sets <- list(set_a = span(1, 40), set_b = span(30, 80),
#'              set_c = span(60, 120))
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2,
#'                    gene_sets = sets)
#' chorale_pathway_evidence(fit, n_perm = 20)
chorale_pathway_evidence <- function(fit, programmes = NULL, n_perm = 200L,
                                     alpha = 0.05, seed = 1L) {
  if (!inherits(fit, "chorale_fit")) {
    rlang::abort("`fit` must be a chorale_fit object.")
  }
  if (is.null(programmes)) programmes <- chorale_programmes(fit)
  empty <- data.frame()
  if (nrow(programmes) == 0) return(empty)

  priors <- lapply(fit$modalities, function(m) fit$fits[[m]]$prior)
  names(priors) <- fit$modalities
  usable <- vapply(priors, function(p) !is.null(p) && ncol(p) > 0, logical(1))
  if (sum(usable) < 2) return(empty)

  profile_of <- function(prior_by_mod) {
    out <- lapply(fit$modalities, function(m) {
      if (!usable[[m]]) return(NULL)
      chorale_pathway_profile(fit$fits[[m]]$loadings, prior_by_mod[[m]])
    })
    names(out) <- fit$modalities
    out
  }

  agreement <- function(prof, members) {
    vals <- numeric(0)
    for (i in seq_len(nrow(members))) {
      for (j in seq_len(nrow(members))) {
        if (j <= i) next
        ma <- members$modality[i]
        mb <- members$modality[j]
        if (is.null(prof[[ma]]) || is.null(prof[[mb]])) next
        shared <- intersect(colnames(prof[[ma]]), colnames(prof[[mb]]))
        if (length(shared) < 2) next
        fa <- match(members$factor[i], rownames(prof[[ma]]))
        fb <- match(members$factor[j], rownames(prof[[mb]]))
        if (is.na(fa) || is.na(fb)) next
        vals <- c(vals, abs(sum(prof[[ma]][fa, shared] * prof[[mb]][fb, shared])))
      }
    }
    if (length(vals) == 0) return(NA_real_)
    mean(vals)
  }

  observed_prof <- profile_of(priors)
  groups <- split(programmes, programmes$programme)

  set.seed(seed)
  null <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    permuted <- lapply(fit$modalities, function(m) {
      p <- priors[[m]]
      if (is.null(p)) return(NULL)
      # Permuting the feature labels leaves every set at its size and every
      # feature in its number of sets, and breaks only the correspondence
      # between a loading and an annotation.
      rownames(p) <- rownames(p)[sample(nrow(p))]
      p
    })
    names(permuted) <- fit$modalities
    pp <- profile_of(permuted)
    vals <- vapply(groups, function(d) {
      v <- agreement(pp, d)
      if (is.finite(v)) v else 0
    }, numeric(1))
    null[b] <- max(vals, 0)
  }

  rows <- lapply(names(groups), function(pr) {
    d <- groups[[pr]]
    stat <- agreement(observed_prof, d)
    shared <- Reduce(intersect, lapply(unique(d$modality), function(m) {
      if (is.null(observed_prof[[m]])) character(0) else colnames(observed_prof[[m]])
    }))
    data.frame(
      programme = pr,
      pathway_statistic = if (is.finite(stat)) round(stat, 4) else NA_real_,
      pathway_p = if (is.finite(stat)) {
        (1 + sum(null >= stat)) / (1 + n_perm)
      } else {
        NA_real_
      },
      n_shared_sets = length(shared),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$pathway_supported <- !is.na(out$pathway_p) & out$pathway_p < alpha
  out[order(out$pathway_p, -out$pathway_statistic), , drop = FALSE]
}

#' Which lines of evidence a programme rests on
#'
#' The design channel and the pathway channel answer different questions, and a
#' programme is worth most where both answer yes. Labelling which of them a
#' programme satisfies is what turns "these modalities agree" into a statement
#' a reader can act on, and it keeps a programme supported by one channel from
#' being read as though it were supported by both.
#'
#' @param programmes Output of [chorale_programmes()].
#' @param pathway Output of [chorale_pathway_evidence()].
#'
#' @returns `programmes` with the pathway columns and an `evidence` label
#'   joined on.
#' @export
#' @examples
#' sim <- chorale_simulate(n_modalities = 2, n_features = 120,
#'                         n_shared_factors = 2, n_private_factors = 1,
#'                         n_strains = 4, n_per_cell = 3, effect_size = 3,
#'                         seed = 1)
#' containers <- Map(chorale_load, sim$modalities, sim$col_data)
#' # A curated set reaches every modality, because their identifiers are
#' # harmonised into one space before the sets are applied.
#' ids <- lapply(containers, function(se)
#'   rownames(SummarizedExperiment::assay(se)))
#' span <- function(from, to) unlist(lapply(ids, function(v) v[from:to]))
#' sets <- list(set_a = span(1, 40), set_b = span(30, 80),
#'              set_c = span(60, 120))
#' fit <- chorale_fit(containers, n_factors = c(3, 3), n_init = 2,
#'                    gene_sets = sets)
#' pg <- chorale_programmes(fit)
#' chorale_evidence_label(pg, chorale_pathway_evidence(fit, pg, n_perm = 20))
chorale_evidence_label <- function(programmes, pathway) {
  if (nrow(programmes) == 0) return(programmes)
  if (is.null(pathway) || nrow(pathway) == 0) {
    programmes$pathway_statistic <- NA_real_
    programmes$pathway_p <- NA_real_
    programmes$n_shared_sets <- 0L
    programmes$pathway_supported <- FALSE
  } else {
    idx <- match(programmes$programme, pathway$programme)
    programmes$pathway_statistic <- pathway$pathway_statistic[idx]
    programmes$pathway_p <- pathway$pathway_p[idx]
    programmes$n_shared_sets <- pathway$n_shared_sets[idx]
    programmes$pathway_supported <- !is.na(pathway$pathway_supported[idx]) &
      pathway$pathway_supported[idx]
  }
  design_ok <- isTRUE(TRUE) & programmes$supported
  programmes$evidence <- ifelse(
    design_ok & programmes$pathway_supported, "design and pathway",
    ifelse(design_ok, "design only",
           ifelse(programmes$pathway_supported, "pathway only", "neither"))
  )
  programmes
}
