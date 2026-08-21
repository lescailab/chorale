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

#' Place a lipidome's features in the vocabulary
#'
#' The counterpart of [chorale_geneset_matrix()] for a modality whose features
#' are lipids. The columns are the same concepts, so a lipidome and a
#' transcriptome end up in one vocabulary rather than in two that would have to
#' be reconciled afterwards.
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
