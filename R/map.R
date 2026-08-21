#' Harmonise feature identifiers to Entrez gene identifiers
#'
#' Maps a vector of feature identifiers (gene symbols, Ensembl or protein
#' identifiers) to Entrez gene identifiers through an organism annotation
#' package. One-to-many mappings are retained rather than resolved arbitrarily:
#' every matching Entrez identifier is returned, each carrying a fractional
#' weight of `1 / n`, where `n` is the number of Entrez identifiers that input
#' identifier maps to.
#'
#' This is how a gene or protein modality reaches the vocabulary of named
#' concepts that connects modalities measured on different individuals: a
#' concept is a set of features, and a feature reaches it through an identifier
#' the two have in common.
#'
#' The organism is an argument rather than a fixture. Any Bioconductor
#' `OrgDb` will do, so a study of human, rat or any other annotated organism
#' maps its features without changing the package. Metabolite identifiers reach
#' the same concepts through their class rather than through gene identifiers,
#' and are handled by [chorale_metabolite_matrix()].
#'
#' @param ids Character vector of feature identifiers to map.
#' @param from Character scalar, the key type of `ids` in the annotation
#'   package, for example `"SYMBOL"`, `"ENSEMBL"` or `"UNIPROT"`.
#' @param orgdb The organism annotation package to map through, named or
#'   supplied as an `OrgDb` object. Defaults to mouse.
#'
#' @returns A data frame with one row per matched (input identifier, Entrez
#'   identifier) pair, columns `id` (the input identifier), `ENTREZID`, and
#'   `weight` (`1 / n` for an input identifier with `n` matches). Input
#'   identifiers with no match are omitted; call with `unique(ids)` first and
#'   compare to the input to build a full per-modality mapping report.
#' @export
#' @examplesIf rlang::is_installed("org.Mm.eg.db")
#' chorale_map(c("Bdnf", "Trem2", "not_a_real_gene"), from = "SYMBOL")
#' # A human study maps through the human annotation instead.
#' # chorale_map(c("BDNF", "TREM2"), from = "SYMBOL", orgdb = "org.Hs.eg.db")
chorale_map <- function(ids, from = "SYMBOL", orgdb = "org.Mm.eg.db") {
  if (!is.character(ids)) {
    rlang::abort("`ids` must be a character vector.")
  }

  db <- if (is.character(orgdb)) {
    rlang::check_installed(orgdb)
    getExportedValue(orgdb, orgdb)
  } else {
    orgdb
  }

  matched <- AnnotationDbi::select(
    db,
    keys = unique(ids),
    keytype = from,
    columns = "ENTREZID"
  )
  matched <- matched[!is.na(matched$ENTREZID), , drop = FALSE]
  names(matched)[names(matched) == from] <- "id"

  if (nrow(matched) == 0) {
    return(data.frame(id = character(), ENTREZID = character(),
                      weight = numeric()))
  }

  n_matches <- stats::ave(matched$id, matched$id, FUN = length)
  matched$weight <- 1 / as.numeric(n_matches)
  matched[, c("id", "ENTREZID", "weight")]
}
