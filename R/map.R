#' Harmonise feature identifiers to mouse Entrez gene identifiers
#'
#' Maps a vector of feature identifiers (gene symbols, Ensembl or protein
#' identifiers) to mouse Entrez gene identifiers via
#' [org.Mm.eg.db::org.Mm.eg.db]. One-to-many mappings are retained rather
#' than resolved arbitrarily: every matching Entrez identifier is returned,
#' each carrying a fractional weight of `1 / n`, where `n` is the number of
#' Entrez identifiers that input identifier maps to (`AGENT_PLAN.md`
#' Section 8.2). Metabolite identifiers, which map through enzyme and
#' reaction annotations rather than gene identifiers, are out of scope for
#' this function.
#'
#' @param ids Character vector of feature identifiers to map.
#' @param from Character scalar, the key type of `ids` in
#'   [org.Mm.eg.db::org.Mm.eg.db] (for example `"SYMBOL"`, `"ENSEMBL"`,
#'   `"UNIPROT"`).
#'
#' @returns A data frame with one row per matched (input identifier, Entrez
#'   identifier) pair, columns `id` (the input identifier), `ENTREZID`, and
#'   `weight` (`1 / n` for an input identifier with `n` matches). Input
#'   identifiers with no match are omitted; call with `unique(ids)` first and
#'   compare to the input to build a full per-modality mapping report.
#' @export
#' @examples
#' chorale_map(c("Bdnf", "Trem2", "not_a_real_gene"), from = "SYMBOL")
chorale_map <- function(ids, from = "SYMBOL") {
  rlang::check_installed("org.Mm.eg.db")

  if (!is.character(ids)) {
    rlang::abort("`ids` must be a character vector.")
  }

  matched <- AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = unique(ids),
    keytype = from,
    columns = "ENTREZID"
  )
  matched <- matched[!is.na(matched$ENTREZID), , drop = FALSE]
  names(matched)[names(matched) == from] <- "id"

  if (nrow(matched) == 0) {
    return(data.frame(id = character(), ENTREZID = character(), weight = numeric()))
  }

  n_matches <- stats::ave(matched$id, matched$id, FUN = length)
  matched$weight <- 1 / as.numeric(n_matches)
  matched[, c("id", "ENTREZID", "weight")]
}
