#' Registry of curated gene-set collections
#'
#' Each entry names one collection and records how to retrieve it from
#' MSigDB. The registry is the exchange point of the gene-set layer: adding an
#' entry, or passing a replacement registry to [chorale_genesets()], changes
#' which vocabulary the factors are described and compared in, without
#' touching the estimator.
#'
#' The default registry is deliberately organism-agnostic and disease-agnostic.
#' Collections tied to one disease would name factors in terms of the answer
#' being sought, and the species is a parameter rather than a fixed choice, so
#' the same registry serves a human analysis unchanged.
#'
#' Each collection records its MSigDB code under each source database, since
#' the mouse-native and human collections are numbered differently. Retrieving
#' mouse sets from the mouse database avoids the ortholog mapping that
#' retrieving them from the human database would impose. A collection registered
#' under one database only is retrievable from that database alone, and a
#' request for another is refused by name: KEGG is distributed in the human
#' collections and reaches other organisms through orthologs.
#'
#' @returns A named list, one element per collection. Each element holds
#'   `codes`, a list keyed by MSigDB source database (`"MM"` for mouse-native,
#'   `"HS"` for human) giving the `collection` and `subcollection` understood
#'   by [msigdbr::msigdbr()], and a human-readable `description`.
#' @export
#' @examples
#' names(chorale_geneset_registry())
chorale_geneset_registry <- function() {
  list(
    hallmark = list(
      codes = list(
        MM = list(collection = "MH", subcollection = NULL),
        HS = list(collection = "H", subcollection = NULL)
      ),
      description = "Hallmark: broad, non-redundant biological processes"
    ),
    reactome = list(
      codes = list(
        MM = list(collection = "M2", subcollection = "CP:REACTOME"),
        HS = list(collection = "C2", subcollection = "CP:REACTOME")
      ),
      description = "Reactome: curated pathway granularity"
    ),
    kegg = list(
      codes = list(
        HS = list(collection = "C2", subcollection = "CP:KEGG_LEGACY")
      ),
      description = "KEGG: a separate curated partition of the same biology"
    ),
    cell_type = list(
      codes = list(
        MM = list(collection = "M8", subcollection = NULL),
        HS = list(collection = "C8", subcollection = NULL)
      ),
      description = "Cell type signatures: composition shifts pathways cannot express"
    )
  )
}

#' Retrieve curated sets for description and cross-modality comparison
#'
#' Returns curated gene sets as Entrez identifiers, the space
#' [chorale_map()] harmonises features into. The sets serve two purposes: they
#' supply the vocabulary a factor is named in, and they place factors from
#' different modalities in one biological coordinate system, which is what lets
#' their agreement be assessed on biology as well as on the design.
#'
#' Sets far outside the size window are dropped. Very small sets are matched by
#' chance and very large ones name nothing specific, so both weaken a factor
#' definition rather than sharpening it.
#'
#' @param collections Character vector naming entries of `registry`.
#' @param species Species the identifiers are returned for, as understood by
#'   [msigdbr::msigdbr()], for example `"Mus musculus"` or `"Homo sapiens"`.
#' @param db_species MSigDB source database, `"MM"` for the mouse-native
#'   collections or `"HS"` for the human ones. Switching an analysis from
#'   mouse to human is this argument and `species` together.
#' @param min_size,max_size Retain sets whose membership falls within this
#'   range, counted after mapping to Entrez identifiers.
#' @param registry A registry as returned by [chorale_geneset_registry()].
#' @param custom An optional named list of character vectors of Entrez
#'   identifiers, appended to the retrieved collections. Supplying `custom`
#'   alone, with `collections = character()`, replaces MSigDB entirely.
#'
#' @returns A named list of character vectors of Entrez gene identifiers, with
#'   an attribute `collection` recording which collection each set came from.
#' @export
#' @examplesIf rlang::is_installed("msigdbr")
#' sets <- chorale_genesets("hallmark", species = "Mus musculus")
#' length(sets)
chorale_genesets <- function(collections = c("hallmark", "reactome", "cell_type"),
                             species = "Mus musculus",
                             db_species = c("MM", "HS"),
                             min_size = 10L,
                             max_size = 500L,
                             registry = chorale_geneset_registry(),
                             custom = NULL) {
  db_species <- match.arg(db_species)
  if (length(collections) > 0) {
    rlang::check_installed("msigdbr")
    unknown <- setdiff(collections, names(registry))
    if (length(unknown) > 0) {
      rlang::abort(paste0(
        "Unknown collection(s): ", paste(unknown, collapse = ", "),
        ". Registered: ", paste(names(registry), collapse = ", "), "."
      ))
    }
  }

  sets <- list()
  origin <- character()

  for (nm in collections) {
    spec <- registry[[nm]]$codes[[db_species]]
    if (is.null(spec)) {
      rlang::abort(paste0(
        "Collection '", nm, "' has no code registered for database '",
        db_species, "'."
      ))
    }
    tbl <- msigdbr::msigdbr(
      db_species = db_species,
      species = species,
      collection = spec$collection,
      subcollection = spec$subcollection
    )
    id_col <- if ("ncbi_gene" %in% names(tbl)) "ncbi_gene" else "entrez_gene"
    tbl <- tbl[!is.na(tbl[[id_col]]), , drop = FALSE]
    split_sets <- split(as.character(tbl[[id_col]]), tbl$gs_name)
    split_sets <- lapply(split_sets, unique)
    sets <- c(sets, split_sets)
    origin <- c(origin, rep(nm, length(split_sets)))
  }

  if (!is.null(custom)) {
    if (!is.list(custom) || is.null(names(custom))) {
      rlang::abort("`custom` must be a named list of character vectors.")
    }
    custom <- lapply(custom, function(x) unique(as.character(x)))
    sets <- c(sets, custom)
    origin <- c(origin, rep("custom", length(custom)))
  }

  sizes <- lengths(sets)
  keep <- sizes >= min_size & sizes <= max_size
  sets <- sets[keep]
  origin <- origin[keep]

  if (length(sets) == 0) {
    rlang::abort("No gene sets retained; widen `min_size` and `max_size`.")
  }

  attr(sets, "collection") <- stats::setNames(origin, names(sets))
  sets
}

#' Represent gene sets as a feature-by-set indicator matrix
#'
#' Builds the feature-by-set matrix the pathway layer reads: rows are
#' the features of one modality, columns are gene sets, and an entry is 1 where
#' the feature belongs to the set. Features carry fractional weight where an
#' identifier maps to several genes, so a one-to-many mapping contributes
#' proportionally rather than counting once per gene.
#'
#' @param feature_ids Character vector of Entrez identifiers, one per feature,
#'   in the row order of the modality's assay.
#' @param sets A named list of gene sets, as returned by [chorale_genesets()].
#' @param weights Optional numeric vector, one per feature, giving the
#'   fractional weight from [chorale_map()]. Defaults to 1 for every feature.
#' @param mapping Optional data frame from [chorale_map()]. When supplied,
#'   `feature_ids` are the original assay identifiers and every mapped target is
#'   aggregated with its fractional weight. This preserves one-to-many maps.
#' @param min_features Drop sets matching fewer than this many features of the
#'   modality, since a set that barely intersects the measured features cannot
#'   define a factor in it.
#'
#' @returns A numeric matrix with one row per feature and one column per
#'   retained gene set.
#' @export
#' @examples
#' sets <- list(set_a = c("1", "2", "3"), set_b = c("3", "4"))
#' chorale_geneset_matrix(c("1", "2", "3", "4"), sets, min_features = 2)
chorale_geneset_matrix <- function(feature_ids, sets, weights = NULL,
                                   mapping = NULL,
                                   min_features = 5L) {
  feature_ids <- as.character(feature_ids)
  if (!is.null(mapping)) {
    required <- c("id", "ENTREZID", "weight")
    if (!is.data.frame(mapping) || !all(required %in% names(mapping))) {
      rlang::abort("`mapping` must contain id, ENTREZID and weight columns.")
    }
    mat <- matrix(0, nrow = length(feature_ids), ncol = length(sets),
                  dimnames = list(feature_ids, names(sets)))
    by_id <- split(mapping, as.character(mapping$id))
    for (i in seq_along(feature_ids)) {
      rows <- by_id[[feature_ids[i]]]
      if (is.null(rows) || nrow(rows) == 0L) next
      target <- as.character(rows$ENTREZID)
      weight <- as.numeric(rows$weight)
      for (j in seq_along(sets)) {
        mat[i, j] <- sum(weight[target %in% as.character(sets[[j]])],
                         na.rm = TRUE)
      }
    }
    keep <- colSums(mat > 0) >= min_features
    out <- mat[, keep, drop = FALSE]
    attr(out, "mapping_provenance") <- data.frame(
      n_input = length(feature_ids),
      n_mapped = sum(feature_ids %in% mapping$id),
      n_mapping_rows = nrow(mapping),
      one_to_many = sum(table(mapping$id) > 1L),
      stringsAsFactors = FALSE)
    return(out)
  }
  if (is.null(weights)) weights <- rep(1, length(feature_ids))
  if (length(weights) != length(feature_ids)) {
    rlang::abort("`weights` must have one entry per feature.")
  }

  mat <- vapply(
    sets,
    function(s) ifelse(feature_ids %in% s, weights, 0),
    numeric(length(feature_ids))
  )
  if (is.null(dim(mat))) {
    mat <- matrix(mat, nrow = length(feature_ids),
                  dimnames = list(feature_ids, names(sets)))
  }
  rownames(mat) <- feature_ids

  keep <- colSums(mat > 0) >= min_features
  mat[, keep, drop = FALSE]
}
