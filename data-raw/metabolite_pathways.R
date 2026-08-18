# Builds inst/extdata/metabolite_pathways.tsv, the bridge that places a
# lipidome in the same pathway vocabulary the transcriptome and proteome are
# described in.
#
# Nothing published maps lipid shorthand classes onto the Reactome sets used
# for genes. The two halves exist separately: Goslin parses shorthand names
# into LIPID MAPS classes (Kopczynski et al., Anal Chem 92:10957, 2020), and
# Reactome annotates pathways with ChEBI compounds. This script joins them
# through the ChEBI ontology, which is what the pathway layer needs and what
# neither source provides on its own.
#
# Sources, all redistributable:
#   ChEBI ontology, EBI, CC BY 4.0
#     https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo.gz
#   ChEBI-to-Reactome mapping, Reactome, CC0
#     https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt
#   MSigDB Reactome collection, via msigdbr, for the target vocabulary
#
# Run from the package root with the project environment active. It reaches the
# network, so it is not run during checks; the table it writes is committed.

library(msigdbr)

dest <- "inst/extdata/metabolite_pathways.tsv"
cache <- Sys.getenv("CHORALE_SOURCE_CACHE",
                    file.path(tempdir(), "chorale-metabolite-sources"))
dir.create(cache, showWarnings = FALSE, recursive = TRUE)

fetch <- function(url, file) {
  path <- file.path(cache, file)
  if (!file.exists(path)) {
    utils::download.file(url, path, mode = "wb", quiet = FALSE)
  }
  path
}

# ---- The curated half: lipid shorthand class to its ChEBI class term --------
#
# One row per class abbreviation as it appears in MS-DIAL, LIPID MAPS and
# Goslin output. The ChEBI term is the class, not a species: its descendants in
# the ontology are the molecules the class contains. Every identifier is
# checked against the ontology below, so a wrong or retired one fails loudly.

classes <- read.table(text = "
abbreviation	chebi_id	chebi_name	lipidmaps_category
FA	CHEBI:35366	fatty acid	FA
CAR	CHEBI:17387	O-acylcarnitine	FA
MG	CHEBI:17408	monoacylglycerol	GL
DG	CHEBI:17815	1,2-diacyl-sn-glycerol	GL
TG	CHEBI:17855	triglyceride	GL
MGDG	CHEBI:17615	1,2-diacyl-3-beta-D-galactosyl-sn-glycerol	GL
DGDG	CHEBI:28396	3-[alpha-D-galactosyl-(1->6)-beta-D-galactosyl]-1,2-diacyl-sn-glycerol	GL
PA	CHEBI:16337	phosphatidic acid	GP
PC	CHEBI:64482	phosphatidylcholine	GP
PE	CHEBI:16038	phosphatidylethanolamine	GP
PS	CHEBI:18303	phosphatidyl-L-serine	GP
PG	CHEBI:17517	phosphatidylglycerol	GP
PI	CHEBI:28874	phosphatidylinositol	GP
CL	CHEBI:28494	cardiolipin	GP
BMP	CHEBI:60815	lysobisphosphatidic acid	GP
LPA	CHEBI:132742	lysophosphatidic acid	GP
LPC	CHEBI:60479	lysophosphatidylcholine	GP
LPE	CHEBI:64574	lysophosphatidylethanolamine	GP
LPS	CHEBI:68510	lysophosphatidylserine	GP
LPI	CHEBI:64931	lysophosphatidyl-1D-myo-inositol	GP
LPG	CHEBI:90454	lysophosphatidylglycerol	GP
Cer	CHEBI:17761	ceramide	SP
SM	CHEBI:64583	sphingomyelin	SP
HexCer	CHEBI:23079	cerebroside	SP
AHexCer	CHEBI:23079	cerebroside	SP
SHexCer	CHEBI:36477	sulfoglycosphingolipid	SP
GM3	CHEBI:84118	ganglioside GM3	SP
GM1	CHEBI:61048	ganglioside GM1	SP
GD1a	CHEBI:28892	ganglioside	SP
So	CHEBI:16393	sphingosine	SP
Sph	CHEBI:16393	sphingosine	SP
SoP	CHEBI:37550	sphingosine 1-phosphate	SP
CE	CHEBI:17002	cholesteryl ester	ST
ST	CHEBI:15889	sterol	ST
BA	CHEBI:3098	bile acid	ST
CoQ	CHEBI:16389	ubiquinones	PR
VitE	CHEBI:33234	vitamin E	PR
", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

# ---- ChEBI ontology: class term to every molecule beneath it ---------------

obo <- fetch(
  "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo.gz",
  "chebi_lite.obo.gz"
)
lines <- readLines(gzfile(obo))

# The ontology has upwards of two hundred thousand terms, so it is parsed by
# position rather than by walking the file: every line belongs to the block
# opened by the most recent [Term], which one pass of cumsum records.
block <- cumsum(lines == "[Term]")
is_id <- startsWith(lines, "id: CHEBI:")
is_name <- startsWith(lines, "name: ")
is_isa <- startsWith(lines, "is_a: CHEBI:")

block_id <- stats::setNames(substring(lines[is_id], 5), block[is_id])
term_id <- unname(block_id)

name_block <- as.character(block[is_name])
name_value <- substring(lines[is_name], 7)
keep <- !duplicated(name_block) & name_block %in% names(block_id)
term_name <- stats::setNames(name_value[keep], block_id[name_block[keep]])
term_name <- term_name[term_id[term_id %in% names(term_name)]]

isa_block <- as.character(block[is_isa])
edge_child <- unname(block_id[isa_block])
edge_parent <- sub(" .*$", "", substring(lines[is_isa], 7))
ok <- !is.na(edge_child)
edge_child <- edge_child[ok]
edge_parent <- edge_parent[ok]

message("ChEBI terms read: ", length(term_name))

missing <- setdiff(classes$chebi_id, term_id)
if (length(missing) > 0) {
  stop("ChEBI identifiers absent from the ontology: ",
       paste(missing, collapse = ", "))
}
mismatch <- classes$chebi_name != unname(term_name[classes$chebi_id])
if (any(mismatch)) {
  stop("ChEBI names disagree with the ontology for: ",
       paste(classes$abbreviation[mismatch], collapse = ", "), " (ontology says ",
       paste(term_name[classes$chebi_id[mismatch]], collapse = ", "), ")")
}

`%||%` <- function(a, b) if (is.null(a)) b else a
children <- split(edge_child, edge_parent)
descendants <- function(root) {
  seen <- root
  stack <- root
  while (length(stack) > 0) {
    node <- stack[[1]]
    stack <- stack[-1]
    kids <- setdiff(children[[node]] %||% character(), seen)
    seen <- c(seen, kids)
    stack <- c(stack, kids)
  }
  seen
}
# ---- Reactome: which pathways those molecules take part in ----------------

r2c <- fetch(
  "https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt",
  "ChEBI2Reactome_All_Levels.txt"
)
reactome <- utils::read.delim(r2c, header = FALSE, quote = "",
                              stringsAsFactors = FALSE)
names(reactome)[1:6] <- c("chebi", "pathway_id", "url", "pathway_name",
                          "evidence", "species")
reactome <- reactome[reactome$species == "Mus musculus", , drop = FALSE]
reactome$chebi <- paste0("CHEBI:", reactome$chebi)
message("mouse ChEBI-to-Reactome rows: ", nrow(reactome))

# ---- The target vocabulary: sets the gene side is described in -------------
#
# A pathway only earns a row if the gene side also carries it as a set within
# the size window chorale_genesets() applies. That is what makes the two
# modalities land in one coordinate system: a metabolite is placed in a set
# genes are placed in, not in a vocabulary of its own.

gs <- msigdbr(db_species = "MM", species = "Mus musculus",
              collection = "M2", subcollection = "CP:REACTOME")
gs_size <- table(gs$gs_name)
vocabulary <- names(gs_size)[gs_size >= 10 & gs_size <= 500]
as_msigdb <- function(x) paste0("REACTOME_", toupper(gsub("[^A-Za-z0-9]+", "_", x)))

# How many distinct compounds each pathway carries at all. A class that
# accounts for a large share of a pathway's compounds is telling us something
# about that pathway; one generic compound inside a pathway annotated with
# hundreds is not, and the share is what separates them.
pathway_total <- tapply(reactome$chebi, reactome$pathway_id,
                        function(v) length(unique(v)))

rows <- list()
for (i in seq_len(nrow(classes))) {
  d <- descendants(classes$chebi_id[i])
  hit <- reactome[reactome$chebi %in% d, , drop = FALSE]
  if (nrow(hit) == 0) next
  hit$msigdb_name <- as_msigdb(hit$pathway_name)
  hit <- hit[hit$msigdb_name %in% vocabulary, , drop = FALSE]
  if (nrow(hit) == 0) next
  agg <- aggregate(list(n_compounds = hit$chebi),
                   by = list(pathway_id = hit$pathway_id,
                             pathway_name = hit$pathway_name,
                             msigdb_name = hit$msigdb_name),
                   FUN = function(v) length(unique(v)))
  rows[[length(rows) + 1]] <- data.frame(
    abbreviation = classes$abbreviation[i],
    lipidmaps_category = classes$lipidmaps_category[i],
    chebi_id = classes$chebi_id[i],
    chebi_name = classes$chebi_name[i],
    agg[, c("pathway_id", "pathway_name", "msigdb_name")],
    n_compounds = agg$n_compounds,
    n_pathway_compounds = as.numeric(pathway_total[agg$pathway_id]),
    n_descendants = length(d),
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)
out$specificity <- round(out$n_compounds / out$n_pathway_compounds, 4)
out <- out[order(out$abbreviation, -out$specificity, -out$n_compounds,
                 out$msigdb_name), ]
utils::write.table(out, dest, sep = "\t", row.names = FALSE, quote = FALSE)
message("wrote ", dest, ": ", nrow(out), " rows, ",
        length(unique(out$abbreviation)), " classes, ",
        length(unique(out$msigdb_name)), " pathways")
