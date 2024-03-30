#' read_gmt
#'
#' Used to read in gmt file downloaded form [MSigDB](https://www.gsea-msigdb.org/gsea/index.jsp)
#'
#' @param gmtfile_ Path to gmt formatted file.
#' @param filteron_ Filter terms without keywords like \code{GOBP}.
#'
#' @export
#'
#' @examples
#' db = read_gmt('yourpath/h.all.v2023.2.Hs.symbols.gmt')
read_gmt = function(gmtfile_, filteron_ = NULL) {

  x = readLines(gmtfile_)
  res = strsplit(x, "\t")
  names(res) = vapply(res, function(y) y[1], character(1))
  res = lapply(res, "[", -c(1:2))
  ont2gene = stack(res)
  ont2gene = ont2gene[, c("ind", "values")]
  colnames(ont2gene) = c("term", "gene")

  if (!is.null(filteron_)) ont2gene = ont2gene[grepl(filteron_, ont2gene$term), ]

  return(as.data.frame(ont2gene))

}

#' enrichGraphInner
#'
#' Used to perform enrichment of a single gene set in a prior gene interaction pair.
#'
#' @keywords internal
enrichGraphInner = function(inputGenes_, geneSet_, bgNet_, bgGeneNum_) {

  geneEnriched_ = intersect(geneSet_, inputGenes_)

  if (length(geneEnriched_) == 0) return(c(ES = 0, p = 1))

  inputNet_ = bgNet_ |>
    filter(from %in% inputGenes_ & to %in% inputGenes_)
  pathwayNet_ = bgNet_ |>
    filter(from %in% geneSet_ & to %in% geneSet_)
  disturbedNet_ = bgNet_ |>
    filter(from %in% geneEnriched_ & to %in% geneEnriched_)

  ES_ = length(geneEnriched_)*nrow(disturbedNet_)/nrow(pathwayNet_); if (is.infinite(ES_)) ES_ = NA
  pHyper_ = phyper(nrow(disturbedNet_)-1, nrow(pathwayNet_), nrow(bgNet_)-nrow(pathwayNet_), nrow(inputNet_), lower.tail = F)

  return(c(ES = ES_, p = pHyper_))

}

#' enrichGraph
#'
#' @param inputGenes_ A string of genes you want to enrich.
#' @param db_ Data frame of the enriched item, the first column is term and the second column is gene.
#' @param pathways_ The specific path term used to perform enrichment, pass \code{all} to enrich all terms.
#' @param bgNet_ Background network type. Currently, prior gene-gene interactions from Reactome and STRING databases are included.
#' You can specify to use either 'reactome' or 'string' by passing the respective argument.
#' If you wish to use custom gene-gene interaction pairs, please pass a dataframe containing two columns of genes.
#' @param thres_ FDR threshold of the output result.
#'
#' @importFrom dplyr mutate filter relocate
#'
#' @export
#'
#' @examples
#' db = read_gmt('yourpath/h.all.v2023.2.Hs.symbols.gmt')
#' result = enrichGraph(geneString, db, 'all')
enrichGraph = function(inputGenes_, db_, pathways_ = 'all', bgNet_ = 'reactome', thres_ = 0.05) {

  if (is.character(bgNet_)) {
    match.arg(bgNet_, c('reactome', 'string'))
    bgNet_ = system.file("extdata", paste0(bgNet_, '.csv'), package = "transGI") |>
      read.csv()
  }

  colnames(bgNet_)[1:2] = c('from', 'to')
  if (pathways_ == 'all') pathways_ = db_$term |> unique()

  bgGeneNum_ = length(unique(c(bgNet_$from, bgNet_$to)))

  terms_ = db_ |> filter(term %in% pathways_); terms_ = split(terms_$gene, terms_$term)

  resultLst_ = lapply(terms_, enrichGraphInner, inputGenes_ = inputGenes_, bgNet_ = bgNet_, bgGeneNum_ = bgGeneNum_)

  result_ = do.call(rbind, resultLst_) |>
    as.data.frame() |>
    mutate(pathway = names(terms_)) |>
    filter(ES != 0) |>
    mutate(p.adj = p.adjust(p), q = p.adjust(p, "fdr")) |>
    relocate(pathway, ES) |>
    filter(q <= thres_)

  return(result_)

}
