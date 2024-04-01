#' pToLabel
#'
#' Used to perform enrichment of a single gene set in a prior gene interaction pair.
#'
#' @keywords internal
pToLabel = function(vec__) {
  vec__[is.na(vec__)] = 1
  cut(
    vec__,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("FDR<0.001", "FDR<0.01", "FDR<0.05", "FDR>0.05")
  )
}

#' splitTerms
#'
#' Used to perform enrichment of a single gene set in a prior gene interaction pair.
#'
#' @importFrom stringr str_split str_pad str_detect
#'
#' @keywords internal
splitTerms = function(Des__) {

  DesVec__ = str_split(Des__, " ", simplify = F)[[1]]
  cutPoint__ = ceiling(length(DesVec__)/2)

  DesVec1__ = DesVec__[1:cutPoint__]; DesVec2__ = DesVec__[(cutPoint__+1):length(DesVec__)]
  DesVec__ = c(paste(DesVec1__, collapse = " "), paste(DesVec2__, collapse = " "))

  if(nchar(DesVec__[1]) > nchar(DesVec__[2])) {
    DesVec__ = str_pad(DesVec__, nchar(DesVec__[1]), 'both')
  } else {
    DesVec__ = str_pad(DesVec__, nchar(DesVec__[2]), 'both')
  }

  return(paste(DesVec__, collapse = '\n'))

}

#' net2gene
#'
#' @param genepair_ Genepair with id
#' @param idNet_ Selected genepair ids
#'
#' @importFrom dplyr filter pull
#'
#' @export
net2gene = function(genepair_, idNet_) {

  colnames(genepair_)[1:3] = c('from', 'to', 'id')

  genes_ = genepair_ |>
    filter(id %in% idNet_) |>
    pull(from, to)

  result_ = unique(c(genes_, names(genes_)))

  return(result_)

}

#' autoBar
#'
#' @param df_
#'
#' @import ggplot2
#' @importFrom dplyr mutate arrange slice_max
#' @importFrom purrr map_chr
#' @importFrom stringr str_replace_all str_sub str_detect
#'
#' @export
autoBar = function(df_, n_ = 5) {

  dataPlot_ = df_ |>
    slice_max(ES, n = n_, with_ties = F) |>
    arrange(ES) |>
    mutate(pathway = pathway |> str_sub(10) |> str_replace_all('_', ' '),
           pathway = map_chr(pathway, ~ ifelse(str_detect(.x, ' '), splitTerms(.x), .x)),
           pathway = factor(pathway, levels = pathway),
           qLabel  = pToLabel(q))

  dataPlot_ |>
    ggplot() +
    geom_col(aes(ES, pathway, fill = qLabel)) +
    geom_text(aes(ES, pathway, label = pathway)) +
    scale_fill_manual(values = c("FDR<0.001" = "#4e62ab", "FDR<0.01" = "#479db4", "FDR<0.05" = "#87d0a6", 'FDR>0.05' = "#cbe99d")) +
    scale_x_continuous(expand = c(0, 0, 0, 0.1)) +
    scale_y_discrete(breaks = NULL) +
    labs(fill = '') +
    xlab('Enrichment Score Based on Graph') +
    ylab('') +
    theme_classic() +
    theme(
      legend.position = "top"
    )

}

#' autoBubble
#'
#' @param df_
#'
#' @import ggplot2
#' @importFrom dplyr mutate arrange slice_max
#' @importFrom purrr map_chr
#' @importFrom stringr str_replace_all str_sub
#'
#' @export
autoBubble = function(df_, n_ = 5) {

  dataPlot_ = df_ |>
    slice_max(ES, n = n_, with_ties = F) |>
    arrange(ES) |>
    mutate(pathway = pathway |> str_sub(10) |> str_replace_all('_', ' '),
           pathway = map_chr(pathway, ~ ifelse(str_detect(.x, ' '), splitTerms(.x), .x)),
           pathway = factor(pathway, levels = pathway),
           qLabel  = pToLabel(q))

  dataPlot_ |>
    ggplot() +
    geom_point(aes(ES, pathway, fill = qLabel), shape = 21, size = 5) +
    geom_text(aes(ES, pathway, label = pathway)) +
    scale_fill_manual(values = c("FDR<0.001" = "#4e62ab", "FDR<0.01" = "#479db4", "FDR<0.05" = "#87d0a6", 'FDR>0.05' = "#cbe99d")) +
    scale_x_continuous(expand = c(0.05, 0, 0, 0.1)) +
    scale_y_discrete(breaks = NULL) +
    labs(fill = '') +
    xlab('Enrichment Score Based on Graph') +
    ylab('') +
    theme_classic() +
    theme(
      legend.position = "top"
    )

}

#' visNet
#'
#' @param inputGenes_ A string of genes you want to enrich.
#' @param bgNet_ Background network type. Currently, prior gene-gene interactions from Reactome and STRING databases are included.
#' You can specify to use either 'reactome' or 'string' by passing the respective argument.
#' If you wish to use custom gene-gene interaction pairs, please pass a dataframe containing two columns of genes.
#' @param degree_ An integer. Label nodes with node degree more than it.
#'
#' @import tidygraph
#' @import ggraph
#' @importFrom dplyr mutate filter
#'
#' @export
visNet = function(inputGenes_, bgNet_ = 'reactome', degree_ = 3) {

  if (is.character(bgNet_)) {
    match.arg(bgNet_, c('reactome', 'string'))
    bgNet_ = system.file("extdata", paste0(bgNet_, '.csv'), package = "transGI") |>
      read.csv()
  }

  colnames(bgNet_)[1:2] = c('from', 'to')

  bgNetSig_ = bgNet_ |>
    filter(from %in% inputGenes_ & to %in% inputGenes_)

  graphGI_ = as_tbl_graph(bgNetSig_) |>
    activate(nodes) |>
    mutate(deg = centrality_degree(mode = 'in'))

  ggraph(graphGI_, layout = 'kk') +
    geom_edge_fan(color = "#394c81", alpha = 0.8, show.legend = FALSE) +
    geom_node_point(aes(size = deg), shape = 21, fill = '#94697a', color = 'white') +
    geom_node_label(aes(filter = deg > degree_, label = name), size = 3, repel = T, max.overlaps = Inf) +
    scale_color_discrete() +
    scale_edge_width(range=c(0.2,3)) +
    guides(size = 'none', fill = 'none') +
    theme_graph()

}
