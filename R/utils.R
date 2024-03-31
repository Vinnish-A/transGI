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

