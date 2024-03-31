#' transGI
#'
#' Nonparametric transformation method of transcriptome based on prior gene interactions
#'
#' @param testMat_ Expression matrix, which is expected to have Gene Symbol represented by rownames and samples represented by colnames.
#' @param method_ Transformation method, currently consisting of either DeltaRank or Pairwise.
#' For detailed information, please refer to \code{\link{delta.rank}} or \code{\link{pair.wise}}.
#' @param bgNet_ Background network type. Currently, prior gene-gene interactions from Reactome and STRING databases are included.
#' You can specify to use either 'reactome' or 'string' by passing the respective argument.
#' If you wish to use custom gene-gene interaction pairs, please pass a dataframe containing two columns of genes.
#' @param nThreads_ Threads to use for transformations, the recommended number of  is between 3 and 6.
#' @param maskMat_ A binary matrix, with the same dimension as \code{testMat_}, is integrated-term-by-term with \code{testMat_},
#' which means that the values of the corresponding positions in \code{testMat_} are penalized to the minimum value for function-based integration.
#' @param controlMat_ Expression matrix with the same dimensions as \code{testMat_}, usually normal tissue.
#' The mean value of expression in the matrix and the corresponding transformation result are subtracted from the final result
#' as the level of gene expression perturbation in the control group.
#'
#' @importFrom future plan
#'
#' @export
#'
#' @examples
#' testMat = system.file('extdata', 'inputMatTest.csv', package = 'transGI') |>
#'   read.csv(row.names = 'symbol') |>
#'   as.matrix()
#' result = transGI(testMat, 'deltarank', 'reactome')
transGI = function(testMat_, method_ = 'deltarank', bgNet_ = 'reactome', nThreads_ = 1, maskMat_ = NULL, controlMat_ = NULL) {

  match.arg(method_, c('deltarank', 'pairwise'))
  if (is.character(bgNet_)) {
    match.arg(bgNet_, c('reactome', 'string'))
    bgNet_ = system.file("extdata", paste0(bgNet_, '.csv'), package = "transGI") |>
      read.csv()
  }
  bgNet_ = bgNet_[(bgNet_[[1]] %in% rownames(testMat_)) & (bgNet_[[2]] %in% rownames(testMat_)), ]

  if (nrow(bgNet_) == 0) stop('pity')

  if (length(maskMat_) != 0) testMat_[matMask_ == 1] = min(testMat_)

  future::plan('multisession', workers = nThreads_)
  res_ = switch(
    method_,
    "deltarank" = delta.rank(testMat_, bgNet_, nThreads_),
    "pairwise"  = pair.wise(testMat_, bgNet_, nThreads_)
  )
  future::plan('sequential')

  if (!is.null(controlMat_)) {

    controlMat_ = controlMat_ |>
      apply(1, mean) |>
      as.matrix()

    colnames(controlMat_) = 'meanPerturbation'

    res_control_ = switch(
      method_,
      "deltarank" = delta.rank(controlMat_, bgNet_),
      "pairwise"  = pair.wise(controlMat_, bgNet_)
    )

    res_[[2]] = res_[[2]] - matrix(rep(unlist(res_control_[[2]]), each = nrow(res_[[2]])), nrow = nrow(res_[[2]]))

  }

  names(res_) = c('genepair', 'converted')
  res_ = lapply(res_, as.data.frame)

  return(res_)

}


