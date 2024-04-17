#' rank.matrix
#'
#' Used to convert the expression matrix to the rank matrix of the expression
#'
#' @param mat_ Expression matrix, which is expected to have Gene Symbol represented by rownames and samples represented by colnames.
#'
#' @keywords internal
rank.matrix = function(mat_) {

  rankmatrix_ = sapply(mat_, rank)

  colnames(rankmatrix_) = colnames(mat_)
  rownames(rankmatrix_) = rownames(mat_)

  return(rankmatrix_)

}

#' delta.rank
#'
#' Used to convert the expression matrix to rank matrix first, and then calculate deltarank according to prior gene pairs
#'
#' @param mat_ Expression matrix, which is expected to have Gene Symbol represented by rownames and samples represented by colnames.
#' @param net_ Background network
#' @param nThreads_ Threads to use for transformations, the recommended number of sessions is between 3 and 6.
#'
#' @importFrom future.apply future_lapply
#'
#' @keywords internal
delta.rank = function(mat_, net_, nThreads_ = 1) {

  mat_ = rank.matrix(mat_)

  calculator_ = function(i__) {
    r1__ = which(rownames(mat_) == net_[[1]][i__])
    r2__ = which(rownames(mat_) == net_[[2]][i__])
    r__  = mat_[r1__, ] - mat_[r2__, ]

    return(list(as.matrix(net_[i__, ]), r__))
  }

  if (nThreads_ != 1) {
    future::plan('multisession', workers = nThreads_)
    deltarank_ = future.apply::future_lapply(1:nrow(net_), calculator_)
    future::plan('sequential')
  } else {
    deltarank_ = lapply(1:nrow(net_), calculator_)
  }

  genePair_ = do.call(rbind, sapply(deltarank_, '[', 1))
  rankValue_ = do.call(cbind, sapply(deltarank_, '[', 2))
  genePair_ = cbind(genePair_, paste0('V', 1:nrow(genePair_)))

  return(list(genePair_, rankValue_))

}
