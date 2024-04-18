#' rank.matrix
#'
#' Used to convert the expression matrix to the rank matrix of the expression
#'
#' @param mat_ Expression matrix, which is expected to have Gene Symbol represented by rownames and samples represented by colnames.
#'
#' @importFrom future.apply future_sapply
#'
#' @keywords internal
rank.matrix = function(mat_) {

  myApply_ = ifelse(nbrOfWorkers() > 1, future_sapply, sapply)
  rankmatrix_ = myApply_(1:ncol(mat_), \(i__) rank(mat_[, i__]))

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
#'
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#'
#' @keywords internal
delta.rank = function(mat_, net_) {

  message('Rank transforming...')
  mat_ = rank.matrix(mat_)
  message('Processed! Continue deltaRanking...')

  myApply_ = ifelse(nbrOfWorkers() > 1, future_lapply, lapply)
  deltarank_ = myApply_(
    X = 1:nrow(net_),
    FUN = function(i__) {
      r1__ = which(rownames(mat_) == net_[[1]][i__])
      r2__ = which(rownames(mat_) == net_[[2]][i__])
      r__  = mat_[r1__, ] - mat_[r2__, ]

      return(list(as.matrix(net_[i__, ]), r__))
    }
  )
  message('Processed')

  genePair_ = do.call(rbind, sapply(deltarank_, '[', 1))
  rankValue_ = do.call(cbind, sapply(deltarank_, '[', 2))
  genePair_ = cbind(genePair_, paste0('V', 1:nrow(genePair_)))

  return(list(genePair_, rankValue_))

}
