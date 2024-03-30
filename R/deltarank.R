rank.matrix = function(mat_) {

  rankmatrix = c()
  for (i in 1: ncol(mat_)){
    temp = rank(mat_[, i])
    rankmatrix = cbind(rankmatrix, temp)
  }
  colnames(rankmatrix) = colnames(mat_)
  rownames(rankmatrix) = rownames(mat_)

  return(rankmatrix)

}

delta.rank = function(mat_, net_, nThreads_ = 1) {

  mat_ = rank.matrix(mat_)

  calculator_ = function(i__) {
    r1__ = which(rownames(mat_) == net_[[1]][i__])
    r2__ = which(rownames(mat_) == net_[[2]][i__])
    r__  = mat_[r1__, ] - mat_[r2__, ]

    return(list(as.matrix(net_[i__, ]), r__))
  }

  if (nThreads_ != 1) {
    deltarank = future.apply::future_lapply(1:nrow(net_), calculator_)
  } else {
    deltarank = lapply(1:nrow(net_), calculator_)
  }

  genePair = do.call(rbind, sapply(deltarank, '[', 1))
  rankValue = do.call(cbind, sapply(deltarank, '[', 2))
  genePair = cbind(genePair, paste0('V', 1:nrow(genePair)))

  return(list(genePair, rankValue))

}