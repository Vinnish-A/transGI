pair.wise = function(mat_, net_, nThreads_ = 1) {

  calculator_ = function(i__) {
    r1__ = which(rownames(mat_) == net_[[1]][i])
    r2__ = which(rownames(mat_) == net_[[2]][i])
    r__  = as.numeric(mat_[r1, ] > mat_[r2, ])

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
