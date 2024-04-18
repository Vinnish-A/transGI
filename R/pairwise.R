#' pair.wise
#'
#' Used to output a binary matrix according to the expression level of the prior gene pairs
#'
#' @param mat_ Expression matrix, which is expected to have Gene Symbol represented by rownames and samples represented by colnames.
#' @param net_ Background network
#'
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#'
#' @keywords internal
pair.wise = function(mat_, net_, nThreads_ = 1) {

  message('Ready to process pairwising...')
  myApply_ = ifelse(nbrOfWorkers() > 1, future_lapply, lapply)
  pairwise_ = myApply_(
    X = 1:nrow(net_),
    FUN = function(i__) {
      r1__ = which(rownames(mat_) == net_[[1]][i])
      r2__ = which(rownames(mat_) == net_[[2]][i])
      r__  = as.numeric(mat_[r1, ] > mat_[r2, ])

      return(list(as.matrix(net_[i__, ]), r__))
    }
  )
  message('Processed')

  genePair_ = do.call(rbind, sapply(pairwise_, '[', 1))
  rankValue_ = do.call(cbind, sapply(pairwise_, '[', 2))
  genePair_ = cbind(genePair_, paste0('V', 1:nrow(genePair_)))

  return(list(genePair_, rankValue_))

}
