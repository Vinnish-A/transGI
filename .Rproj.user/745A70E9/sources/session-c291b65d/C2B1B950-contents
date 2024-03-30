transGI = function(testMat_, method_ = 'deltarank', bgNet_ = 'reactome', nThreads_ = 1, maskMat_ = NULL, controlMat_ = NULL) {

  match.arg(method_, c('deltarank', 'pairwise'))
  if (is.character(bgNet_)) {
    match.arg(bgNet_, c('reactome', 'string'))
    bgNet_ = system.file("extdata", paste0(bgNet_, '.csv'), package = "transGI") |>
      read.csv()
  }
  bgNet_ = bgNet_[(bgNet_[[1]] %in% rownames(testMat_)) & (bgNet_[[2]] %in% rownames(testMat_)), ]

  if (nrow(bgNet_) == 0) stop('pity')

  if (length(matsMask_) != 0) testMat_[matMask_ == 1] = min(testMat_)

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


