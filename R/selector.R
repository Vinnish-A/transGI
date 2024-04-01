
#' @export
selector = c('cal_OR', 'cal_AUROC', 'cal_variance', 'cal_cor', 'cal_HR', 'cal_cindex')

#' cal_OR
#'
#' @param fea_ Features for calculation
#' @param target_ Outcome variable used for calculation, which should be a binary variable for OR and AUROC
#'
#' @export
cal_OR = function(fea_, target_) {

  data_ = data.frame(x = fea_, y = target_)
  fit_  = tryCatch(
    expr  = glm(y ~ x, data = data_, family = "binomial"),
    error = function(e) {c(1, 1, 1, 1)},
    warning = function(w) {c(1, 1, 1, 1)}
  )

  if (is.vector(fit_)) return(setNames(fit_, c("OR", "L95", "H95", "p")))

  summary_ = summary(fit_)
  coef_ = summary_$coefficients[2, ]

  result_ = c(
    OR  = exp(coef_[[1]]),
    L95 = exp(coef_[[1]] - 1.96 * coef_[[2]]),
    H95 = exp(coef_[[1]] + 1.96 * coef_[[2]]),
    p   = coef_[[4]]
  )

  return(result_)

}

#' cal_AUROC
#'
#' @param fea_ Features for calculation
#' @param target_ Outcome variable used for calculation, which should be a binary variable for OR and AUROC
#'
#' @importFrom pROC roc
#'
#' @export
cal_AUROC = function(fea_, target_) {

  result_ = suppressMessages(auc(target_, fea_))
  names(result_) = 'AUROC'

  return(result_)

}

#' cal_variance
#'
#' @param fea_ Features for calculation
#' @param method_ Evaluation methods for dispersion, including variance, p-values of Kruskal or Wilcoxon rank sum test.
#' If using Kruskal or Wilcoxon rank sum test, an additional grouping variable needs to be provided.
#'
#' @export
cal_variance = function(fea_, method_ = c('variance', 'kw', 'wilcoxon'), ...) {

  data_ = data.frame(x = fea_)

  if (method_ %in% c('kw', 'wilcoxon')) {
    data_ = cbind(data_, list(...))
    colnames(data_)[2] = 'group'
  }

  result_ = switch(
    method_,
    variance = c(var = sd(fea_)**2),
    kw = c(p = kruskal.test(x ~ group, data = data_)$p.value),
    wilcoxon = c(p = wilcox.test(x ~ group, data = data_)$p.value)
  )

  return(result_)

}

#' cal_cor
#'
#' @param fea_ Features for calculation
#' @param target_ Outcome variable used for calculation,
#' according to different calculation methods, data types can be integers, continuous numeric, or factors.
#' @param method_ Evaluation methods for correlation, including pearson spearman, and kendall correlation coefficient
#' @param pvalue_ If the test would be executed. Output the p-value if T.
#' @param ... Further arguments to be passed to or from methods.
#'
#' @export
cal_cor = function(fea_, target_, method_ = 'pearson', pvalue_ = F, ...) {

  match.arg(method_, c('pearson', 'spearman', 'kendall'))

  if (pvalue_) {
    result_ = do.call(cor.test, list(x = fea_, y = target_, method = method_, ...))
    result_ = c(result_$estimate, result_$p.value)
    names(result_) = c(method_, 'p')
  } else {
    result_ = do.call(cor, list(x = fea_, y = target_, method = method_, ...))
    names(result_) = method_
  }

  return(result_)

}


#' cal_HR
#'
#' @param fea_ Features for calculation
#' @param time_ Time when the event occurred or the experiment ended
#' @param event_  Whether the termination event occurs or not
#'
#' @import survival
#'
#' @export
cal_HR = function(fea_, time_, event_) {

  data_ = data.frame(x = fea_, time = time_, event = event_)
  fit_  = tryCatch(
    expr  = coxph(Surv(time, event) ~ x, data = data_),
    error = function(e) {c(1, 1, 1, 1)},
    warning = function(w) {c(1, 1, 1, 1)}
  )

  if (is.vector(fit_)) return(setNames(fit_, c("HR", "L95", "H95", "p")))

  summary_ = summary(fit_)
  coef_ = summary_$coefficients[1, ]

  result_ = c(
    HR  = exp(coef_[[1]]),
    L95 = exp(coef_[[1]] - 1.96 * coef_[[3]]),
    H95 = exp(coef_[[1]] + 1.96 * coef_[[3]]),
    p   = coef_[[5]]
  )

  return(result_)

}

#' cal_cindex
#'
#' @param fea_ Features for calculation
#' @param time_ Time when the event occurred or the experiment ended
#' @param event_  Whether the termination event occurs or not
#'
#' @import survival
#' @importFrom Hmisc rcorr.cens
#'
#' @export
cal_cindex = function(fea_, time_, event_) {

  obj_surv_ = Surv(time_, event_)
  result_ = c(cindex = rcorr.cens(fea_, obj_surv_)["C Index"])

  return(result_)

}
