
pToLabel = function(vec__) {
  vec__[is.na(vec__)] = 1
  cut(
    vec__,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("FDR<0.001", "FDR<0.01", "FDR<0.05", "FDR>0.05")
  )
}

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

autoBar = function(df_) {



}

autoBubble = function(df_) {



}
