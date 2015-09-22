

sampW <- function(uProdBeta, xiDay, pregDayBool, pregCycIdx) {
  
  lamDay <- xiDay[pregDayBool] * exp( uProdBeta[pregDayBool] )

  # Calculate in 2 steps b/c we need lamList later
  lamList <- lapply(X=pregCycIdx, FUN=function(x) lamDay[x])
  meanVec <- sapply(X=lamList, FUN=sum)
  
  wSumVec <- sampPoisZeroTr(lambda=meanVec)
  wVals <- sampWDay(nDraws=wSumVec, probs=lamList)
}




# Sample from a zero-truncated Poisson dist ------------------------------------

sampPoisZeroTr <- function(lambda) {
  qpois(runif(n=length(lambda), min=exp(-lambda), max=1), lambda)
}




# Sample W from multinom probs -------------------------------------------------
#
# Expects length q vector 'nDraws' of number of draws per multinomial sample, and
# length q list 'probs' of (unnormalized) probabilities for each bin
#
#   * 'rmultinom' normalizes probs, thus we don't have to beforehand
#   * 'rmultinom' handles the 1 day case as desired

sampWDay <- function(nDraws, probs) {
  unlist( mapply(FUN=rmultinom, n=1, size=nDraws, prob=probs) )
}