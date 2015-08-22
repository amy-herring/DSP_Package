
# Sample xi from conditional dist ----------------------------------------------

sampXi <- function(W, uProdBeta, phi, idIdx, idPregIdx, pregCycIdx, n) {
  
  a <- phi + replace(integer(n), list=idPregIdx, values=sapply(pregCycIdx, function(x) sum(W[x])))
  b <- phi + sapply(idIdx, function(x) sum(exp(uProdBeta[x])))
  
  return ( rgamma(n=n, shape=a, rate=b) )
}
