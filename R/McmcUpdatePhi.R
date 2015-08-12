
# Sample phi from proposal dist ------------------------------------------------

sampPhiProp <- function(phi, delta, phiPropDist="uniform") {
  
  if (phiPropDist == "normal")
    phiProp <- abs( rnorm(n=1, mean=phi, sd=delta) )
  else
    phiProp <- abs( runif(n=1, min=(phi - delta), max=(phi + delta)) )
  
  return (phiProp)
}




# Calculate log(r) -------------------------------------------------------------
#
# r is given by
#
#           ( prod_i pi(xi_i | phiProposal) ) * pi(phiProposal)
#           ---------------------------------------------------
#            ( prod_i pi(xi_i | phiCurrent) ) * pi(phiCurrent)
#

getPhiLogR <- function(xi, phiCurr, phiProp, hypPhi) {
  
  logNumerTerm1 <- sum( dgamma(x=xi, shape=phiProp, rate=phiProp, log=TRUE) )
  logNumerTerm2 <- dgamma(x=phiProp, shape=hypPhi$c1, rate=hypPhi$c2, log=TRUE)
  
  logDenomTerm1 <- sum( dgamma(x=xi, shape=phiCurr, rate=phiCurr, log=TRUE ) )
  logDenomTerm2 <- dgamma(x=phiCurr, shape=hypPhi$c1, rate=hypPhi$c2, log=TRUE)
  
  return( logNumerTerm1 + logNumerTerm2 - logDenomTerm1 - logDenomTerm2 )
}
