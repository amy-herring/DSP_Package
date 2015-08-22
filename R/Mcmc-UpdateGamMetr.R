

# Sample logical M from full conditional dist ----------------------------------
#
# Full conditional distribution for M_h is given by
#
#       a          1               /  a = p(W | gamma, theta_1, xi) * P(M = 1)
#     -----  =  -------    where  |
#     a + b     1 + b/a            \  b = p(W | gamma, theta_0, xi) * (PM = 0)
#

sampM <- function(Wfull, xiDay, uProdTheta, p) {
  wMean <- list( point = xiDay * exp( uProdTheta$point ),
                 cont  = xiDay * exp( uProdTheta$cont ) )
  
  logTerm <- list( point = getLikW(Wfull, wMean$point) + log(p),
                   cont  = getLikW(Wfull, wMean$cont) + log(1 - p) )
  
  mProb <- 1 / (1 + exp(logTerm$cont - logTerm$point))
  return ( sampBool(mProb) )
}




# Calculate log gamma acceptance ratio r ---------------------------------------
#
# Acceptance ratio for gamma_h is given by
#
#
#            p(W | gamma*, xi, data) * p(gamma_h*)
#         -------------------------------------------
#         p(W | gamma^(s), xi, data) * p(gamma_h^(s))
# 
#
# where gamma* denotes the gamma vector with the h-th term replaced by the proposal 
# value and similarly for gamma^(s)

getGamLogR <- function(Wfull, xiDay, uProdTheta, currTheta, propTheta, hypGamH) {

  meanVec <- list( curr = xiDay * exp( uProdTheta$cont ),
                   prop = xiDay * exp( uProdTheta$prop ) )
  
  wLogLik <- list( curr = getLikW(Wfull, meanVec$curr),
                   prop = getLikW(Wfull, meanVec$prop) )
  
  gamLogLik <- list( curr = dgammaH(currTheta, hypGamH),
                     prop = dgammaH(propTheta, hypGamH) )

  return (wLogLik$prop + gamLogLik$prop - wLogLik$curr - gamLogLik$curr)
}




# Calculate log-lik value for p(W | gamma, xi) ---------------------------------
#
# Note that this is the joint density function, not a vector of the marginal probs

getLikW <- function(W, lambda, logLik=TRUE) {
 logLikVal <- sum( dpois(W, lambda=lambda, log=TRUE) )
 
 if (logLik)
   return (logLikVal)
 else
   return ( exp(logLikVal) )
}




# Calculate likelihood value for p(gamma_h) ------------------------------------

dgammaH <- function(gamVal, hypGamH, logLik=TRUE) {
  list2env(hypGamH, envir=environment())
  
  if (gamVal == 1)
    logLikVal <- log(p)
  else {
    logNrmVal <- log( pgamma(bndU, shape=a, rate=b) - pgamma(bndL, shape=a, rate=b) )
    logTrGamma <- dgamma(gamVal, shape=a, rate=b, log=TRUE) - logNrmVal
    logLikVal <- log(1 - p) + logTrGamma
  }
  
  if (logLik)
    return (logLikVal)
  else
    return ( exp(logLikVal) )
}


