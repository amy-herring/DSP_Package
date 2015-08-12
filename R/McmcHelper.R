
# Combine user input gamma hyperpars with default ------------------------------

getHypGam <- function(varNames, userHypGam) {
  parNames <- c("a","b","p","bndL","bndU")
  defaultVals <- list(a=1, b=1, p=0.5, bndL=0, bndU=Inf)
  
  # Combine user input hyperparameters (when provided) with default hyperparameters
  # Uses 'parNames' and 'defaultVals' from parent environment
  getHypGam <- function(thisUserHyp) {
    if ( is.null(thisUserHyp) )
      return (defaultVals)
    
    # else: custom params provided for this variable
    # Combine default values with custom input values (when given)
    combineHyp <- defaultVals
    for (thisName in parNames)
      if (thisName %in% names(thisUserHyp))
        combineHyp[[thisName]] <- thisUserHyp[[thisName]]
    
    return (combineHyp)
  }
  
  hypGam <- lapply(varNames, function(x) getHypGam(userHypGam[[x]]))
  names(hypGam) <- varNames
  
  return (hypGam)
}

# keys <- unique(c(names(lst1), names(lst2)))
# setNames(mapply(c, lst1[keys], lst2[keys]), keys)       # Try adapting this code




# Combine user input phi hyperpars with default --------------------------------

getHypPhi <- function(userHypPhi) {
  hypPhi <- list(c1=1, c2=1)
  
  if ("c1" %in% names(userHypPhi)) hypPhi$c1 <- userHypPhi$c1
  if ("c2" %in% names(userHypPhi)) hypPhi$c2 <- userHypPhi$c2
  
  return (hypPhi)
}




# Initial values for gamma set to mean of priors -------------------------------

getGamInit <- function(hypGam, gamIsTrunBool) {
  
  # Calculate the mean for a prior distribution of gamma_h
  getPriorMean <- function(thisHypGam, thisGamIsTrBool) {
    if (!thisGamIsTrBool)
      return (thisHypGam$a / thisHypGam$b)
    
    # else: prior dist is truncated and we have to integrate to obtain mean
    
    # The term inside the integral for calculating the mean
    expecFcn <- function(x) x * dgamma(x, shape=thisHypGam$a, rate=thisHypGam$b)
    integralTerm <- integrate(expecFcn, lower=thisHypGam$bndL, upper=thisHypGam$bndU)$value
    normalizeConst <- with(pgamma(bndU, shape=a, rate=b) - pgamma(bndL, shape=a, rate=b), 
                           data=thisHypGam)
    return (integralTerm / normalizeConst)
  }
  
  gamInit <- sapply(1:length(hypGam), function(j) getPriorMean(hypGam[[j]], gamIsTrunBool[j]))
  return (gamInit)
}



