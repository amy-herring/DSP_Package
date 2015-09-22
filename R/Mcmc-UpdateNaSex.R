
sampIdx <- function(uProdBeta, xiDay, naSexDat, nObs) {
  
  # 'cycFullIdx': indexes of observations with sex, partitioned by cycle
  cycFullIdx <- sampCycIdx(uProdBeta, xiDay, naSexDat)

  sexIdx <- unlist(cycFullIdx)
  sexBool <- replace(logical(nObs), list=sexIdx, values=TRUE)
  
  sexPregBool <- (sexBool & naSexDat$pregDayBool)
  pregCycIdx <- getPregCycIdx(cycFullIdx, naSexDat$pregCycBool)
  idIdx <- lapply(naSexDat$idIdx, function(j) j[ sexBool[j] ])
  #subjIdIdx <- seq(1:n)[ sapply(idIdx, function(j) TRUE %in% sexBool$all[j]) ]
  uBool <- lapply(naSexDat$uBool, function(x) x & sexBool)
  pregUBool <- lapply(naSexDat$pregUBool, function(x) x & sexBool)
  #pregDayBool <- naSexDat$pregDayBool[sexBool$all]
  
#   cat(length(naSexDat$idIdx), length(idIdx))
  
  # Construct index data structure
  idx <- list( preg    = which(sexPregBool),
               U       = list( all  = uBool,
                               preg = pregUBool ),
               pregCyc = pregCycIdx,
               subj    = list( obs  = idIdx,
                               preg = naSexDat$pregIdx ) )
  
  return (idx)
}




# Sample intercourse values for missing ----------------------------------------

#sampCycIdx <- function(pregCycBool, uBeta, xiDay, cycPermsIdx, sexPriorLik) {
sampCycIdx <- function(uBeta, xiDay, naSexDat) {
  postSexProbs <- getPostSexProbs(uBeta, xiDay, naSexDat)
  mapply(FUN=sample, x=naSexDat$cycPermsIdx, size=1, prob=postSexProbs) 
}




# Unnormalized probabilities for sex imputation --------------------------------

#getPostSexProbs <- function(pregCycBool, uBeta, xiDay, cycPermsIdx, sexPriorLik) {
getPostSexProbs <- function(uBeta, xiDay, naSexDat) {

  getThisProb <- function(thisY, thisPermIdx, thisPermPrior) {
    getExpTerm <- function(idx) {
      if (identical(length(idx), 0L)) 1
      else exp(-xiDay[head(idx, 1)] * sum(exp(uBeta[idx])))
    }
      
    # No missing so prob of the provided sex perm is 1
    if (identical(length(thisPermIdx), 1L)) 
      return (1)
    # else: cycle has a missing and we have to calculate probs
    expTerm <- sapply(thisPermIdx, getExpTerm)
    yLikVec <- if (thisY) 1 - expTerm else expTerm
    return (yLikVec * thisPermPrior)
  }
  
  #return ( with( mapply(getThisProb, pregCycBool, cycPermsIdx, sexPriorLik), naSexDat ) )
  return ( with(naSexDat, mapply(getThisProb, pregCycBool, cycPermsIdx, sexPriorLik)) )
  #return ( mapply(getThisProb, naSexDat$pregCycBool, naSexDat$cycPermsIdx, naSexDat$sexPriorLik) )
}




# All possible permutations of cycle indices -----------------------------------
#
# I.e. all known sex days in a cycle combined with all sets of T/F for missing sex

getCycPermsIdx <- function(cycIdx, naSexBool) {
  
  # Combine a vector with each vector in a list
  combVecWithList <- function(theVec, theList) {
    if (identical(length(theList), 0L)) theVec
    else lapply(theList, function(x) sort(c(x, theVec)))
  }
  
  cycSexPartit <- getCycSexPartit(cycIdx, naSexBool)
  return ( mapply(combVecWithList, cycSexPartit$known, cycSexPartit$perms) )
}




# Enumerate all possible intercourse sets --------------------------------------
#
# Split all cycles into known sex and missing sex, enumerate all possible missing sex sets, 
# and then for each cycle combine the known with each possible enumeration of missing

getCycSexPartit <- function(cycIdx, naSexBool) {
  
  # 'cycSexPartit': for each cycle, a list of indices corresp to the known (or missing) obs
  cycSexPartit <- list( known = lapply(cycIdx, function(j) Filter(function(j) !naSexBool[j], j)),
                        miss  = lapply(cycIdx, function(j) Filter(function(j) naSexBool[j], j)) )
  # 'cycSexPartit$perms': power set for indices of each cycle of missing sex
  cycSexPartit$perms <- lapply(cycSexPartit$miss, getPowSet)
  
  return (cycSexPartit)
}




# Prior likelihoods for all possible sex sets ----------------------------------

getSexPriorLik <- function(cycSexPriors) {
  
  # Calc prior probs for a list of T/F perms of missing sex
  getCycLiks <- function(x) {
    n <- length(x)
    if (identical(n, 0L)) 
      return (1)
    else {
      # :tfPerms': 2^n x n matrix of all T/F perms of length n
      tfPerms <- getTfPerms(length(x))
      xComp <- 1 - x
      # 'permProbs': vectors with elements corresp to P(X_ijk = m) for choice of m=0 or m=1
      permProbs <- lapply(1:2^n, function(j) c( x[tfPerms[j, ]], xComp[!tfPerms[j, ]] ))
      return ( sapply(permProbs, prod) )
    }
  }
  
  return ( lapply(cycSexPriors, getCycLiks) )
}




# Partition the priors for missing into cycles ---------------------------------

getCycSexPriors <- function(cycIdx, naSexBool, naSexPriors) {
  priorFull <- replace(numeric(length(naSexBool)), list=naSexBool, values=naSexPriors)
  lapply(cycIdx, function(j) Filter(as.logical, priorFull[j]))
}




# Compute the power set for vector input ---------------------------------------
# 
# the 'times' parameter is chosen because 2^n / (2 * k) = 2^(n - 1) / k and the 
# product of twice 'each' and 'times' must be to equal 2^n

getPowSet <- function(set) {
  if (identical(length(set), 0L)) 
    return (list(numeric(0)))
  
  # else: non-empty set
  n <- length(set)
  keepBool <- getTfPerms(n)
  return ( lapply(1:2^n, function(j) set[keepBool[j, ]]) )
}




# 2^n x n matrix of all size n perms of T/F ------------------------------------
#
# PRE: n is integer, n >= 1

getTfPerms <- function(n) {
  sapply(2^(1:n - 1), function(k) 
    rep(c(FALSE, TRUE), each=k, times=(2^(n - 1) / k)))
}




# Index partitioned by subject -------------------------------------------------

getCurrIdx <- function(fullIdx, sexBool) {
  idx <- lapply(fullIdx, function(currSubjIdx) Filter(function(k) sexBool[k], currSubjIdx))
  #return ( Filter(length, idx) )
}




# Update 'pregCycIdx' after inputing missing sex -------------------------------

getPregCycIdx <- function(cycFullIdx, pregCycBool) {
  pregCycFullIdx <- cycFullIdx[pregCycBool]
  pregFullIdx <- unlist(pregCycFullIdx)
  nObs <- length(pregFullIdx)
  newIdx <- replace(integer(nObs), pregFullIdx, 1:nObs)
  return (lapply(pregCycFullIdx, function(j) newIdx[j]))
}

