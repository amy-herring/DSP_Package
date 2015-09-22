
# Tracks observation indices for id, sex, preg, U ------------------------------

getIdx <- function() {
  
  cycIdx <- sampCycIdx(pregCycBool, uBeta, xiDay, cycPermsIdx, sexPriorLik)
  pregCycIdx <- getPregCycIdx(cycIdx, pregCycBool)
  
  sexIdx <- list( all=unlist(cycFullIdx) )
  sexBool <- replace(logical(nObs), list=sexIdx, values=TRUE)
  sexPregIdx <- which(sexBool & pregDayBool)
  
  subjIdx <- lapply(naSexDat$idIdx, function(j) j[ sexBool[j] ])
  uIdx <- list( all  = lapply(naSexDat$uBool$all, function(x) which(x & sexBool)),
                preg = lapply(naSexDat$uBool$preg, function(x) which(x & sexBool)) )
  
  outObj <- list( sex = list( all  = sexIdx,
                              preg = sexPregIdx ),
                  pregCyc = pregCycIdx,
                  U = uIdx,
                  subj = subjIdx )
  
  return (outObj)
}





