
# Create objects for use in the MCMC sampler -----------------------------------

getSamplerObj <- function(modelObj, fwLen) {
  # Contains 'Y', 'X', 'id', 'U'
  list2env(modelObj, envir=environment())
  
  pregDayBool <- rep(convToBool(Y), each=fwLen)
  sexBool <- convToBool(X)
  sexPregBool <- (sexBool & pregDayBool)
  
  n <- length(unique(id[sexBool])) # number of individuals
  q <- ncol(U) # number of covariates
  
  # Reduce data to days in which intercourse occured ---------------------------
  
  # Correspond to the rows after reducing data via sexBool / sexPregBool
  subSexRows <- replace(rep(0, length(id)), list=which(sexBool), values=1:sum(sexBool))
  subSexPregRows <- replace(rep(0, length(id)), list=which(sexPregBool), values=1:sum(sexPregBool))
  
  # Elements are indices corresponding to a subject
  idIdx <- Filter(length, lapply(X=unique(id), FUN=function(x) subSexRows[(id == x) & sexBool]))
  idDayExpan <- rep(1:length(idIdx), times=sapply(idIdx, length))
  
  # Elements are indices of cycles that have pregnancy
  pregCycIdx <- Filter(length, lapply(X=unique(id), FUN=function(x) 
    subSexPregRows[(id == x) & sexBool & pregDayBool]))
  
  # Reduce objects to intercourse days
  U <- U[sexBool, ]
  pregDayBool <- pregDayBool[sexBool]
  # Indexes the subjects who have a pregnancy; used for updating xi
  idPregIdx <- which( tapply(pregDayBool, INDEX=id[sexBool], FUN=function(x) TRUE %in% x) )
  
  # Convert binary cols of U to boolean
  covIsBinVec <- apply(U, MARGIN=2, FUN=function(x) isTRUE(all.equal(names(table(x)), c("0","1"))))
  uBool <- lapply(1:q, function(j) if (!covIsBinVec[j]) NULL else (U[, j] == 1))
  pregUBool <- lapply(uBool, function(x) x[pregDayBool])
  
  # ----------------------------------------------------------------------------
  
  samplerObj <- list( U = U,
                      pregDayBool = pregDayBool,
                      pregCycIdx = pregCycIdx,
                      uBool = uBool,
                      pregUBool = pregUBool,
                      idIdx = idIdx,
                      idDayExpan = idDayExpan,
                      idPregIdx = idPregIdx,
                      n = n,
                      q = q,
                      subjId = unique(id[sexBool]),
                      varNames = colnames(U) )
  return (samplerObj)
}