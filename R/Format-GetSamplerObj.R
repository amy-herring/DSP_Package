
# Create objects for use in the MCMC sampler -----------------------------------

getSamplerObj <- function(modelObj, cycVec, fwLen, useNA) {
  # Contains 'Y', 'X', 'id', 'U'
  list2env(modelObj, envir=environment())
  
  pregCycBool <- convToBool(Y)
  pregDayBool <- rep(pregCycBool, each=fwLen)
  sexBool <- convToBool(X)
  sexMissBool <- is.na(X)
  sexPregBool <- (sexBool & pregDayBool)
  
  n <- length(unique(id[sexBool]))  # number of individuals
  q <- ncol(U)                      # number of covariates
  
  # Reduce to days in which intercourse occured --------------------------------
  
  # Correspond to the rows after reducing data via sexBool / sexPregBool
  subSexRows <- replace(rep(0, length(id)), list=which(sexBool), values=1:sum(sexBool))
  subSexPregRows <- replace(rep(0, length(id)), list=which(sexPregBool), values=1:sum(sexPregBool))
  
  # Elements are indices corresponding to a subject
  idIdx <- Filter(length, lapply(X=unique(id), FUN=function(x) subSexRows[(id == x) & sexBool]))
  subjId <- id[ sapply(idIdx, head, 1) ]
  idDayExpan <- rep(1:length(idIdx), times=sapply(idIdx, length))
  nObs <- length(idDayExpan)
  
  # Elements are indices of cycles
  sexCycVec <- cycVec[sexBool]
  cycIdx <- lapply(idIdx, function(j) 
    lapply(unique(sexCycVec[j]), function(x) j[sexCycVec[j] == x]))
  cycIdx <- unlist(cycIdx, recursive=FALSE)
  
  # Elements are indices of cycles that have pregnancy
  pregCycIdx <- Filter(length, lapply(X=unique(id), FUN=function(x) 
    subSexPregRows[(id == x) & sexBool & pregDayBool]))
  
  # Reduce objects to intercourse days
  U <- U[sexBool, ]
  pregDayBool <- pregDayBool[sexBool]
  naSexIdx <- which(sexMissBool[sexBool])
  # 'idPregIdx': indexes the subjects who have a pregnancy; used for updating xi
  idPregIdx <- which( tapply(pregDayBool, INDEX=id[sexBool], FUN=function(x) TRUE %in% x) )
  pregCycBool <- sapply(cycIdx, function(j) pregDayBool[j[1]])
  
  # Convert binary cols of U to boolean
  gamIsBinBool <- apply(U, MARGIN=2, FUN=function(x) isTRUE(all.equal(names(table(x)), c("0","1"))))
  uBool <- lapply(1:q, function(j) if (!gamIsBinBool[j]) NULL else (U[, j] == 1))
  pregUBool <- lapply(uBool, function(x) if (is.null(x)) NULL else (x & pregDayBool))
  
  # Convert uBool, pregUBool to index
  uIdx <- lapply(uBool, function(x) if (is.null(x)) NULL else which(x))
  pregUIdx <- lapply(pregUBool, function(x) if (is.null(x)) NULL else which(x))
  
  # Missing data objects -------------------------------------------------------

#   useNaSexBool <- FALSE
#   if (identical(useNA, "sex")) {
#     useNaSexBool <- TRUE
#     # 'naSexDat': 
#     naSexDat <- list( idIdx = idIdx,
#                       pregCycIdx = cycIdx$preg,
#                       #idDayExpan = idDayExpan,
#                       uBool = uBool,
#                       pregUBool = pregUBool,
#                       pregDayBool = pregDayBool,
#                       subjId = subjId )
#     # Objects for use when sampling missing intercourse data
#     naSexBool <- replace(logical(nObs), list=naSexIdx, values=TRUE)
#     cycPermsIdx <- getCycPermsIdx(cycIdx$all, naSexBool)
#     cycSexPriors <- getCycSexPriors(cycIdx$all, naSexBool, naSexPriors)
#     sexPriorLik <- getSexPriorLik(cycSexPriors)
#     #rm(naSexBool, cycPermsIdx)
#   }
  
  
  if (!identical(useNA, "sex")) {
    useNaSexBool <- FALSE
    
    idx <- list( preg    = which(pregDayBool),
                 U       = list( all  = uIdx,
                                 preg = pregUIdx ),
                 pregCyc = pregCycIdx,
                 subj    = list( obs  = idIdx,
                                 preg = idPregIdx ) )
  }
  else {
    useNaSexBool <- TRUE
    
    naSexBool <- replace(logical(nObs), list=naSexIdx, values=TRUE)
    cycPermsIdx <- getCycPermsIdx(cycIdx, naSexBool)
    
    naSexPriors <- getNaSexProbs( which(sexMissBool) )
    cycSexPriors <- getCycSexPriors(cycIdx, naSexBool, naSexPriors)
    sexPriorLik <- getSexPriorLik(cycSexPriors)
    
    naSexDat <- list( pregCycBool = pregCycBool,
                      pregDayBool = pregDayBool,
                      uBool       = uBool,
                      pregUBool   = pregUBool,
                      idIdx       = idIdx,
                      pregIdx     = idPregIdx,
                      cycPermsIdx = cycPermsIdx,
                      sexPriorLik = sexPriorLik )
  }
  
  
  #cat(length(idIdx), length(naSexDat$idIdx))

  
  # Construct return object ----------------------------------------------------
  
#   samplerObj <- list( U = U,
#                       uBool = uBool,
#                       pregUBool = pregUBool,
#                       pregDayBool = pregDayBool,
#                       idIdx = idIdx,
#                       idDayExpan = idDayExpan,
#                       cycIdx = list( all = cycIdx,
#                                      preg = pregCycIdx ),
#                       pregCycBool = pregCycBool,
#                       idPregIdx = idPregIdx,
#                       gamIsBinBool = gamIsBinBool,
#                       naSexIdx = naSexIdx,
#                       naSexPriors = getNaSexProbs(naSexIdx),
#                       n = n,
#                       nObs = length(idDayExpan),
#                       sexIdx = list(all=seq(1, length(idDayExpan))),
#                       q = q,
#                       subjId = unique(id[sexBool]),
#                       varNames = colnames(U) )
  
#   idx <- list( cyc = list( all = cycIdx,
#                            uPreg = pregUBool)
                           
               
  
  samplerObj <- list( U            = U,
                      idDayExpan   = idDayExpan,
                      gamBinBool   = gamIsBinBool,
                      varNames     = colnames(U),
                      useNaSexBool = useNaSexBool,
                      subjId       = subjId,
                      n            = n,
                      nObs         = nObs,
                      q            = q )
  
  if (!useNaSexBool) {
    samplerObj$idx <- idx
  }
  else {
#     samplerObj$naSexDat <- naSexDat
#     samplerObj$sexPriors <- sexPriors
    samplerObj$naSexDat <- naSexDat
  }

  
  return (samplerObj)
}




# Calculate prob of sex for missing obs ----------------------------------------

getNaSexProbs <- function(sexMissIdx) {
  rep(0.36, length(sexMissIdx))
}