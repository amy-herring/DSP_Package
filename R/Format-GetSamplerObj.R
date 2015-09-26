
#' DSP input data structure creation
#' 
#' \code{getSamplerObj} creates data structures used by the DSP sampler 
#' algorithm.  In fact, the only data created by \code{dspDat} that is used by
#' the sampler is the data returned by this function.
#' 
#' 
#' @inheritParams dspDat
#'   
#' @param modelObj The return object from \code{\link{getModelObj}}.
#'   
#' @param cycVec The return object from \code{\link{getCommonCyc}}.
#'   
#' @details Add details about the software
#'   
#' @return \code{getSamplerObj} returns a \code{list} containing the following 
#'   components:
#'   
#'   \describe{
#'   
#'   \item{\code{U}}{Design matrix of model covariates, after removing 
#'   observations for which intercourse did not occur (and for missing values 
#'   when \code{useNA} is specified as \code{"none"}).  See 
#'   \code{\link{getModelObj}} for details regarding construction of \code{U} 
#'   (before non-intercourse observations are removed). }
#'   
#'   \item{\code{idDayExpan}}{An \code{nObs}-\code{integer} vector that is used 
#'   to expand an \code{n}-vector to an \code{nObs}-vector. In more detail, it 
#'   is a vector of the form (1,1,1,2,2,2,2,...,n,n), where the number of 1's 
#'   corresponds to the number of study observations for which intercourse 
#'   occured for the first subject, the number of 2's corresponds to the number 
#'   of study observations for which intercourse occured for the second subject,
#'   etc. }
#'   
#'   \item{\code{gamBinBool}}{A \code{q}-\code{logical} vector such that the 
#'   \code{p}-th element states the logical status of whether the \code{q}-th 
#'   predictor variable in the design matrix assumes values of only either 0 or 
#'   1. }
#'   
#'   \item{\code{varNames}}{A \code{q}-\code{character} vector such that the 
#'   \code{p}-th element provides the name of the \code{q}-th predictor variable
#'   in the design matrix. }
#'   
#'   \item{\code{useNaSexBool}}{A \code{logical} value supplying the logical 
#'   status of whether \code{useNA} has value \code{"sex"}. }
#'   
#'   \item{\code{subjId}}{An \code{n}-vector of type determined at runtime 
#'   supplying the subject IDs. }
#'   
#'   \item{\code{n}}{An \code{integer} value supplying the number of subjects in
#'   the study (after removing subjects without any days of intercourse in the 
#'   fertile window). }
#'   
#'   \item{\code{nObs}}{An \code{integer} value supplying the number of 
#'   fertile-window study days during which intercourse occured. }
#'   
#'   \item{\code{q}}{An \code{integer} value supplying the number of predictor 
#'   variables in the design matrix. }
#'   
#'   \strong{idx}: if \code{useNA} is specified as \code{"none"}, then a 
#'   \code{list} with name \code{idx} is included as an element in the 
#'   \code{getSamplerObj} return object; if \code{useNA} is specified as 
#'   \code{"sex"}, then \code{idx} is not included in the \code{getSamplerObj} 
#'   return object, and instead has to be created each scan through the sampler.
#'   
#'   \describe{
#'   
#'   \item{\code{preg}}{An \code{integer} vector supplying the observation 
#'   indices (after removing non-intercourse days) such that the observation day
#'   was a day within a cycle for which a successful pregnancy occured. }
#'   
#'   \item{\code{U}}{A \code{list} containing elements \code{all} and 
#'   \code{preg}.
#'   
#'   \code{all} is a \code{list} containing with \code{q} elements.  Each 
#'   element is such that the \code{p}-th element is either \code{NULL} if the 
#'   \code{p}-th predictor variable takes on values not just 0 or 1, or 
#'   otherwise is an \code{integer} vector supplying the observation indices 
#'   (after removing non-intercourse days) such that the \code{p}-th predictor 
#'   variable for that observation has a value of 1.
#'   
#'   \code{preg} is constructed in the same way as \code{all}, with the 
#'   exception that the indices are further restricted to belong to observations
#'   such that the observation day was a day within a cycle for which a 
#'   successful pregnancy occured. }
#'   
#'   \item{\code{pregCyc}}{A \code{list} that partitions the observations 
#'   (indices) by cycle, restricted to observations such that the observation 
#'   day was a day within a cycle for which a successful pregnancy occured.  In 
#'   more detail, the process is conceptually equivalent to (i) subsetting the 
#'   data to observations such that the observation day was a day within a cycle
#'   for which a successful pregnancy occured, and (ii) partitioning this 
#'   subsetted data into vectors of indices such that each vector contains the 
#'   observation indices belonging to a single cycle.   }
#'   
#'   \item{\code{subj}}{An \code{integer} vector such that the elements are 
#'   values corresponding to the indices of subjects that became pregnant during
#'   the study. }
#'   
#'   }
#'   
#'   \strong{naSexDat}: if \code{useNA} is specified as \code{"none"}, then the 
#'   following object is included as in element in the return \code{list} 
#'   object:
#'   
#'   }


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
  idDayExpan <- as.integer( rep(1:length(idIdx), times=sapply(idIdx, length)) )
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
  
  # Construct idx or naSexDat object -------------------------------------------

  if (!identical(useNA, "sex")) {
    useNaSexBool <- FALSE
    
    idx <- list( preg    = which(pregDayBool),
                 U       = list( all  = uIdx,
                                 preg = pregUIdx ),
                 pregCyc = pregCycIdx,
                 subj    = list( obs  = idIdx,
                                 preg = idPregIdx ) )
  } # End indexing data structure creation
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
  } # End NA sex data structure creation
 
  # Construct samplerObj data structure ----------------------------------------
  
  samplerObj <- list( U            = U,
                      idDayExpan   = idDayExpan,
                      gamBinBool   = gamIsBinBool,
                      varNames     = colnames(U),
                      useNaSexBool = useNaSexBool,
                      subjId       = subjId,
                      n            = n,
                      nObs         = nObs,
                      q            = q )

  if (!useNaSexBool) 
    samplerObj$idx <- idx 
  else 
    samplerObj$naSexDat <- naSexDat

  return (samplerObj)
}




# Calculate prob of sex for missing obs ----------------------------------------

getNaSexProbs <- function(sexMissIdx) {
  rep(0.36, length(sexMissIdx))
}