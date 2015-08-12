
# Reduce data sets to pertinant observations -----------------------------------

getRedDat <- function(cleanDat, varNames, idVec, cycList) {
  idName <- varNames$id
  cycName <- varNames$cyc
  baselineRed <- cycleRed <- NULL
  
  # Create reduced baseline data
  if (!is.null(cleanDat$bas)) {
    basObsBool <- (cleanDat$bas[[idName]] %in% idVec)
    baselineRed <- cleanDat$bas[basObsBool, varNames$basIncl, drop=FALSE]
  }
  
  # Create reduced cycle data
  if (!is.null(cleanDat$cyc)) {
    cycObsIdx <- getObsIdx(cleanDat$cyc[[idName]], cleanDat$cyc[[cycName]], idVec, cycList)
    cycleRed <- cleanDat$cyc[cycObsIdx, varNames$cycIncl, drop=FALSE]
  }
  
  # Create reduced daily data
  dayObsIdx <- getObsIdx(cleanDat$day[[idName]], cleanDat$day[[cycName]], idVec, cycList)
  dailyRed <- cleanDat$day[dayObsIdx, varNames$dayIncl]
  
  return ( list( bas = baselineRed, 
                 cyc = cycleRed, 
                 day = dailyRed ) )
}





# Obtain idx of common pairs of id's + cycles ----------------------------------

getObsIdx <- function(datId, datCyc, idVec, cycList) {
  n <- length(idVec)
  idIdx <- lapply(idVec, function(x) which(datId == x))
  
  getCycBool <- function(j) ( datCyc[idIdx[[j]]] %in% cycList[[j]] )
  obsIdx <- lapply(1:n, function(j) subset(idIdx[[j]], getCycBool(j)))
  
  return ( unlist(obsIdx) )
}
