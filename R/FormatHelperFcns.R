
# Partition the model variables by dataset -------------------------------------

getVarNames <- function(formula, baseline, cycle, daily, idName, cycName, sexName, fwName) {
  explNames <- all.vars(formula)
  varInSetBool <- function(set) explNames %in% names(set)
  
  varNames <- list( id = idName,
                    cyc = cycName,
                    preg = explNames[1],
                    sex = sexName,
                    fw = fwName,
                    basIncl = c(idName, explNames[varInSetBool(baseline)]),
                    cycIncl = c(idName, cycName, explNames[varInSetBool(cycle)]),
                    dayIncl = c(idName, cycName, sexName, fwName, explNames[varInSetBool(daily)]) )
  return (varNames)
}




# Obtain vector of common id's that are needed ---------------------------------
#
# NULL[[idName]] == NULL    <-- behavior of NULL object for when bas / cyc are null
# unique(NULL) == NULL
# length(NULL) == 0

getCommonId <- function(cleanDat, idName) {
  basId <- cleanDat$bas[[idName]]
  cycId <- unique( cleanDat$cyc[[idName]] )
  dayId <- unique( cleanDat$day[[idName]] )
  
  intrAllowNull <- function(a, b) if (!is.null(b)) intersect(a, b) else a  
  idVec <- Reduce(f=intrAllowNull, x=list(dayId, cycId, basId))
  
  return (idVec)
}




# Obtain list of common cycles that are needed ---------------------------------

getCommonCyc <- function(cleanDat, varNames, idVec) {
  cycId <- cleanDat$cyc[[varNames$id]]
  dayId <- cleanDat$day[[varNames$id]]
  cycName <- varNames$cyc
  n <- length(idVec)
  
  if (is.null(cleanDat$cyc))
    cycList <- lapply(idVec, function(x) unique(cleanDat$day[x == dayId, cycName]))
  else {
    cycById <- list( cyc = lapply(idVec, function(x) cleanDat$cyc[x == cycId, cycName]),
                     day = lapply(idVec, function(x) cleanDat$day[x == dayId, cycName]) )
    cycList <- lapply(1:n, function(j) intersect(cycById$cyc[[j]], cycById$day[[j]]) )
  }
  return (cycList)
}




# Convert possible 0/1 or "yes"/"no" to boolean --------------------------------
#
# PRE: already checked that all values are logical or numeric or start with "y","Y","n", or "N"

convToBool <- function(x) {
  if (is.logical(x))
    return (x)
  else if (is.numeric(x))
    return (as.logical(x))
  else
    sapply(substr(x, 1, 1), function(l) if ((l == "n") | (l == "N")) FALSE else TRUE)
}
