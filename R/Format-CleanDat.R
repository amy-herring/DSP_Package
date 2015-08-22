
# Data preperataion cleaning step ----------------------------------------------
#
# Remove non-fertile window days, cycles that have wrong number of days, and
# observations / cycles that have missing data in the model variables.

getCleanDat <- function(baseline, cycle, daily, varNames, fwLen) {
  # 'NULL' value for 'fwInd' indicates that daily data has already been restricted to FW
  if (!is.null(varNames$fw))
    dayFw <- daily[convToBool(daily[[varNames$fw]]), ]
  else
    dayFw <- daily
  
  # Rows with missing data (in the needed cols) for baseline, cycle, daily datasets
  noMissBool <- getNoMissBool(baseline, cycle, dayFw, varNames)

  # Finds in daily data FW cycles of the right length and w/o missing (in the needed cols)
  cycInDailyIdx <- getCycInDailyIdx(dayFw[[varNames$id]], dayFw[[varNames$cyc]])
  cycInDailyNoMissIdx <- Filter(unlist(cycInDailyIdx, recursive=FALSE), f=function(j) 
    ((length(j) == fwLen) & !(FALSE %in% noMissBool$day[j])))

  cleanDat <- list( bas = baseline[noMissBool$bas, ],
                    cyc = cycle[noMissBool$cyc, ],
                    day = dayFw[unlist(cycInDailyNoMissIdx), ] )

  return (cleanDat)
}




# Indices in daily data for each cycle -----------------------------------------

getCycInDailyIdx <- function(idVec, cycVec) {
  idIdx <- lapply(unique(idVec), function(x) which(idVec == x))
  
  getCycByIdIdx <- function(idx) {
    thisCyc <- cycVec[idx]
    lapply(unique(thisCyc), function(j) idx[thisCyc == j])
  }
  return ( lapply(idIdx, getCycByIdIdx) )
}




# Boolean for missing in bas, cyc, day datasets --------------------------------

getNoMissBool <- function(baseline, cycle, daily, varNames) {
  
  getBool <- function(dat, inclNames) {
    if (is.null(dat)) 
      return (NULL)
    else
      return ( apply(is.na(dat[, inclNames]), MARGIN=1, FUN=function(x) !(TRUE %in% x)) )
  }
  
  noMissBool <- list( bas = getBool(baseline, varNames$basIncl),
                      cyc = getBool(cycle, varNames$cycIncl),
                      day = getBool(daily, varNames$dayIncl) )
  
  return (noMissBool)
}
