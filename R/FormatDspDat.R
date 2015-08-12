
# Mung data into format for use by mcmc sampler ================================

dspDat <- function(formula, baseline=NULL, cycle=NULL, daily, 
                   idName, cycName, sexName, fwName=NULL, fwLen) {
  source("Format_Data/FormatCheckValidInput.R")
  source("Format_Data/FormatHelperFcns.R")
  source("Format_Data/FormatCleanDat.R")
  source("Format_Data/FormatGetRedDat.R")
  source("Format_Data/FormatGetModelObj.R")
  source("Format_Data/FormatGetSamplerObj.R")
  source("Format_Data/FormatGetDatInfo.R")
  
  # TODO: check valid input
  
  # Partition the model variable names by dataset
  varNames <- getVarNames(formula, baseline, cycle, daily, idName, cycName, sexName, fwName)
  
  # Sort data by id/cyc
  if (!is.null(baseline)) baseline <- baseline[order(baseline[[idName]]), ]
  if (!is.null(cycle)) cycle <- cycle[order(cycle[[idName]], cycle[[cycName]]), ]
  daily <- daily[order(daily[[idName]], daily[[cycName]]), ]
  
  # For daily data: remove non-FW days, cycles that have wrong number of FW days or include
  # missing in the cycle (in the model variables).  For baseline / cycle: remove observations 
  # that have missing data (in the model variables).
  cleanDat <- getCleanDat(baseline, cycle, daily, varNames, fwLen)
  
  # Reduce data to subjects and cycles that are common to all datasets (after possibly removing
  # missing)
  idVec <- getCommonId(cleanDat, idName)
  cycList <- getCommonCyc(cleanDat, varNames, idVec)
  redDat <- getRedDat(cleanDat, varNames, idVec, cycList)

  # Create X, Y, and U (from the Dunson and Stanford paper)
  modelObj <- getModelObj(formula, redDat, varNames, fwLen, cycList)
  
  # Create objects for use in MCMC sampler (see 'FormatGetSamplerObj.R' for more details)
  samplerObj <- getSamplerObj(modelObj, fwLen)
  
  # Stats related to munging process for use by summary fcn
  datInfo <- getDatInfo(formula, baseline, cycle, daily, cleanDat, 
                        redDat, modelObj, varNames, fwLen, idVec, cycList)
  
  # Construct dspDat object
  dspDat <- list( baseline = baseline,
                  cycle = cycle,
                  daily = daily,
                  cleanDat = cleanDat,
                  redDat = redDat,
                  modelObj = modelObj,
                  samplerObj = samplerObj,
                  datInfo = datInfo )
  structure(dspDat, class="dspDat")
}




# Describes the data munging process ===========================================

summary.dspDat <- function(dspDat) {
  datInfo <- dspDat$datInfo
  hline <- paste0(rep("-", 60), collapse="")
  
  # internal functions to assist printing --------------------------------------
  
  printStats <- function(title, dat) {
    cat("\n", hline, "\n", title, ":\n\n",
        "    baseline data:  ", dat$bas$sub, " subjects\n",
        "    cycle data:     ", dat$cyc$sub, " subjects with ", dat$cyc$cyc, " cycles\n",
        "    daily data:     ", dat$day$sub, " subjects with ", dat$day$cyc, " cycles and ", 
        dat$day$day, " days\n", sep="")
  }
  
  printVars <- function(charVec, returnWidth=40) {
    if (is.null(charVec))
      return (NULL)
    
    # else
    currLen <- 0
    outVec <- NULL
    
    for (i in 1:length(charVec)) {
      
      if (currLen >= returnWidth) {
        outVec <- paste0(outVec, "\n", paste(rep("", 22), collapse=" "), charVec[i], ", ")
        currLen <- 0
      }
      else {
        outVec <- paste0(outVec, charVec[i], ", ")
        currLen <- currLen + nchar(charVec[i])
      }
    }
    
    return ( substr(outVec, start=1, stop=(nchar(outVec) - 2)) )
  }
  
  # ----------------------------------------------------------------------------
  
  printStats("Raw data", datInfo$numRaw)
  printStats("Clean data", datInfo$numClean)

  numRed <- datInfo$numRed
  cat("\n", hline, "\nCombining the data:\n\n",
      "    common observations:  ", numRed$sub, 
      " subjects with ", numRed$cyc, 
      " cycles and ", numRed$day, " days\n", sep="")
  
  cat("\n", hline, "\nThe model variables:\n\n",
      "    variable names:  ", printVars(datInfo$modelVars), "\n",
      "    design matrix:   ", printVars(datInfo$designMatVars), "\n\n", sep="")
  
  # TODO: ave number of cycles in study (tot, preg, not preg), num pregnant, num sex
  
  cat(hline, "\n", sep="")
}


