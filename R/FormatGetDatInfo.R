
# Obtain info about the data munging process -----------------------------------

getDatInfo <- function(formula, baseline, cycle, daily, cleanDat, 
                       redDat, modelObj, varNames, fwLen, idVec, cycList) {
  
  idName <- varNames$id
  cycName <- varNames$cyc
  numCycInDaily <- list( raw = sum( sapply(getCycInDailyIdx(daily[[idName]], 
                                                            daily[[cycName]]), length) ),
                         cln = sum( sapply(getCycInDailyIdx(cleanDat$day[[idName]], 
                                                            cleanDat$day[[cycName]]), length) ) )
  
  # Number of subjects / cycles / days in raw data
  numRaw <- list( bas = list( sub = nrow(baseline) ),
                  cyc = list( sub = length(unique(cycle[[idName]])),
                              cyc = nrow(cycle) ),
                  day = list( sub = length(unique(daily[[idName]])),
                              cyc = numCycInDaily$raw,
                              day = nrow(daily) ) )    
  
  # Number of subjects / cycles / days in clean data (i.e. after removing missing)
  numClean <- list( bas = list( sub = nrow(cleanDat$bas) ),
                    cyc = list( sub = length(unique(cleanDat$cyc[[idName]])),
                                cyc = nrow(cleanDat$cyc) ),
                    day = list( sub = length(unique(cleanDat$day[[idName]])),
                                cyc = numCycInDaily$cln,
                                day = nrow(cleanDat$day) ) )

  # Number of subjects / cycles / days in clean data (i.e. after reducing to common subj + cycs)
  numRed <- list( sub = length(idVec),
                  cyc = length(unlist(cycList)),
                  day = nrow(redDat$day) )
  
  return ( list( numRaw = numRaw,
                 numClean = numClean,
                 numRed = numRed,
                 modelVars = all.vars(formula)[-1],
                 designMatVars = colnames(modelObj$U),
                 numSex = sum(convToBool(modelObj$X)),
                 numPreg= sum(convToBool(modelObj$Y)) ) )
}