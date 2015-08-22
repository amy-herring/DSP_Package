
# Combine baseline, cycle, and daily -------------------------------------------

getCombDat <- function(formula, redDat, varNames, fwLen, cycList) {
  # Number of cycles for each subject (after removing subjects with 0 common cycles)
  niVec <- Filter(sapply(cycList, length), f=as.logical)
  numOf <- list( subj = length(niVec),
                 cyc  = sum(niVec) )
  
  # Expansion vectors to convert bas / cyc to daily dimension
  basExpan <- rep(1:numOf$subj, times=(niVec * fwLen))
  cycExpan <- rep(1:numOf$cyc, each=fwLen)
  
  # Convert FW indicator to FW day (assumes that days are consistently ordered across cycles)
  if (!is.null(varNames$fw))
    redDat$day[[varNames$fw]] <- as.factor( rep(1:fwLen, times=numOf$cyc) )
  
  getPregOrNull <- function(set) if (varNames$preg %in% names(set)) varNames$preg else NULL
  # Rows from top to bottom: id/cyc, preg, preg, day vars, bas vars, cyc vars
  combDat <- with(redDat, list( day[, c(varNames$id, varNames$cyc)],
                                cyc[cycExpan, getPregOrNull(cyc), drop=FALSE],
                                day[, getPregOrNull(day), drop=FALSE],
                                day[, c(varNames$sex, varNames$modelVars$day), drop=FALSE], 
                                bas[basExpan, varNames$modelVars$bas, drop=FALSE],
                                cyc[cycExpan, varNames$modelVars$cyc, drop=FALSE] ))
  
  # Combine datasets into a daily dataset that still contains factors
  return( data.frame( Filter(length, combDat) ) )
}




# Create model objects corresponding to Dunson paper ---------------------------

getModelObj <- function(combDat, formula, varNames, fwLen) {
  list( Y  = combDat[seq(from=fwLen, to=nrow(combDat), by=fwLen), varNames$preg],
        X  = combDat[[varNames$sex]],
        id = combDat[[varNames$id]],
        U  = model.matrix(formula, data=combDat) )
}