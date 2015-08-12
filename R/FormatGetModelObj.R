
# Combine base, cyc, day data into desing matrix -------------------------------

getModelObj <- function(formula, redDat, varNames, fwLen, cycList) {
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
  
  # Combine datasets into a daily dataset that still contains factors
  covDat <- Filter( length, list( redDat$bas[basExpan, , drop=FALSE],
                                  redDat$cyc[cycExpan, , drop=FALSE],
                                  redDat$day[, , drop=FALSE] ) )
  # Covariate matrix converted to design matrix
  U <- model.matrix(formula, data=data.frame(covDat))
  
  # Convert Y to cycle format (as in Dunson and Stanford paper)
  if (varNames$preg %in% names(redDat$cyc))
    Y <- redDat$cyc[[varNames$preg]]
  else
    Y <- redDat$day[seq(from=fwLen, to=nrow(redDat$day), by=fwLen), varNames$preg]
  
  modelObj <- list( Y = Y,
                    X = redDat$day[[varNames$sex]],
                    id = redDat$day[[varNames$id]],
                    U = U )
  return (modelObj)
}