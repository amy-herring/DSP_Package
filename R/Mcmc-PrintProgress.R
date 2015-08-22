
# Print burn progress bar ------------------------------------------------------

printBurn <- function(nStar, barLen) {
  prevStar <- nStar - 1
  
  cat(pasteC(rep("\b", barLen - prevStar + 1)), "*", 
      pasteC(rep(" ",  barLen - prevStar - 1)), "|", sep="")
}




# Print progress / verbose info ------------------------------------------------

printProg <- function(trackProg, nSamp, s, gamOut, varNames, gamIsBinBool, metCtr) {
  format4 <- function(x) format(round(x, 4), nsmall=4)
  perc <- function(x) formatC(round(100 * x), width=3)
  contBool <- (FALSE %in% gamIsBinBool)
  
  if ( identical(trackProg, "percent") )
    cat(round(100 * s / nSamp), "%..  ", sep="")
  
  else {
    meanGam <- sapply(gamOut[1:s, ], mean)
    quantGam <- t( apply(gamOut[1:s, ], MARGIN=2, FUN=quantile, probs=c(0.025, 0.975)) )
    spaceVec <- sapply(max(nchar(varNames)) - nchar(varNames) + 4,  function(x) 
      paste(rep(" ", x), collapse=""))
    headerSpace <- paste(rep(" ", 8 + max(nchar(varNames)), collapse=""))
    
    cat("\nCompletion percentage: ", perc(s / nSamp), "%\n", sep="")
    cat("phi acceptance rate:   ", perc(metCtr$phiAccept / s), "%\n", sep="")
    cat("Exponentiated coefficient statistics:\n\n")
    cat(headerSpace, "  Mean      2.5%     97.5%", if (contBool) "    Accept", "\n", 
        headerSpace, "------    ------    ------", if (contBool) "    ------", "\n", sep="") 
    for (i in 1:length(varNames))
      cat("    ", varNames[i], spaceVec[i], format4(meanGam[i]), "    ", format4(quantGam[i, 1]),
          "    ", format4(quantGam[i, 2]), 
          if (!gamIsBinBool[i]) paste0("    ", perc(metCtr$gamAccept[i] / metCtr$gamTotal[i]), "%"),
          "\n", sep="")
    cat(paste(c(rep("-", length(headerSpace) + 26 + ifelse(contBool, 10, 0)), "\n"), collapse=""))
  }
}