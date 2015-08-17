
#' MCMC sampler for day-specific probabilities model
#' 
#' \code{dsp} is an MCMC sampler for the methodology proposed by Dunson and 
#' Stanford in \emph{Bayesian inferences on predictors of conception 
#' probabilities} (2005).  The function either writes the samples to file or 
#' returns them as data, as specified by the \code{saveToFile} flag.
#' 
#' @param dspDat An object of \code{\link[base]{class}} \code{\link{dspDat}}.
#'   
#' @param nSamp The number of scans for which to perform the sampler.  
#'   Includes possible burn-in samples (as specified by the \code{nBurn}
#'   parameter), so that specifying e.g. \code{nSamp = 15000} and
#'   \code{nBurn = 5000} results in 15,000 total sampler scans, of which
#'   5,000 are designated as burn-in scans and 10,000 are recorded for
#'   posterior inference.
#'   
#' @param hypGam Either \code{NULL} or a \code{list} containing hyperparameters 
#'   to be specified for the exponentiated model regression coefficients.  None,
#'   some, or all of the hyperparameters can be or need be specified.
#'   
#'   Each exponentiated regression coefficient has a prior defined in terms of 5
#'   hyperparameters. These hyperparameters are the (i) prior probability of the
#'   point mass state, the (ii) shape and (iii) rate of the gamma distribution 
#'   state, and the (iv) lower (v) and upper bounds of the gamma distribution 
#'   state.
#'   
#'   If not specified by the function input, then a default value is provided 
#'   for each of the hyperparameters.  These default parameters, correponding to
#'   their description in the preceeding paragraph, are (i) \code{0.5}, (ii) 
#'   \code{1}, (iii) \code{1}, (iv) \code{0}, and (v) \code{Inf}.
#'   
#'   Exponentiated regression coefficient hyperparameter specifications must be 
#'   provided as follows.  If the input to \code{hypGam} is \code{NULL}, then 
#'   every hyperparameter is taken to be the default value.  If some of the 
#'   hyperparameters are to be specified, then \code{hypGam} must be a 
#'   \code{list} containing a sub-hierarchy of \code{list}s; each of these 
#'   second-level \code{list}s must have the name of one of model design matrix 
#'   variables.  Thus if non-\code{NULL}, then \code{hypGam} is a \code{list} 
#'   containing between \code{1} and \code{q} \code{list}s, where \code{q} is 
#'   the number of covariates in the model (after recoding categorical variables
#'   to dummy-variable form).
#'   
#'   Each second-level \code{list} in \code{hypGam} must contain between 
#'   \code{1} and \code{5} \code{numeric} values with possible names (i) 
#'   \code{p}, (ii) \code{a}, (iii) \code{b}, (iv) \code{bndL}, or (v) 
#'   \code{bndU} corresponding to the hyperparameter description from before 
#'   with matching index.  The order of the objects in either level of 
#'   \code{hypGam} does not matter.
#'   
#' @param tuningGam *******
#'   
#' @param hypPhi Either (i) \code{NULL} or (ii) a \code{list} containing one or 
#'   two \code{numeric} objects with names \code{c1} and/or \code{c2}.  These 
#'   values correspond (respectively) to the shape and rate parameters of the 
#'   prior (gamma) distribution for the variance parameter of the woman-specific
#'   fecundability mulitpliers.
#'   
#' @param tuningPhi The value of the tuning parameter for the Metropolis step 
#'   for the variance parameter of the woman-specific fecundability mulitpliers.
#'   The proposal value for this variance parameter for the (s+1)-th scan in the
#'   sampler algorithm, is given by sampling from (the absolute value of) a 
#'   uniform distribution with lower/upper bound given by adding/subtracting the
#'   tuning parameter to/from the s-th value of the parameter.
#'   
#' @param trackProg **********
#'   
#' @param progQuants **********
#'   
#' @param saveToFile **********
#'   
#' @param nBurn ***********
#'   
#' @param nThin **********
#'   
#' @details ****
#'   
#' @return *****
#'   
#'   
#' @author David A. Pritchard and Sam Berchuck, 2015
#'   
#' @references Dunson, David B., and Joseph B. Stanford. "Bayesian inferences on
#'   predictors of conception probabilities." \emph{Biometrics} 61.1 (2005): 
#'   126-133.


dsp <- function(dspDat, nSamp=1e4, hypGam=NULL, tuningGam=NULL, hypPhi=NULL, tuningPhi=0.3, 
                trackProg="percent", progQuants=seq(0.1, 1.0, 0.1), saveToFile=FALSE,
                nBurn=0, nThin=1) {

  
  # TODO: check if valid input
  
  list2env(dspDat$samplerObj, envir=environment())
  rm(dspDat)
  
  # Objects related to burn-in and thinning
  nBurn <- as.integer(nBurn)
  burnPhaseBool <- !identical(nBurn, 0L)
  nThin <- as.integer(nThin)
  thinIsOneBool <- identical(nThin, 1L)
  
  # Add a backslash to 'outPath' if necessary
  #outPath <- format_outPath(outPath)
  
  # Objects for progress statistics
  nKeep <- nSamp - nBurn
  printProgBool <- !identical(trackProg, "none")
  if (printProgBool) {
    # 'trackVals': sampler iterations at which we print the percentage of progress
    trackVals <- sapply(progQuants, function(x) tail(which(1:nKeep <= x * nKeep), 1))
    # Write first line of 'trackProg' == "percent" option
    if (identical(trackProg, "percent"))
      cat("Progress:  ")
  }
  
  # Combine default hyperparameters with custom user input hyperparameters
  hypPhi <- getHypPhi(hypPhi)
  hypGam <- getHypGam(varNames, hypGam)
  tuningGam <- getTuningGam(q)
  gamIsTrunBool <- sapply(1:q, function(j) 
    with(!isTRUE(all.equal(c(bndL, bndU), c(0, Inf))), data=hypGam[[j]]))
  
  # Set initial values:  uses mean of prior dists for phi and gamma
  Wfull <- integer(nrow(U))
  phi <- hypPhi$c1 / hypPhi$c2
  gamCoef <- theta <- getGamInit(hypGam, gamIsTrunBool)
  uProdBeta <- drop( U %*% log(gamCoef) )
  xi <- rep(1, n)
  xiDay <- xi[idDayExpan]
  
  # Metropolis acceptance rate counters
  metCtr <- list( phiAccept = 0L,
                  gamAccept = integer(q),
                  gamTotal  = integer(q) )
  
  # Inititalize MCMC output files / objects --------------------------------------
  if (saveToFile) {
    write(varNames, file=paste0(outPath, "GAMMA.csv"), sep=",", ncolumns=q)
    write(subjId, file=paste0(outPath, "XI.csv"), sep=",", ncolumns=n)
    write("phi", file=paste0(outPath, "PHI.csv"), sep=",", ncolumns=1)
  }
  else {
    phiOut <- numeric(nKeep)
    xiOut <- setNames(data.frame(matrix(nrow=nKeep, ncol=n)), subjId)
    gamOut <- setNames(data.frame(matrix(nrow=nKeep, ncol=q)), varNames)
  }

    
  # Begin MCMC sampler ==========================================================
  
  # Subtract 'nBurn' from 's' to prevent work later checking for thinning
  for (s in (1 - nBurn):nKeep) {
    
    # Sample latent variable W
    W <- sampW(uProdBeta, xiDay, pregDayBool, pregCycIdx)
    # 'Wfull': W's for every sex day, even those that are always 0
    Wfull[pregDayBool] <- W
    
    # Sample regression coefficients gamma
    for (h in 1:q) {
      
      # Binary case has closed-form full conditional
      if (gamIsBinBool[h]) {
        uProdBetaNoH <- getUProdBetaNoH(uProdBeta, drop(U[, h]), gamCoef[h])
        gamCoef[h] <- sampGam(W, uProdBetaNoH, xiDay, 
                              hypGam[[h]], uBool[[h]], pregUBool[[h]], TRUE)
        uProdBeta <- getUProdBeta(uProdBetaNoH, drop(U[, h]), gamCoef[h])
      }
      
      # Continuous case requires sampling via Metropolis algorithm
      else {
        # 'uProdTheta': list with uProdBeta for point mass / continuous value of theta
        uProdTheta <- getUProdTheta(uProdBeta, drop(U[, h]), gamCoef[h], theta[h])
    
        # M is the state part of the mixture distribution
        Mbool <- sampM(Wfull, xiDay, uProdTheta, hypGam[[h]]$p)
        
        # Corresponds to the point mass part of the mixture distribution
        if (Mbool) {
          gamCoef[h] <- 1
          theta[h] <- with(rgamma(1, shape=a, rate=b), data=hypGam[[h]])
          uProdBeta <- uProdTheta$point
        }
        # Corresponds to the continuous part of the mixture distribution
        else {
          propTheta <- abs( runif(1, theta[h] - tuningGam[h], theta[h] + tuningGam[h]) )
          uProdTheta$prop <- uProdTheta$point + drop(U[, h] * log(propTheta))
          
          logR <- getGamLogR(Wfull, xiDay, uProdTheta, theta[h], propTheta, hypGam[[h]])
          if (log(runif(1)) < logR) {
            theta[h] <- propTheta
            uProdBeta <- uProdTheta$prop
            metCtr$gamAccept[h] <- metCtr$gamAccept[h] + 1L
          }
          else {
            uProdBeta <- uProdTheta$cont
          }
          gamCoef[h] <- theta[h]
          metCtr$gamTotal[h] <- metCtr$gamTotal[h] + 1L
        }
      }
    } # End gamma update
    
    # Sample woman-specific fecundability multiplier xi
    xi <- sampXi(W, uProdBeta, phi, idIdx, idPregIdx, pregCycIdx, n)
    xiDay <- xi[idDayExpan]
    
    # Metroplois step for phi, the variance parameter for xi
    phiProp <- sampPhiProp(phi, tuningPhi)
    phiLogR <- getPhiLogR(xi=xi, phiCurr=phi, phiProp=phiProp, hypPhi=hypPhi)
    if (log(runif(1)) < phiLogR) {
      phi <- phiProp
      if (!burnPhaseBool)
        metCtr$phiAccept <- metCtr$phiAccept + 1L
    }
    
    # Write samples to output
    if (burnPhaseBool) {
      if (identical(s, 0L))
        burnPhaseBool <- FALSE
    }
    else if (thinIsOneBool || identical(s %% nThin, 0L)) {
      if (saveToFile) {
        write(phi, file=paste0(outPath, "PHI.csv"), sep=",", ncolumns=1, append=TRUE)
        write(xi, file=paste0(outPath, "XI.csv"), sep=",", ncolumns=n, append=TRUE)
        write(gamCoef, file=paste0(outPath, "GAMMA.csv"), sep=",", ncolumns=q, append=TRUE)
      }
      else {
        phiOut[s] <- phi
        xiOut[s, ] <- xi
        gamOut[s, ] <- gamCoef
      }
    }
    
    # Print progress / verbose info
    if (printProgBool && (s %in% trackVals))
      printProg(trackProg, nKeep, s, gamOut, varNames, gamIsBinBool, metCtr)

  } # End DSP sampler ----------------------------------------------------------
  
  # Construct and return sampler output ----------------------------------------  
  outObj <- list( formula = formula,
                  hypGam = hypGam, 
                  tuningGam = tuningGam, 
                  hypPhi = hypPhi, 
                  tuningPhi = tuningPhi,
                  nSamp = nSamp,
                  nBurn = nBurn,
                  nThin = nThin )  
  if (!saveToFile) {
    outObj$phi <- phiOut
    outObj$xi  <- xiOut
    outObj$gam <- gamOut
  }
  
  return (outObj)
}




getUProdTheta <- function(uProdBeta, UH, gamCoefH, thetaH) {
  if (identical(gamCoefH, 1))
    uProdTheta <- list( point = uProdBeta,
                        cont  = uProdBeta + drop(UH * log(thetaH)) )
  else
    uProdTheta <- list( point = uProdBeta - drop(UH * log(thetaH)),
                        cont  = uProdBeta )
  return (uProdTheta)
}




# Combine user input gam tuning vals with default ------------------------------

getTuningGam <- function(q) {
  rep(0.1, q)
}




# Calculate U[, -h] %*% beta ---------------------------------------------------

getUProdBetaNoH <- function(uProdBeta, UH, gamH) {
  if (identical(gamH, 1))
    uProdBeta
  else
    uProdBeta - drop(UH * log(gamH))
}




# Calculate U %*% beta ---------------------------------------------------------

getUProdBeta <- function(uProdBetaNoH, UH, gamH) {
  if (identical(gamH, 1))
    uProdBetaNoH
  else
    uProdBetaNoH + drop(UH * log(gamH))
}



# Print progress / verbose info ------------------------------------------------

printProg <- function(trackProg, nKeep, s, gamOut, varNames, gamIsBinBool, metCtr) {
  format4 <- function(x) format(round(x, 4), nsmall=4)
  perc <- function(x) formatC(round(100 * x), width=3)
  contBool <- (FALSE %in% gamIsBinBool)
  
  if ( identical(trackProg, "percent") )
    cat(round(100 * s / nKeep), "%..  ", sep="")
  
  else {
    meanGam <- sapply(gamOut[1:s, ], mean)
    quantGam <- t( apply(gamOut[1:s, ], MARGIN=2, FUN=quantile, probs=c(0.025, 0.975)) )
    spaceVec <- sapply(max(nchar(varNames)) - nchar(varNames) + 4,  function(x) 
      paste(rep(" ", x), collapse=""))
    headerSpace <- paste(rep(" ", 8 + max(nchar(varNames)), collapse=""))
    
    cat("\nCompletion percentage: ", perc(s / nKeep), "%\n", sep="")
    cat("phi acceptance rate:   ", perc(metCtr$phiAccept / s), "%\n", sep="")
    cat("Exponentiated coefficient statistics:\n\n")
    cat(headerSpace, "  Mean      2.5%     97.5%", if (contBool) "    Accept", "\n", sep="")
    cat(headerSpace, "------    ------    ------", if (contBool) "    ------","\n", sep="")  
    for (i in 1:length(varNames))
      cat("    ", varNames[i], spaceVec[i], format4(meanGam[i]), "    ", format4(quantGam[i, 1]),
          "    ", format4(quantGam[i, 2]), 
          if (!gamIsBinBool[i]) paste0("    ", perc(metCtr$gamAccept[i] / metCtr$gamTotal[i]), "%"),
          "\n", sep="")
    cat(paste(c(rep("-", length(headerSpace) + 26 + ifelse(contBool, 10, 0)), "\n"), collapse=""))
  }
}




# Ensure proper format for 'outPath' -------------------------------------------

format_outPath <- function(path) {
  pathLen <- nchar(path)
  if ( !identical(substr(path, pathLen, pathLen), "/") )
    path <- paste0(path, "/")

  return (path)
}
