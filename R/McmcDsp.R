
#' MCMC sampler for day-specific probabilities of conception methodology
#' 
#' \code{dsp} is an MCMC sampler for the methodology proposed by Dunson and 
#' Stanford in \emph{Bayesian Inferences on Predictors of Conception 
#' Probabilities} (2005).
#' 
#' @param dspDat An object of \code{class} \code{\link{dspDat}}.
#'   
#' @param nSamp The number of post-burn-in scans for which to perform the 
#'   sampler.
#'   
#' @param nBurn Number of sampler scans included in the burn-in phase.
#'   
#' @param nThin Value such that during the post-burn-in phase, only every 
#'   \code{nThin}-th scan is recorded for use in posterior inference.  Default 
#'   of \code{1} corresponds to every scan being retained.
#'   
#' @param hypGam Either \code{NULL} or a \code{list} containing hyperparameters 
#'   to be specified for the exponentiated model regression coefficients.  None,
#'   some, or all of the hyperparameters can be/need be specified.
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
#' @param tuningGam Either \code{NULL} or a list containing one or more 
#'   \code{numeric} objects each with the name of one of the model design matrix
#'   variables; the values are used as tuning parameters for the Metropolis step
#'   for regression coeffients corresponding to continuous predictor variables. 
#'   Categorical variables do not have a Metropolis step and values provided for
#'   them are ignored.  If \code{NULL}, then default values are provided for 
#'   every regression coefficient (that corresponds to a continuous predictor 
#'   variable).  If some but not all of the tuning parameters for regression 
#'   coefficients corresponding to continuous predictor variables are provided, 
#'   then default values are provided for the remaining. The default tuning 
#'   value for any regression coefficient is \code{0.25}.
#'   
#' @param hypPhi Either \code{NULL} or a \code{list} containing one or two 
#'   \code{numeric} objects with names \code{c1} and/or \code{c2}.  If supplied,
#'   then these values correspond (respectively) to the shape and rate 
#'   parameters of the prior (gamma) distribution for the variance parameter of 
#'   the woman-specific fecundability multipliers.  If \code{NULL}, then default
#'   values are provided for both \code{c1} and \code{c2}, and if exactly one of
#'   either \code{c1} or \code{c2} are provided, then a default value is 
#'   provided for the other.  The default values for \code{c1} and \code{c2} are
#'   \code{1} and \code{1}, respectively.
#'   
#' @param tuningPhi Metropolis tuning parameter for the variance parameter of 
#'   the woman-specific fecundability mulitpliers. The proposal value for this 
#'   variance parameter is sampled from a uniform distribution with support as 
#'   determined by the tuning parameter.
#'   
#' @param trackProg One of either \code{"none"}, \code{"percent"}, or 
#'   \code{"verbose"}; partial matching is supported.
#'   
#' @param progQuants Vector with values in (0,1].  Ignored if \code{trackProg} 
#'   is specified as \code{"none"}.  If \code{trackProg} is one of 
#'   \code{"percent"} or \code{"verbose"} then the specified output is printed 
#'   whenever the post-burn-in progress of the sampler reaches one of the 
#'   quantiles specified by \code{progQuants}.
#'   
#' @param saveToFile \code{logical} specifying whether the samples from the 
#'   post-burn-in phase are to be either written to file or returned as 
#'   \code{data.frame} objects.  Note that in either case output characterizing 
#'   the model is returned by the sampler.
#'   
#' @param outPath String specifying the local path into which output files 
#'   containing the MCMC samples are to be placed.  Ignored if \code{saveToFile}
#'   is \code{FALSE}.
#'   
#'   
#' @details Takes preprocessed fertility data in the form of a 
#'   \code{\link{dspDat}} object and performs an MCMC sampling algorithm for the
#'   Dunson and Stanford day-specific probabilities of conception methodology.
#'   
#'   Selection of the covariates to include in the model is performed when
#'   creating the \code{dspDat} object.
#'   
#'   
#' @return \code{dsp} returns a list containing the following objects
#'   
#'   \describe{ \item{\code{formula}}{Model \code{formula}, as passed to the 
#'   \code{dsp} sampler through the input to the \code{dspDat} parameter.}
#'   
#'   \item{\code{hypGam}}{\code{list} containing a sub-hierarchy of 
#'   \code{list}s, each containing the hyperparameter values used for the 
#'   sampler for the regression coefficients.}
#'   
#'   \item{\code{tuningGam}}{********}
#'   
#'   \item{\code{hypPhi}}{Hyperparameters for the variance parameter of the 
#'   woman-specific fecundability multipliers.}
#'   
#'   \item{\code{tuningPhi}}{Metropolis tuning parameter used for sampling the 
#'   variance parameter of the woman-specific fecundability multipliers.}
#'   
#'   \item{\code{nSamp}}{Input to \code{nSamp} parameter.}
#'   
#'   \item{\code{nBurn}}{Input to \code{nBurn} parameter.}
#'   
#'   \item{\code{nThin}}{Input to \code{nThin} parameter.}
#'   
#'   \item{\code{outPath}}{Input to \code{outPath} parameter.}
#'   
#'   \item{\code{phi}}{If \code{saveToFile} specified as \code{FALSE}, then a 
#'   vector containing the post-burn-in samples for the variance parameter of 
#'   the woman-specific fecundability multipliers.}
#'   
#'   \item{\code{xi}}{If \code{saveToFile} specified as \code{FALSE}, then a 
#'   \code{data.frame} containing the post-burn-in samples for the 
#'   woman-specific fecundability multipliers.}
#'   
#'   \item{\code{gam}}{If \code{saveToFile} specified as \code{FALSE}, then a 
#'   \code{data.frame} containing the post-burn-in samples for the regression 
#'   coefficients.}
#'   
#'   }
#'   
#' @author David A. Pritchard and Sam Berchuck, 2015
#'   
#' @references Dunson, David B., and Joseph B. Stanford. "Bayesian inferences on
#'   predictors of conception probabilities." \emph{Biometrics} 61.1 (2005): 
#'   126-133.


dsp <- function(dspDat, nSamp=1e4, nBurn=0, nThin=1, hypGam=NULL, tuningGam=NULL, 
                hypPhi=NULL, tuningPhi=0.3, trackProg="percent", progQuants=seq(0.1, 1.0, 0.1), 
                saveToFile=FALSE, outPath) {

  # TODO: check if valid input
  
  # TODO: work out 'format_outPath' function
  # TODO: verbose print doesn't work for 'saveToFile == TRUE'
  # TODO: 'getHypGam' not yet written
  
  list2env(dspDat$samplerObj, envir=environment())
  rm(dspDat)
  
  # Objects related to burn-in and thinning
  nBurn <- as.integer(nBurn)
  burnPhaseBool <- !identical(nBurn, 0L)
  nThin <- as.integer(nThin)
  thinIsOneBool <- identical(nThin, 1L)
  
  # Add a backslash to 'outPath' if necessary.  How does this work for Windows OS?
  if (saveToFile)
    outPath <- format_outPath(outPath)
  
  # Objects for progress statistics
  printProgBool <- !identical(trackProg, "none")
  bbl <- 50  # short for "burn bar length"; choice of 50 is arbitrary
  burnQuants <- seq(1 / bbl, 1, 1 / bbl)
  if (printProgBool) {
    # 'trackVals': sampler iterations at which we print the percentage of progress
    trackVals <- list() 
    trackVals$prog <- sapply(progQuants, function(x) tail(which(1:nSamp <= x * nSamp), 1))
    trackVals$burn <- sapply(burnQuants, function(x) tail(which(1:nBurn <= x * nBurn), 1) - nBurn)
    
    if (burnPhaseBool) {
      # Write first line of burn-in progress bar
      cat("Burn-in progress:  |", pasteC(rep(" ", bbl)), "|", sep="")
    }
    else if (identical(trackProg, "percent")) {
      # Write first line of 'trackProg' == "percent" option
      cat("Sampler progress:  ")
    }
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
    phiOut <- numeric(nSamp)
    xiOut <- setNames(data.frame(matrix(nrow=nSamp, ncol=n)), subjId)
    gamOut <- setNames(data.frame(matrix(nrow=nSamp, ncol=q)), varNames)
  }

    
  # Begin MCMC sampler ==========================================================
  
  # Subtract 'nBurn' from 's' to prevent work later checking for thinning
  for (s in (1 - nBurn):nSamp) {
    
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
      if (identical(s, 0L)) {
        burnPhaseBool <- FALSE
        printBurn(which(s == trackVals$burn), bbl)
        if (identical(trackProg, "percent")) 
          cat("\nPost-burn-in progress:  ")
        else if (identical(trackProg, "verbose"))
          cat("\nEntering post-burn-in phase of sampler..\n")
      }
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
    if (printProgBool) {
      if (burnPhaseBool) {
        if (s %in% trackVals$burn) {
          printBurn(which(s == trackVals$burn), bbl)
        }
      }
      else if (s %in% trackVals$prog)
        printProg(trackProg, nSamp, s, gamOut, varNames, gamIsBinBool, metCtr)
    }

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
  else
    outObj$outPath <- outPath
  
  return (outObj)
}




# Print burn progress bar ------------------------------------------------------

printBurn <- function(nStar, barLen) {
  prevStar <- nStar - 1
  
  cat(pasteC(rep("\b", barLen - prevStar + 1)), "*", 
      pasteC(rep(" ",  barLen - prevStar - 1)), "|", sep="")
}




# Calculate U %*% beta for theta_0 and theta_1 ---------------------------------

getUProdTheta <- function(uProdBeta, UH, gamCoefH, thetaH) {
  if (identical(gamCoefH, 1))
    uProdTheta <- list( point = uProdBeta,
                        cont  = uProdBeta + drop(UH * log(thetaH)) )
  else
    uProdTheta <- list( point = uProdBeta - drop(UH * log(thetaH)),
                        cont  = uProdBeta )
  return (uProdTheta)
}




# Convenience function for collapsing a vector ---------------------------------
pasteC <- function(x) {
  paste(x, collapse="")
}




# Combine user input gam tuning vals with default ------------------------------

getTuningGam <- function(q) {
  rep(0.25, q)
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




# Ensure proper format for 'outPath' -------------------------------------------

format_outPath <- function(path) {
  ifelse(identical(substr(path, nchar(path), nchar(path)), "/"), path,  paste0(path, "/"))
}
