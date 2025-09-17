library(FunCZIDM)

getCenterScaledCovProfile <- function(covVals, centerScaleList) {
  for (i in 1:length(centerScaleList)) {
    if (is.character(centerScaleList[[i]]) || is.null(centerScaleList[[i]])) {
      covVals[i] <- ifelse(covVals[i] == 1, 1, 0)
    } else {
      covVals[i] <- round((covVals[i] - 
                     centerScaleList[[i]]$center)/centerScaleList[[i]]$scale, 4)
    }
  }

  return(covVals)
}

getCenterScaledChange <- function(change, centerScaleList, covIdx, covVals) {
  if (is.character(centerScaleList[[covIdx]]) 
      || is.null(centerScaleList[[covIdx]])) {
    # If the covariate is a categorical variable, then the change is not scaled.
    change <- ifelse(covVals[covIdx] == 1, -1, 1)
    return(change) 
  } else {
    # If the covariate is a continuous variable, then scale the change.
    return(change / centerScaleList[[covIdx]]$scale)
  }
}

RAData <- function(output, covariates, xAxisRange, interceptBetas = TRUE) {
  minTimepoint <- min(output$varyingCov)
  maxTimepoint <- max(output$varyingCov)
  if (minTimepoint < xAxisRange[1]) {
    minTimepoint <- xAxisRange[1]
  }
  if (maxTimepoint > xAxisRange[2]) {
    maxTimepoint <- xAxisRange[2]
  }
  testPoints <- seq(minTimepoint, maxTimepoint, length.out = 100)
  basisTP <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                              knots=output$interiorKnots,
                              Boundary.knots=output$boundaryKnots)
  basisTP <- cbind(1, basisTP)
  numDf <- output$df + 1
  
  RA <- FunCZIDM:::getRAFits(output$beta, basisTP, covariates,
                             output$XvartoXColMapping)

  if (interceptBetas) {
    intercepts <- FunCZIDM:::getBetaVC(output$beta, c(1, numDf), basisTP) 
    return(list("RA" = RA, "intercepts" = intercepts, "points" = testPoints))
  }
  return(list("RA" = RA, "points" = testPoints))
}

deltaRAData <- function(output, cov, catsToCheckNames, covariates, changes,
                        xAxisRange, betas = TRUE) {
  minTimepoint <- min(output$varyingCov)
  maxTimepoint <- max(output$varyingCov)
  if (minTimepoint < xAxisRange[1]) {
    minTimepoint <- xAxisRange[1]
  }
  if (maxTimepoint > xAxisRange[2]) {
    maxTimepoint <- xAxisRange[2]
  }
  testPoints <- seq(minTimepoint, maxTimepoint, length.out = 100)
  basisTP <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                              knots=output$interiorKnots,
                              Boundary.knots=output$boundaryKnots)
  basisTP <- cbind(1, basisTP)
  numDf <- output$df + 1

  catsToCheck <- c()
  for (cat in catsToCheckNames) {
    catsToCheck <- c(catsToCheck, which(grepl(cat, output$catNames)))
  }

  fit <- FunCZIDM:::getFit(output$beta, basisTP, covariates,
                           output$XvartoXColMapping)
  sumExpFit <- FunCZIDM:::getSumExpFit(fit)

  cov <- which(cov == output$colMapping)
  startTV <- match(cov, output$XvartoXColMapping) # First instance of the cov
  if (cov == output$XvartoXColMapping[startTV + 1]) {
    # If the next instance in the mapping is not the next covariate, then
    # the covariate does not have varying coefficients. So startTV == endTV.
    endTV <- startTV + numDf - 1
    betaVCFits <- FunCZIDM:::getBetaVC(output$beta, c(startTV, endTV), basisTP)
  } else {
    betaVCFits <- output$beta[startTV, , , drop=FALSE]
    betaVCFits <- array(rep(betaVCFits, each = 100), 
                        dim = c(100, dim(betaVCFits)[2], dim(betaVCFits)[3]))
  }

  result <- list()
  if (length(changes) == 1) {
    for (catToCheck in catsToCheck) {
      catName <- output$catNames[catToCheck]
      RAMeanAndCI <- FunCZIDM:::getDeltaRAMeanAndCI(catToCheck - 1,# zero-idxing
                                                    betaVCFits, fit, sumExpFit,
                                                    change=changes)
      raMat <- cbind(
        lowerCI = RAMeanAndCI[, 1],
        upperCI = RAMeanAndCI[, 3],
        mean = RAMeanAndCI[, 2],
        testpoints = testPoints
      )
    
      entry <- list(RA = data.frame(raMat))
    
      if (betas) {
        betaConfInt <- apply(betaVCFits[,catToCheck,], 1,
                             function(x) quantile(x, c(0.025, 0.975)))
        betaMat   <- cbind(
          lowerCI = betaConfInt[1,],
          upperCI = betaConfInt[2,],
          mean = rowMeans(betaVCFits[,catToCheck,]),
          testpoints = testPoints
        )
        entry$beta <- data.frame(betaMat)
      }

      result[[catName]] <- entry
    }
  } else {
    for (catToCheck in catsToCheck) {
      firstPass = TRUE
      for (change in changes) {
        catName <- output$catNames[catToCheck]
        RAMeanAndCI <- FunCZIDM:::getDeltaRAMeanAndCI(catToCheck - 1,# zero-idx
                                                      betaVCFits, fit,
                                                      sumExpFit, change=change)
        if (firstPass) {
          nullInCI = (RAMeanAndCI[, 1] <= 1)*(1 <= RAMeanAndCI[, 3])
          mean = RAMeanAndCI[, 2] 
          firstPass = FALSE
        } else {
          nullInCI = cbind(nullInCI, 
                           (RAMeanAndCI[, 1] <= 1)*(1 <= RAMeanAndCI[, 3]))
          mean = cbind(mean, RAMeanAndCI[, 2]) 
        }
      }
      raMat <- list(
        nullInCI = t(nullInCI),
        means = t(mean),
        testpoints = testPoints
      )

      result[[catName]] <- raMat
    }
  }

  return(result)
}

diversityData <- function(output, covariates, xAxisRange, l = 0) {
  minTimepoint <- min(output$varyingCov)
  maxTimepoint <- max(output$varyingCov)
  if (minTimepoint < xAxisRange[1]) {
    minTimepoint <- xAxisRange[1]
  }
  if (maxTimepoint > xAxisRange[2]) {
    maxTimepoint <- xAxisRange[2]
  }
  testPoints <- seq(minTimepoint, maxTimepoint, length.out = 100)
  basisTP <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                              knots=output$interiorKnots,
                              Boundary.knots=output$boundaryKnots)
  basisTP <- cbind(1, basisTP)
  numDf <- output$df + 1
  
  div <- FunCZIDM:::getHillDiversity(output$beta, basisTP, covariates, 
                                     output$XvartoXColMapping, l =l)

  return(list("div" = div, "testpoints" = testPoints))
}

deltaDivData <- function(output, cov, covariates, changes, l, xAxisRange,
                         betas = TRUE) {
  minTimepoint <- min(output$varyingCov)
  maxTimepoint <- max(output$varyingCov)
  if (minTimepoint < xAxisRange[1]) {
    minTimepoint <- xAxisRange[1]
  }
  if (maxTimepoint > xAxisRange[2]) {
    maxTimepoint <- xAxisRange[2]
  }
  testPoints <- seq(minTimepoint, maxTimepoint, length.out = 100)
  basisTP <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                              knots=output$interiorKnots,
                              Boundary.knots=output$boundaryKnots)
  basisTP <- cbind(1, basisTP)
  numDf <- output$df + 1

  fit <- FunCZIDM:::getFit(output$beta, basisTP, covariates, 
                           output$XvartoXColMapping)
  sumExpFit <- FunCZIDM:::getSumExpFit(fit)

  cov <- which(cov == output$colMapping)
  startTV <- match(cov, output$XvartoXColMapping) # First instance of the cov
  if (cov == output$XvartoXColMapping[startTV + 1]) {
    # If the next instance in the mapping is not the next covariate, then
    # the covariate does not have varying coefficients. So startTV == endTV.
    endTV <- startTV + numDf - 1
    betaVCFits <- FunCZIDM:::getBetaVC(output$beta, c(startTV, endTV), basisTP)
  } else {
    betaVCFits <- output$beta[startTV, , , drop=FALSE]
    betaVCFits <- array(rep(betaVCFits, each = 100), 
                        dim = c(100, dim(betaVCFits)[2], dim(betaVCFits)[3]))
  }

  result <- list()
  if (length(changes) == 1) { 
      divMeanAndCI <- FunCZIDM:::getDeltaHillDivMeanAndCI(betaVCFits, fit, 
                                                          sumExpFit, 
                                                          change=changes, l=l)
      divMat <- cbind(
        lowerCI = divMeanAndCI[, 1],
        upperCI = divMeanAndCI[, 3],
        mean = divMeanAndCI[, 2],
        testpoints = testPoints
      )
    
      result <- list(div = divMat)
  } else {
    firstPass = TRUE
    for (change in changes) {
      divMeanAndCI <- FunCZIDM:::getDeltaHillDivMeanAndCI(betaVCFits, fit,
                                                          sumExpFit, 
                                                          change=change, l=l)
      if (firstPass) {
        nullInCI = (divMeanAndCI[, 1] <= 1)*(1 <= divMeanAndCI[, 3])
        mean = divMeanAndCI[, 2] 
        firstPass = FALSE
      } else {
        nullInCI = cbind(nullInCI, 
                         (divMeanAndCI[, 1] <= 1)*(1 <= divMeanAndCI[, 3]))
        mean = cbind(mean, divMeanAndCI[, 2]) 
      }
    }
    divMat <- list(
      nullInCI = t(nullInCI),
      means = t(mean),
      testpoints = testPoints
    )
    
    result <- list(div = divMat)
  }

  return(result)
}