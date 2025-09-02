# Author: Brody Erlandson
#
# This file contains the helper functions for the wrapper functions in
# FunCZIDM.R. 


getBasisX <- function(X, varyingCov, covWithVC, df = 4, degree=3,
                      basisFunc=bs) {
  # Adds nessary columns to X to include varying coefficients. That is,
  # the basis functions of the varying variable element-wise multiplied by each 
  # column of X seperately. This allows for a linear fit of X with coefficients,
  # then take the last df*lenght(varyingCov) columns of X. Each df columns are  
  # the basis functions for the corresponding column of X in varyingCov 
  # X: data matrix
  # varyingCov: vector of values that the coefficients will be varying with
  # covWithVC: vector of column numbers of X that will have varying coefficients

  # Basis functions for varying coefficients
  basis <- basisFunc(varyingCov, df=df, intercept=FALSE)
  interiorKnots <- attr(basis, "knots")
  boundaryKnots <- attr(basis, "Boundary.knots")
  basis <- cbind(1, basis)
  df <- df + 1

  i = 1
  colIdentifier <- c()
  newX <- matrix(NA, nrow(X), 0)
  for (col in covWithVC) {
    while (i < col) {
      # If the next column is not in covWithVC, then it does not have varying 
      # coefficients.
      newX <- cbind(newX, X[, i])
      colIdentifier <- c(colIdentifier, i)
      i <- i + 1
    }
    basisTimeByX <- basis*X[, col] # element-wise multiplication by column
    newX <- cbind(newX, basisTimeByX)
    colIdentifier <- c(colIdentifier, rep(col, df))

    i <- i + 1
  }
  
  return(list("Xvar"=newX, "XvartoXColMapping"=colIdentifier,
              "interiorKnots"=interiorKnots, "boundaryKnots"=boundaryKnots))
}

getNullInBetaCIProportion <- function(output, XvartoXColMapping, colMapping,
                                      varyingCov, catNames, interiorKnots, 
                                      boundaryKnots, basisFunc=splines::bs, 
                                      df=4, fileName = "") {
  # Extracts the null in CI statistics from the output of FunCZIDM
  # output: output of FunCZIDM
  # nullFileName: file name to save the null in CI statistics

  # getting the points to evaluate
  minVarCov <- min(varyingCov)
  maxVarCov <- max(varyingCov)
  testPoints <- seq(minVarCov,  maxVarCov, length.out = 250)

  # get the basis functions for the test points
  basis <- basisFunc(testPoints, df=df, intercept=FALSE,
                   knots=interiorKnots,
                   Boundary.knots=boundaryKnots)
  basis <- cbind(1, basis)
  numDf <- df + 1
  
  # get the zero functions for checking the null in CI
  zeroFunc <- testPoints*0

  # writing the header of the csv to the file
  cat("category,covariate,Proportion 0 in 95 CI\n", file=fileName)
  
  startTV <- match(2, XvartoXColMapping) # First instance of the first covariate
  numCov <- length(colMapping) # number of covariates
  numCat <- length(catNames) # number of categories
  for (cov in 2:numCov) {
    # if the next instance in the mapping is not the next covariate, then
    # the covariate does not have varying coefficients. So startTV == endTV.
    if (cov == XvartoXColMapping[startTV + 1]) {
      endTV <- startTV + numDf - 1
      betaVCFits <- getBetaVC(output$beta, c(startTV, endTV), basis)
    } else {
      endTV <- startTV
      betaVCFits <- output$beta[startTV, , , drop=FALSE]
    }

    # the 2.5th and 97.5th percentiles of the beta functions
    betaQuant <- apply(betaVCFits, c(1, 2), quantile, probs=c(0.025, 0.975))
    # checking if the zero function is in the 95% CI of the beta functions
    zeroInCI <- betaQuant[1,,] < zeroFunc & zeroFunc < betaQuant[2,,]
    # the proportion of function that have the zero function in the 95% CI
    propZeroInCI <- apply(zeroInCI, 2, mean)
    covariate <- colMapping[cov][[1]] # covariate name
    for (catToCheck in 1:numCat) {
      category <- catNames[catToCheck]
      cat(category, covariate, propZeroInCI[catToCheck], sep = ",", 
          file=fileName, append=TRUE)
      cat("\n", file=fileName, append=TRUE)
    }
    startTV <- endTV + 1 
  }
}

errorChecksFunCZIDM <- function(counts, covariates, ids, varyingCov, 
                                iter, burnIn, thin, adjustFreq, proposalCap,
                                ZIGrouped, returnBurnIn, printProgress, 
                                toReturn, betaInitial, rInitial, priors, 
                                proposalVars, covWithoutVC, df, degree, basisFunc, 
                                saveToFile, fileName, saveNullInBetaCI, 
                                nullInBetaCIFileName) {
  nullsAllowed <- c("toReturn", "betaInitial", "rInitial", "priors", 
                    "proposalVars", "covWithoutVC")

  nonNegInt <- c("counts", "ids", "iter", "burnIn", "thin", 
                 "adjustFreq", "covWithoutVC")
  for (name in nonNegInt) {
    if (name %in% nullsAllowed && is.null(get(name))) {
      next # Allow NULL values for these parameters
    }
    checkNonNegativeInteger(get(name), name)
  }

  singleNumeric <- c("iter", "burnIn", "thin", "adjustFreq", "proposalCap")
  for (name in singleNumeric) {
    checkSingleNumeric(get(name), name)
  }

  bool <- c("ZIGrouped", "returnBurnIn", "printProgress", "saveToFile", 
            "saveNullInBetaCI")
  for (name in bool) {
    checkTrueOrFalse(get(name), name)
  }
  
  matrix <- c("counts", "betaInitial", "rInitial")
  for (name in matrix) {
    if (name %in% nullsAllowed && is.null(get(name))) {
      next # Allow NULL values for these parameters
    }
    checkMatrix(get(name), name, numericOnly = TRUE)
  }

  vector <- c("covWithoutVC", "ids", "varyingCov")
  for (name in vector) {
    if (name %in% nullsAllowed && is.null(get(name))) {
      next # Allow NULL values for these parameters
    }
    checkVector(get(name), name, numericOnly = TRUE)
  }

  list <- c("priors", "proposalVars")
  for (name in list) {
    if (name %in% nullsAllowed && is.null(get(name))) {
      next # Allow NULL values for these parameters
    }
    checkList(get(name), name)
  }

  # checks that only have to be done once
  checkDataFrame(covariates, "covariates")
  checkDataFrameTypes(covariates, "covariates")
  checkSameNumRows(counts, covariates, "counts", "covariates")
  checkNonNegative(proposalCap, "proposalCap")
  if (!is.null(toReturn)) {
    checkVectorStrings(toReturn, c("eta", "c", "u", "r", "phi", "lambda", "nu", 
                                   "tau", "xi", "RA"),
                       "returnable parameters")
  }
  if (!is.null(priors)) {
    checkListNamesAndValues(priors, list("first beta sd" = "non-negative",
                                         "a" = "non-negative", 
                                         "b" = "non-negative",
                                         "alpha" = "non-negative",
                                         "beta" = "non-negative",
                                         "kappaShape" = "non-negative",
                                         "kappaRate" = "non-negative"),
                            "priors")
  }
  if (!is.null(proposalVars)) {
    checkListNamesAndValues(proposalVars, 
                            list("beta proposal sd" = "non-negative",
                            "r proposal sd" = "non-negative"),
                            "proposal variances")
  }
}

checkSingleNumeric <- function(x, name) {
  if (!is.numeric(x)) {
    stop(paste(name, "must be numeric."))
  } 
  if (length(x) != 1) {
    stop(paste(name, "must be a single numeric value."))
  }
}

checkLessThanOrEqualTo <- function(x, limit, name) {
  if (any(x > limit)) {
    stop(paste(name, "must not contain values greater than", limit, "."))
  }
}

checkNonNegative <- function(x, name) {
  if (!is.numeric(x)) {
    stop(paste(name, "must be numeric."))
  }
  if (any(x < 0)) {
    stop(paste(name, "must not contain negative values."))
  }
}

checkNonNegativeInteger <- function(x, name) {
  if (!is.numeric(x)) {
    stop(paste(name, "must be non-negative integers."))
  }
  if (any(x < 0)) {
    stop(paste(name, "must not contain negative values."))
  }
  if (any(x != floor(x))) {
    stop(paste(name, "must be non-negative integers."))
  }
}

checkTrueOrFalse <- function(x, name) {
  if (!is.logical(x)) {
    stop(paste(name, "must be TRUE or FALSE."))
  }
  if (length(x) != 1) {
    stop(paste(name, "must be a single TRUE or FALSE value."))
  }
}

checkList <- function(x, name) {
  if (!is.list(x)) {
    stop(paste(name, "must be a list or NULL."))
  }
  if (length(x) == 0) {
    stop(paste(name, "must not be empty."))
  }
}

checkVector <- function(x, name, numericOnly = FALSE) {
  if (!is.vector(x)) {
    stop(paste(name, "must be a vector."))
  }
  if (length(x) == 0) {
    stop(paste(name, "must not be empty."))
  }
  if (any(is.na(x))) {
    stop(paste(name, "must not contain NA values."))
  }
  if (numericOnly && !is.numeric(x)) {
    stop(paste(name, "must be a numeric vector."))
  }
}

checkMatrix <- function(x, name, numericOnly = FALSE) {
  if (!is.matrix(x)) {
    stop(paste(name, "must be a matrix."))
  }
  if (nrow(x) == 0 || ncol(x) == 0) {
    stop(paste(name, "must not be empty."))
  }
  if (any(is.na(x))) {
    stop(paste(name, "must not contain NA values."))
  }
  if (numericOnly && !is.numeric(x)) {
    stop(paste(name, "must be a numeric matrix."))
  }
}

checkDataFrame <- function(x, name) {
  if (!is.data.frame(x)) {
    stop(paste(name, "must be a data frame."))
  }
  if (nrow(x) == 0 || ncol(x) == 0) {
    stop(paste(name, "must not be empty."))
  }
  if (any(sapply(x, is.na))) {
    stop(paste(name, "must not contain NA values."))
  }
}

checkDataFrameTypes <- function(x, name) {
  if (any(!sapply(x, function(col) is.numeric(col) || is.factor(col)))) {
    stop(paste(name, "must contain only numeric or factor columns."))
  }
} 

checkSameNumRows <- function(x, y, nameX, nameY) {
  if (nrow(x) != nrow(y)) {
    stop(paste(nameX, "and", nameY, "must have the same number of rows."))
  }
}

checkLength <- function(x, expectedLength, xName) {
  if (length(x) != expectedLength) {
    stop(paste(xName, "must have length", expectedLength, "."))
  }
}

checkValuesBetween <- function(x, min, max, name) {
  if (!is.numeric(x)) {
    stop(paste(name, "must be numeric."))
  }
  if (any(x <= min | x >= max)) {
    stop(paste(name, "must be between", min, "and", max, "."))
  }
}

checkVectorStrings <- function(x, expectedStrings, xName) {
  if (!is.character(x)) {
    stop(paste(xName, "must be a character vector."))
  }
  if (length(x) == 0) {
    stop(paste(xName, "must not be empty."))
  }
  if (any(!x %in% expectedStrings)) {
    stop(paste("The following", xName, "are not expected:",
                paste(x[!x %in% expectedStrings], collapse = ", ")))
  }
}

checkListNamesAndValues <- function(x, expected, xName) {
  for (name in names(x)) {
    outString <- paste0("The '", name, "' in ", xName)
    if (expected[[name]] == "numeric") {
      checkSingleNumeric(x[[name]], outString)
    } else if (expected[[name]] == "non-negative integer") {
      checkNonNegativeInteger(x[[name]], outString)
    } else if (expected[[name]] == "non-negative") {
      checkNonNegative(x[[name]], outString)
    } else if (expected[[name]] == NULL) {
      stop(paste(name, "is not a valid name in", xName, "."))
    }
  }
}

makeBlockMatrix <- function(X, RAND_EFF_COLS, UNIQUE_OBS_ID) {
  # Make block matrix for random effects
  # X: data matrix
  # RAND_EFF_COLS: column numbers of X that are random effects
  # UNIQUE_OBS_ID: vector of unique observation IDs (1-n, sequentially)

  NUM_OBS <- length(unique(UNIQUE_OBS_ID))
  NUM_RAND_EFF <- length(RAND_EFF_COLS)
  NUM_TOTAL <- nrow(X)

  blockMatrix <- matrix(0, NUM_TOTAL, NUM_OBS*NUM_RAND_EFF)
  observations <- unique(UNIQUE_OBS_ID)

  for (obsID in observations) {
    obsRows <- which(UNIQUE_OBS_ID == obsID)
    blockMatrix[obsRows, 
                ((obsID-1)*NUM_RAND_EFF+1):(obsID*NUM_RAND_EFF)] <- X[obsRows,
                                                                 RAND_EFF_COLS]
  }

  return(blockMatrix)
}

getCenterScaledCovProfile <- function(covVals, centerScaleList) {
  for (i in 1:length(centerScaleList)) {
    if (is.character(centerScaleList[[i]]) || is.null(centerScaleList[[i]])) {
      covVals[i] <- ifelse(covVals[i] == 1, 1, 0)
    } else {
      covVals[i] <- round((covVals[i] 
                   - centerScaleList[[i]]$center)/centerScaleList[[i]]$scale, 4)
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