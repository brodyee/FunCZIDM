# Author: Brody Erlandson
#
# This file contains the helper functions and wrapper functions for the 
# simulation studies.

###############################################################################
#                   Helper functions for plot and results 
############################################################################### 

getTraceplot <- function(sample, trueValue=NA) {
  plot(sample, type="l", xlab="", ylab="")

  if (!is.na(trueValue)) {
    abline(h=trueValue, col="red")
  }
}

saveTraceplotGrid <- function(samples, filename, width=1000, height=1000, 
                              trueValues=NA) {
  if (length(dim(samples)) == 3) {
    totalSamples <- dim(samples)[1]*dim(samples)[2]
    dims <- dim(samples)[1:2]
    doubleLoop <- TRUE
  } else if (length(dim(samples)) == 2) {
    totalSamples <- dim(samples)[1]
    dims <- dim(samples)[1]
    doubleLoop <- FALSE
  }
  NUM_COL <- 2
  NUM_ROW <- ceiling(totalSamples/NUM_COL)
  height <- height*ceiling((NUM_ROW/12))
  
  filename <- paste0(substr(filename, 1, nchar(filename)-3), "png")
  png(filename, width=width, height=height)

  par(mfrow=c(NUM_ROW, NUM_COL), mar=c(0, 2, 0, 0), oma=c(0,0,0,0), 
      mgp=c(0, 0, -.5), tck=-.01)
  if (doubleLoop) {
    for (j in 1:dims[2]) {
      for (i in 1:dims[1]) {
        if (length(trueValues) > 1) {
          getTraceplot(samples[i, j, ], trueValues[i, j])
        } else {
          getTraceplot(samples[i, j, ])
        }
      }
    }
  } else {
    for (i in 1:dims) {
      if (length(trueValues) > 1) {
        getTraceplot(samples[i, ], trueValues[i])
      } else {
        getTraceplot(samples[i, ])
      }
    }
  }

  dev.off()
}

getCaptureRate <- function(sample, trueValues, credLevel=.95, 
                           checkZeros=FALSE) {
  NUM_ROW <- dim(trueValues)[1]
  NUM_COL <- dim(trueValues)[2]
  captured <- 0
  credProbs <- c((1-credLevel)/2, 1-(1-credLevel)/2)

  nearZero <- 0
  for (i in 1:NUM_ROW) {
    for (j in 1:NUM_COL) {
      credIterval <- quantile(sample[i, j, ], probs = credProbs)
      if (credIterval[1] <= trueValues[i, j] 
          & trueValues[i, j] <= credIterval[2]) {
        captured <- captured + 1
      } else if (checkZeros) {
         if (trueValues[i, j] < .001) {
           nearZero <- nearZero + 1
           if (round(credIterval[1], 2) == 0) {
             captured <- captured + 1
           }
         }
      } 
    }
  }

  cat("------------------- ", credLevel*100,
      "% Credible Interval Capture Rate ------------------------\n",
       sep = "")
  cat("Capture rate: ", captured/(NUM_ROW*NUM_COL), "\n", sep = "")
  if (checkZeros) {
    cat("Number of values near zero that weren't captured: ", nearZero,
        "\n", sep = "")
    cat("Proprtion of total values near zero that weren't captured: ", 
        (nearZero)/(NUM_ROW*NUM_COL), "\n", sep = "")
  }

}

getMAE <- function(sample, trueValues) {
  NUM_ROW <- dim(trueValues)[1]
  NUM_COL <- dim(trueValues)[2]
  mae <- 0

  for (i in 1:NUM_ROW) {
    for (j in 1:NUM_COL) {
      mae <- mae + abs(trueValues[i, j] - mean(sample[i, j, ]))
    }
  }

  cat("----------------------- Mean Absolute Error -------------------------\n")
  cat("MAE: ", mae/(NUM_ROW*NUM_COL), "\n", sep = "")
}

getAitchison <- function(sample, trueValues) {
  NUM_ROW <- dim(trueValues)[1]
  NUM_COL <- dim(trueValues)[2]
  mad <- 0
  eps <- .0000001

  for (i in 1:NUM_ROW) {
    mad <- mad + FunCZIDM:::aitchison_distance(trueValues[i,] + eps,
                                               rowMeans(sample[i,,]) + eps)
  }

  cat("-------------------- Mean Aitchison Distance ------------------------\n")
  cat("MAD: ", mad/NUM_ROW, "\n", sep = "")
}

###############################################################################
#                   Helper functions for alternative models 
############################################################################### 

# This file contains the helper functions for the wrapper functions in
# FunCZIDM.R. 
# Author: Brody Erlandson

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

FunCDM <- function(counts, covariates, ids, varyingCov,
                   iter = 10000, burnIn = 5000, thin = 1, adjustFreq = 250,
                   proposalCap = 0, returnBurnIn = FALSE, printProgress = TRUE, 
                   toReturn = NULL, betaInitial = NULL, rInitial = NULL,
                   priors = NULL, proposalVars = NULL, covWithVC=NULL, df=4, 
                   degree=3, basisFunc=splines::bs, saveToFile = TRUE, 
                   fileName = "output.rds", saveNullInBetaCI = TRUE, 
                   nullInBetaCIFileName = "nullInBetaCI.csv") {               
  # error checks
  # No error check made for this, as this function is only used for simulation 
  # studies

  adjustProp <- adjustFreq != 0 # If adjustFreq is 0, no adjustment of proposals
  capProposals <- proposalCap > 0 # If proposalCap is 0, no capping of proposals

  # processing the data
  centerScaleList <- list() # List to store the center and scale of each cov
  centerScaleList[[1]] <- "no c/s" # Intercept does not have c/s
  listIdx <- 2 # Start from 2 to skip intercept
  for (i in 1:ncol(covariates)) {
    if (is.numeric(covariates[, i])) {
      centerScaleList[[listIdx]] <- list(center = mean(covariates[, i]),
                                       scale = sd(covariates[, i]))
      covariates[, i] <- (covariates[, i] - centerScaleList[[listIdx]]$center) / 
                          centerScaleList[[listIdx]]$scale
    } else {
      centerScaleList[[listIdx]] <- "no c/s" # No center/scale for factors
      if (length(levels(covariates[, i])) > 2) {
        for (j in 3:length(levels(covariates[, i]))) {
          # If the factor has more than 2 levels, then we need to account for 
          # the extra dummy variables for the factor levels
          listIdx <- listIdx + 1
          centerScaleList[[listIdx]] <- "no c/s"
        }
      } 
    }
    listIdx <- listIdx + 1
  }
  X <- model.matrix(~., data=covariates) 
  
  colMapping <- list() # List mapping the covariate numbers to colnames
  i <- 1 
  for (name in colnames(X)) {
    colMapping[[i]] <- name
    i <- i + 1
  }
  numObsPerID <- c(as.matrix(table(ids))) # number of observations per ID
  n <- length(numObsPerID) # number of IDs
  orderedID <- rep(1:n, times=numObsPerID) # Makes IDs 1, 2, ..., n
  idEndIdx <- cumsum(numObsPerID) - 1 # -1 for 0 indexing
  if (is.null(covWithVC)) {
    covWithVC <- 1:ncol(X) # default to all columns
  } else {
    c(1, covWithVC + 1) # + 1 for intercept 
  }
  basisList <- getBasisX(X, varyingCov, covWithVC, df=df, degree=degree,
                         basisFunc=basisFunc) 
  Xvar <- basisList$Xvar # X with basis for varying coefficients
  XvartoXColMapping <- basisList$XvartoXColMapping # Column mapping (Xvar -> X)
  interiorKnots <- basisList$interiorKnots # Interior knots of basis functions
  boundaryKnots <- basisList$boundaryKnots # Boundary knots of basis functions

  # updating rCols to correspond to the Xvar
  rCols <- c(0)

  if (is.null(toReturn)) 
    toReturn <- character(0) # default to all columns
 
  start <- Sys.time() 
  output <- FunCDMSampler(iter, counts, Xvar, idEndIdx, rCols,
                          BURN_IN = burnIn,
                          NUM_THIN = thin,
                          ADJ_FREQ = adjustFreq,
                          PROPOSAL_CAP = proposalCap,
                          ADJ_PROPOSALS = adjustProp,
                          RETURN_BURN_IN = returnBurnIn,
                          CAP_PROPOSALS = capProposals,
                          PRINT_PROGRESS = printProgress,
                          TO_RETRUN = toReturn,
                          betaInitial = betaInitial,
                          rInitial = rInitial,
                          priors = priors,
                          proposalVars = proposalVars)
  end <- Sys.time()
  totalTime <- round(difftime(end, start, units = "secs"), 3)
  print(paste("Total time taken:", totalTime, "seconds"))

  output[["interiorKnots"]] <- interiorKnots
  output[["boundaryKnots"]] <- boundaryKnots
  output[["varyingCov"]] <- varyingCov
  output[["colMapping"]] <- colMapping
  output[["XvartoXColMapping"]] <- XvartoXColMapping
  output[["basisFunc"]] <- basisFunc
  output[["df"]] <- df
  output[["catNames"]] <- colnames(counts)
  output[["centerScaleList"]] <- centerScaleList
  
  if (saveNullInBetaCI) { 
    getNullInBetaCIProportion(output, XvartoXColMapping, colMapping, varyingCov,
                              colnames(counts), interiorKnots, 
                              boundaryKnots, basisFunc=basisFunc, df=df, 
                              fileName = nullInBetaCIFileName)
  }
    
  if (saveToFile) {
    saveRDS(output, file = fileName)
    return(NULL)
  }

  return(output)
}

VarCoDM <- function(counts, covariates, ids, varyingCov, iter = 10000,
                    burnIn = 5000, thin = 1, adjustFreq = 250, proposalCap = 0, 
                    returnBurnIn = FALSE, printProgress = TRUE, toReturn = NULL,
                    betaInitial = NULL, priors = NULL, proposalVars = NULL, 
                    covWithVC=NULL, df=4, degree=3, basisFunc=splines::bs,
                    saveToFile = TRUE, fileName = "output.rds", 
                    saveNullInBetaCI = TRUE, 
                    nullInBetaCIFileName = "nullInBetaCI.csv") {
  # error checks
  # No error check made for this, as this function is only used for simulation 
  # studies

  adjustProp <- adjustFreq != 0 # If adjustFreq is 0, no adjustment of proposals
  capProposals <- proposalCap > 0 # If proposalCap is 0, no capping of proposals

  # processing the data
  centerScaleList <- list() # List to store the center and scale of each cov
  centerScaleList[[1]] <- "no c/s" # Intercept does not have c/s
  listIdx <- 2 # Start from 2 to skip intercept
  for (i in 1:ncol(covariates)) {
    if (is.numeric(covariates[, i])) {
      centerScaleList[[listIdx]] <- list(center = mean(covariates[, i]),
                                       scale = sd(covariates[, i]))
      covariates[, i] <- (covariates[, i] - centerScaleList[[listIdx]]$center) / 
                          centerScaleList[[listIdx]]$scale
    } else {
      centerScaleList[[listIdx]] <- "no c/s" # No center/scale for factors
      if (length(levels(covariates[, i])) > 2) {
        for (j in 3:length(levels(covariates[, i]))) {
          # If the factor has more than 2 levels, then we need to account for 
          # the extra dummy variables for the factor levels
          listIdx <- listIdx + 1
          centerScaleList[[listIdx]] <- "no c/s"
        }
      } 
    }
    listIdx <- listIdx + 1
  }
  X <- model.matrix(~., data=covariates) 
  
  colMapping <- list() # List mapping the covariate numbers to colnames
  i <- 1 
  for (name in colnames(X)) {
    colMapping[[i]] <- name
    i <- i + 1
  }
  numObsPerID <- c(as.matrix(table(ids))) # number of observations per ID
  n <- length(numObsPerID) # number of IDs
  orderedID <- rep(1:n, times=numObsPerID) # Makes IDs 1, 2, ..., n
  idEndIdx <- cumsum(numObsPerID) - 1 # -1 for 0 indexing
  if (is.null(covWithVC)) {
    covWithVC <- 1:ncol(X) # default to all columns
  } else {
    c(1, covWithVC + 1) # + 1 for intercept 
  }
  basisList <- getBasisX(X, varyingCov, covWithVC, df=df, degree=degree,
                         basisFunc=basisFunc) 
  Xvar <- basisList$Xvar # X with basis for varying coefficients
  XvartoXColMapping <- basisList$XvartoXColMapping # Column mapping (Xvar -> X)
  interiorKnots <- basisList$interiorKnots # Interior knots of basis functions
  boundaryKnots <- basisList$boundaryKnots # Boundary knots of basis functions

    if (is.null(toReturn)) 
    toReturn <- character(0) # default to all columns
 
  start <- Sys.time() 
  output <- VarCoDMSampler(iter, counts, Xvar, idEndIdx, 
                           BURN_IN = burnIn,
                           NUM_THIN = thin,
                           ADJ_FREQ = adjustFreq,
                           PROPOSAL_CAP = proposalCap,
                           ADJ_PROPOSALS = adjustProp,
                           RETURN_BURN_IN = returnBurnIn,
                           CAP_PROPOSALS = capProposals,
                           PRINT_PROGRESS = printProgress,
                           TO_RETRUN = toReturn,
                           betaInitial = betaInitial,
                           priors = priors,
                           proposalVars = proposalVars)
  end <- Sys.time()
  totalTime <- round(difftime(end, start, units = "secs"), 3)
  print(paste("Total time taken:", totalTime, "seconds"))

  output[["interiorKnots"]] <- interiorKnots
  output[["boundaryKnots"]] <- boundaryKnots
  output[["varyingCov"]] <- varyingCov
  output[["colMapping"]] <- colMapping
  output[["XvartoXColMapping"]] <- XvartoXColMapping
  output[["basisFunc"]] <- basisFunc
  output[["df"]] <- df
  output[["catNames"]] <- colnames(counts)
  output[["centerScaleList"]] <- centerScaleList
  
  if (saveNullInBetaCI) { 
    getNullInBetaCIProportion(output, XvartoXColMapping, colMapping, varyingCov,
                              colnames(counts), interiorKnots, 
                              boundaryKnots, basisFunc=basisFunc, df=df, 
                              fileName = nullInBetaCIFileName)
  }
    
  if (saveToFile) {
    saveRDS(output, file = fileName)
    return(NULL)
  }

  return(output)
}

VarCoZIDM <- function(counts, covariates, ids, varyingCov, iter = 10000, 
                      burnIn = 5000, thin = 1, adjustFreq = 250, 
                      proposalCap = 0, ZIGrouped = TRUE, returnBurnIn = FALSE,
                      printProgress = TRUE, toReturn = NULL, betaInitial = NULL,
                      rInitial = NULL, priors = NULL, proposalVars = NULL, 
                      covWithVC=NULL, df=4, degree=3, basisFunc=splines::bs,
                      saveToFile = TRUE, fileName = "output.rds", 
                      saveNullInBetaCI = TRUE,
                      nullInBetaCIFileName = "nullInBetaCI.csv") {
  # error checks
  # No error check made for this, as this function is only used for simulation 
  # studies

  adjustProp <- adjustFreq != 0 # If adjustFreq is 0, no adjustment of proposals
  capProposals <- proposalCap > 0 # If proposalCap is 0, no capping of proposals

  # processing the data
  centerScaleList <- list() # List to store the center and scale of each cov
  centerScaleList[[1]] <- "no c/s" # Intercept does not have c/s
  listIdx <- 2 # Start from 2 to skip intercept
  for (i in 1:ncol(covariates)) {
    if (is.numeric(covariates[, i])) {
      centerScaleList[[listIdx]] <- list(center = mean(covariates[, i]),
                                       scale = sd(covariates[, i]))
      covariates[, i] <- (covariates[, i] - centerScaleList[[listIdx]]$center) / 
                          centerScaleList[[listIdx]]$scale
    } else {
      centerScaleList[[listIdx]] <- "no c/s" # No center/scale for factors
      if (length(levels(covariates[, i])) > 2) {
        for (j in 3:length(levels(covariates[, i]))) {
          # If the factor has more than 2 levels, then we need to account for 
          # the extra dummy variables for the factor levels
          listIdx <- listIdx + 1
          centerScaleList[[listIdx]] <- "no c/s"
        }
      } 
    }
    listIdx <- listIdx + 1
  }
  X <- model.matrix(~., data=covariates) 
  
  colMapping <- list() # List mapping the covariate numbers to colnames
  i <- 1 
  for (name in colnames(X)) {
    colMapping[[i]] <- name
    i <- i + 1
  }
  numObsPerID <- c(as.matrix(table(ids))) # number of observations per ID
  n <- length(numObsPerID) # number of IDs
  orderedID <- rep(1:n, times=numObsPerID) # Makes IDs 1, 2, ..., n
  idEndIdx <- cumsum(numObsPerID) - 1 # -1 for 0 indexing
  if (is.null(covWithVC)) {
    covWithVC <- 1:ncol(X) # default to all columns
  } else {
    c(1, covWithVC + 1) # + 1 for intercept 
  }
  basisList <- getBasisX(X, varyingCov, covWithVC, df=df, degree=degree,
                         basisFunc=basisFunc) 
  Xvar <- basisList$Xvar # X with basis for varying coefficients
  XvartoXColMapping <- basisList$XvartoXColMapping # Column mapping (Xvar -> X)
  interiorKnots <- basisList$interiorKnots # Interior knots of basis functions
  boundaryKnots <- basisList$boundaryKnots # Boundary knots of basis functions

  if (is.null(toReturn)) 
    toReturn <- character(0) # default to all columns
 
  start <- Sys.time() 
  output <- VarCoZIDMSampler(iter, counts, Xvar, idEndIdx,
                             BURN_IN = burnIn,
                             NUM_THIN = thin,
                             ADJ_FREQ = adjustFreq,
                             PROPOSAL_CAP = proposalCap,
                             ZI_GROUPED = ZIGrouped,
                             ADJ_PROPOSALS = adjustProp,
                             RETURN_BURN_IN = returnBurnIn,
                             CAP_PROPOSALS = capProposals,
                             PRINT_PROGRESS = printProgress,
                             TO_RETRUN = toReturn,
                             betaInitial = betaInitial,
                             priors = priors,
                             proposalVars = proposalVars)
  end <- Sys.time()
  totalTime <- round(difftime(end, start, units = "secs"), 3)
  print(paste("Total time taken:", totalTime, "seconds"))

  output[["interiorKnots"]] <- interiorKnots
  output[["boundaryKnots"]] <- boundaryKnots
  output[["varyingCov"]] <- varyingCov
  output[["colMapping"]] <- colMapping
  output[["XvartoXColMapping"]] <- XvartoXColMapping
  output[["basisFunc"]] <- basisFunc
  output[["df"]] <- df
  output[["catNames"]] <- colnames(counts)
  output[["centerScaleList"]] <- centerScaleList
  
  if (saveNullInBetaCI) { 
    getNullInBetaCIProportion(output, XvartoXColMapping, colMapping, varyingCov,
                              colnames(counts), interiorKnots, 
                              boundaryKnots, basisFunc=basisFunc, df=df, 
                              fileName = nullInBetaCIFileName)
  }
    
  if (saveToFile) {
    saveRDS(output, file = fileName)
    return(NULL)
  }

  return(output)
}