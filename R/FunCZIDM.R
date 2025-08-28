# Author: Brody Erlandson
#
# This file contains the wrapper function for the FunC-ZIDM regression model
# and the functions to extract the results from the model.


#' @title FunCZIDM
#' @description This function is used to sample from the FunC-ZIDM regression
#'              model. 
#' @param counts (integer matrix) Each categories counts for each sample.
#' @param covariates (data.frame) The covariates used for the regression, 
#' should not include the intercept. These should correspond to the columns of
#' \code{counts}.
#' @param ids (integer vector) The index of subject ids in the 'covariates'. 
#' @param varyingCov (numeric vector) The index of the covariate the varying 
#' coefficients will be a function of in \code{covariates}. This is usually 
#' time.
#' @param iter (integer, default=10000) Total number of iterations.
#' @param burnIn (integer, default=5000) Number of burn-in iterations, must be
#'  less than \code{iter}.
#' @param thin (integer, default=1) Number of iterations to thin by.
#' @param adjustFreq (numeric, default=250) Frequency at which 
#' \eqn{\beta_{jpd}^*} proposals are adjusted.
#' @param proposalCap (numeric, default=0) Cap for proposal variances, if 0 
#' there will be no cap.
#' @param ZIGrouped (logical, default=TRUE) Whether to use the zero inflation 
#' indicator for all samples of an id (TRUE), or each sample have their own 
#' indicator (FALSE).
#' @param returnBurnIn (logical, default=FALSE) Whether to return burn-in 
#' samples.
#' @param printProgress (logical, default=TRUE) Whether to print progress bar.
#' @param toReturn (vector, default=NULL) Vector of parameter names to 
#' return in addition to the default parameters. The default parameters are:
#' \itemize{
#' \item{"beta"}: The \eqn{\beta_{jpd}^*} coefficients.
#' \item{"betaAcceptProp"}: The acceptance proportion for the 
#' \eqn{\beta_{jpd}^*} proposals.
#' \item{"rAcceptProp"}: The acceptance proportion for the \eqn{r_{jp}} 
#' proposals.
#' \item{"rMeans"}: The mean of the \eqn{r_{jp}} coefficients.
#' \item{"etaMeanPropZeros"}: The average level of zero-inflation indicated.
#' }
#' Returnable parameters are:
#' \itemize{
#' \item{"eta"}: The zero-inflation indicator.
#' \item{"c"}: The zero-inflation concentration parameter.
#' \item{"u"}: The auxiliary variable.
#' \item{"r"}: The individual- and category-specific coefficients.
#' \item{"phi"}: The individual- and category-specific coefficients sd
#' \item{"lambda"}: The local parameter of the horseshoe prior.
#' \item{"nu"}: The local parameter auxiliary variable of the horseshoe prior.
#' \item{"tau"}: The global parameter of the horseshoe prior.
#' \item{"xi"}: The global parameter auxiliary variable of the horseshoe prior.
#' \item{"RA"}: The mean estimated relative abundance of each sample.
#' }
#' If NULL, only the default parameters are returned.
#' @param betaInitial (matrix, default=NULL) Initialization values for
#' \eqn{beta_{jpd}^*}. If NULL, the initial values will be drawn from 
#' Unif(-.75,.75).
#' @param rInitial (matrix, default=NULL) Initialization values for
#' \eqn{r_{jp}}. If NULL, the initial values will be drawn from Unif(-.05, .05).
#' @param priors (list, default=NULL) List with prior hyperparameters 
#' specifications. The hyperparameters that can be specified are:
#' \itemize{
#' \item{"first beta sd"}: The prior sd for the first 
#' \eqn{\beta_{jpd}^*} coefficient (\eqn{\beta_{j00}^*}).
#' \item{"a"}: The prior shape parameter for the variance of the \eqn{r_{jp}}.
#' \item{"b"}: The prior rate parameter for the variance of the \eqn{r_{jp}}.
#' \item{"alpha"}: The prior fisrt shape parameter for the probability of the 
#' zero-inflation indicator.
#' \item{"beta"}: The prior second shape parameter for the probability of the
#' zero-inflation indicator.
#' \item{"kappaShape"}: The prior shape parameter for the \eqn{\kappa_{jp}}, 
#' i.e. the student-t portion of the regularized horseshoe.
#' \item{"kappaRate"}: The prior rate parameter for the \eqn{\kappa_{jp}}, 
#' i.e. the student-t portion of the regularized horseshoe.
#' }
#' NULL gives the default values, which are:
#' \itemize{
#' \item{"first beta sd"}: 1
#' \item{"a"}: 3
#' \item{"b"}: 9
#' \item{"alpha"}: .01
#' \item{"beta"}: 10
#' }
#' @param proposalVars (list, default=NULL) List with proposal 
#' sd for \eqn{\beta_{jpd}^*} and \eqn{r_{jp}}. If NULL, the default proposal 
#' sd are used, which are:
#' \itemize{
#' \item{"beta proposal sd"}: .3 
#' \item{"r proposal sd"}: 1
#' }
#' @param covWithoutVC (vector, default=NULL) Indices of the \code{covariates}
#' that will not have varying coefficients. If NULL, all covariates will have
#' varying coefficients. Intercept will always have varying coefficients.
#' @param df (numeric, default=4) Degrees of freedom for the basis function.
#' @param degree (numeric, default=3) Degree for the basis function.
#' @param basisFunc (function, default=splines::bs) Function to generate basis 
#' functions.
#' @param saveToFile (logical, default=TRUE) Whether to save output to file.
#' @param fileName (character, default="output.rds") File name for 
#' output.
#' @param saveNullInBetaCI (logical, default=TRUE) Whether to a csv that 
#' contains the proportion of the \eqn{\beta_{jp}(t)} coefficients that contain
#' the null value in their credible interval. Gives an idea of which covariate
#' and category combinations may have meaningful effects.
#' @param nullInBetaCIFileName (character, default="nullInBetaCI.csv") 
#' File name for null CI statistics.
#' @return Depending on saveToFile, either returns NULL after saving or an 
#' output list.
#' @export 
FunCZIDM <- function(counts, covariates, ids, varyingCov, 
                     iter = 10000, burnIn = 5000, thin = 1, adjustFreq = 250, 
                     proposalCap = 0, ZIGrouped = TRUE, returnBurnIn = FALSE,
                     printProgress = TRUE, toReturn = NULL, betaInitial = NULL, 
                     rInitial = NULL, priors = NULL, proposalVars = NULL, 
                     covWithoutVC=NULL, df=4, degree=3, basisFunc=splines::bs,
                     saveToFile = TRUE, fileName = "output.rds", 
                     saveNullInBetaCI = TRUE, 
                     nullInBetaCIFileName = "nullInBetaCI.csv") {
  # error checks
  args <- as.list(environment())
  do.call(errorChecksFunCZIDM, args) 

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

  covWithVC <- 1:ncol(X)
  # dropping those that in covWithoutVC
  for (cov in covWithoutVC) {
    toDrop <- which(grepl(colnames(infantCovariates)[cov], colnames(X)))
    covWithVC <- covWithVC[-toDrop]
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
  output <- FunCZIDMSampler(iter, counts, Xvar, idEndIdx, rCols,
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
    saveRDS(output, file = fileName, compress="xz")
    return(NULL)
  }

  return(output)
}

#' @title Get \eqn{\beta_{jp}(t)} Functions
#' @description Extracts the \eqn{\beta_{jp}(t)} function samples from the 
#' output of \code{FunCZIDM}
#' @param output (list) Output list returned from \code{FunCZIDM}
#' @param cov (integer, default=NULL) The index of \code{covariates} from 
#' \code{FunCZIDM} for which \eqn{\beta_{jp}(t)} functions to extract. 
#' Intercept if NULL.
#' @return A 3-dimensional array of \eqn{\beta_{jp}(t)} function samples for 
#' the specified covariate.
#' @export
getBetaFunctions <- function(output, cov = NULL) {
  checkList(output, "output")
  numCov <- length(output$colMapping)
  if (!is.null(cov)) {
    checkNonNegativeInteger(cov, "cov")
    checkLessThanOrEqualTo(cov, numCov, "cov")
  } else { # cov is NULL
    cov <- 1 # Default to intercept
  } 
  
  minVarCov <- min(output$varyingCov)
  maxVarCov <- max(output$varyingCov)
  testPoints <- seq(minVarCov,  maxVarCov, length.out = 250)
  basis <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                           knots=output$interiorKnots,
                           Boundary.knots=output$boundaryKnots)
  basis <- cbind(1, basis) # Add intercept
  numDf <- output$df + 1 # + 1 for intercept

  startTV <- match(cov, output$XvartoXColMapping) # First instance of the cov
  if (cov == output$XvartoXColMapping[startTV + 1]) {
    # If the next instance in the mapping is not the next covariate, then
    # the covariate does not have varying coefficients. So startTV == endTV.
    endTV <- startTV + numDf - 1
    betaVCFits <- getBetaVC(output$beta, c(startTV, endTV), basis)
  } else {
    betaVCFits <- output$beta[startTV, , , drop=FALSE]
    betaVCFits <- array(rep(betaVCFits, each = 100), 
                        dim = c(100, dim(betaVCFits)[2], dim(betaVCFits)[3]))
  }
  colnames(betaVCFits) <- output$catNames # Set column names to categories

  return(betaVCFits)
}

#' @title Calculates the Relateive Abundance
#' @description This function calculates the \eqn{RA_{j}(t)} from the from the 
#' output of \code{FunCZIDM}
#' @param output (list) Output list from \code{FunCZIDM}
#' @param covProfile (numeric vector, default=NULL) Covariate profile (default 
#' is baseline, i.e. all zeros).
#' @return A 3-dimensional array of \eqn{RA_{j}(t)} samples.
#' @export 
calcRA <- function(output, covProfile=NULL) {
  checkList(output, "output")
  numCov <- length(output$colMapping) 
  if (!is.null(covProfile)) {
    checkVector(covProfile, "covProfile", numericOnly=TRUE)
    checkLength(covProfile, numCov, "covProfile")
    getCenterScaledCovProfile(covProfile, output$centerScaleList)
  } else { # covProfile is NULL
    covProfile <- c(1, rep(0, numCov)) # default to baseline profile
  }

  minVarCov <- min(output$varyingCov)
  maxVarCov <- max(output$varyingCov)
  testPoints <- seq(minVarCov, maxVarCov, length.out = 250)
  basis <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                           knots=output$interiorKnots,
                           Boundary.knots=output$boundaryKnots)
  basis <- cbind(1, basis) # Add intercept
  numDf <- output$df + 1 # + 1 for intercept

  return(getRAFits(output$beta, basis, covProfile, 
                   output$XvartoXColMapping))
}

#' @title Calculate the Multipicative Change in Relateive Abundance
#' @description This function calculates the \eqn{\Delta_vRA_{jp}(t)} statistic
#' from the samples of \code{FunCZIDM}. 
#' @param output (list) Output list from \code{FunCZIDM}
#' @param change (numeric) The change of the covariate.
#' @param covProfile (numeric vector, default=NULL) Covariate profile.
#' If NULL, the baseline profile is used (ie all zeros).
#' @param forCovs (integer vector, default=NULL) Indices of \code{covariates}
#' for which \eqn{\Delta_vRA_{jp}(t)} to be calculated.
#' @param forCats (integer vector, default=NULL) Indices of categories  
#' for which \eqn{\Delta_vRA_{jp}(t)} is calculated.
#' @return A list of 3-dimensional array \eqn{\Delta_vRA_{jp}(t)} samples
#' @export 
calcDeltaRA <- function(output, change, covProfile=NULL, forCovs=NULL, 
                        forCats=NULL) {
  checkList(output, "output")
  numCov <- length(output$colMapping) 
  numCat <- length(output$catNames)
  if (!is.null(covProfile)) {
    checkVector(covProfile, "covProfile", numericOnly=TRUE)
    checkLength(covProfile, numCov, "covProfile")
    getCenterScaledCovProfile(covProfile, output$centerScaleList)
  } else { # covProfile is NULL
    covProfile <- c(1, rep(0, numCov)) # default to baseline profile
  }
  if (!is.null(forCovs)) {
    checkValuesBetween(forCovs, 1, numCov)
  } else {
    forCovs <- 2:numCov # Default to all covariates (except intercept)
  }
  if (!is.null(forCats)) {
    checkValuesBetween(forCats, 1, numCat)
  } else {
    forCats <- 1:numCat # Default to all categories
  }
  for (cov in forCovs) {
    getCenterScaledChange(change, output$centerScaleList, cov, covProfile)
  }

  minVarCov <- min(output$varyingCov)
  maxVarCov <- max(output$varyingCov)
  testPoints <- seq(minVarCov,  maxVarCov, length.out = 250)
  basis <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                            knots=output$interiorKnots,
                            Boundary.knots=output$boundaryKnots)
  basis <- cbind(1, basis) # Add intercept
  numDf <- output$df + 1 # + 1 for intercept

  # getting fit and sumExpFit for statistic calculations
  numCov <- length(output$colMapping) # Number of covariates
  fit <- getFit(output$beta, basis, covProfile, output$XvartoXColMapping)
  sumExpFit <- getSumExpFit(fit)

  if (is.null(forCats)) { # Default to all categories
    toReturn <- list()
    startTV <- match(2, output$XvartoXColMapping) # First instance of 1st cov
    numCat <- length(output$catNames) # number of categories
    i <- 1
    for (cov in forCovs) {
      while (cov > output$XvartoXColMapping[startTV]) {
        # If the cov is not current startTV, then the current startTV is not
        # not a covariate that is wanted to calculate deltaRA for.
        startTV <- startTV + 1
      }

      # if the next instance in the mapping is not the next covariate, then
      # the covariate does not have varying coefficients. So startTV == endTV.
      if (cov == output$XvartoXColMapping[startTV + 1]) {
        endTV <- startTV + numDf - 1
        betaVCFits <- getBetaVC(output$beta, c(startTV, endTV), basis)
      } else {
        endTV <- startTV
        betaVCFits <- output$beta[startTV, , , drop=FALSE]
        betaVCFits <- array(rep(betaVCFits, each = 100), 
                        dim = c(100, dim(betaVCFits)[2], dim(betaVCFits)[3]))
      }

      covariteName <- output$colMapping[[cov]]
      toReturn[[covariteName]] <- getAllCatDeltaRA(betaVCFits, fit, sumExpFit,
                                                   change=change[i])
      colnames(toReturn[[covariteName]]) <- output$catNames
      i <- i + 1
      startTV <- endTV + 1 # Move to next covariate
    }
  } else { # Only for specific categories
    toReturn <- list()
    startTV <- match(2, output$XvartoXColMapping) # First instance of 1st cov
    numCat <- length(output$catNames) # number of categories
    i <- 1
    for (cov in forCovs) {
      while (cov > output$XvartoXColMapping[startTV]) {
        # If the cov is not current startTV, then the current startTV is not
        # not a covariate that is wanted to calculate deltaRA for.
        startTV <- startTV + 1
      }

      # if the next instance in the mapping is not the next covariate, then
      # the covariate does not have varying coefficients. So startTV == endTV.
      if (cov == output$XvartoXColMapping[startTV + 1]) {
        endTV <- startTV + numDf - 1
        betaVCFits <- getBetaVC(output$beta, c(startTV, endTV), basis)
      } else {
        endTV <- startTV
        betaVCFits <- output$beta[startTV, , , drop=FALSE]
        betaVCFits <- array(rep(betaVCFits, each = 100), 
                        dim = c(100, dim(betaVCFits)[2], dim(betaVCFits)[3]))
      }
      covariateName <- output$colMapping[[cov]]

      # combine the deltaRA for the specified categories as an array
      dRA <- abind::abind(lapply(forCats, function(cat)
                                     getDeltaRA(cat, betaVCFits, fit, sumExpFit,
                                                change=change[i])),
                  along=3)
      toReturn[[covariateName]] <- aperm(dRA, c(1, 3, 2))
      colnames(toReturn[[covariateName]]) <- output$catNames[forCats]
      i <- i + 1
      startTV <- endTV + 1 # Move to next covariate
    }
  }

  return(toReturn)
}

#' @title Calculate \eqn{\alpha}-diversity Statistic
#' @description Calculates the \eqn{\alpha_p(t)} statistic using the 
#' Hill diversity index  from the \code{FunCZIDM} samples.
#' @param output (list) Output list from \code{FunCZIDM}
#' @param l (numeric, default=0) Parameter l the diversity weight parameter. 
#' If 0, the Shannon diversity index is calculated, if 1, the count diversity 
#' index is calculated. Closer to 1 gives more weight to rare categories. 
#' @param covProfile (numeric vector, default=NULL) Covariate profile. If NULL, 
#' the baseline profile is used (ie all zeros).
#' @return A 3-dimensional array of \eqn{\alpha_p(t)} samples
#' @export 
calcAlphaDiv <- function(output, l=0, covProfile=NULL) {
  checkList(output, "output")
  numCov <- length(output$colMapping) 
  if (!is.null(covProfile)) {
    checkVector(covProfile, "covProfile", numericOnly=TRUE)
    checkLength(covProfile, numCov, "covProfile")
  } else { # covProfile is NULL
    covProfile <- c(1, rep(0, numCov)) # default to baseline profile
  }
  if (l < 0 | l > 1)
    stop(paste("l must be in [0, 1)."))

  minVarCov <- min(output$varyingCov)
  maxVarCov <- max(output$varyingCov)
  testPoints <- seq(minVarCov,  maxVarCov, length.out = 250)
  basis <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                            knots=output$interiorKnots,
                            Boundary.knots=output$boundaryKnots)
  basis <- cbind(1, basis) # Add intercept
  numDf <- output$df + 1 # + 1 for intercept

  return(getHillDiversity(output$beta, basis, covProfile, 
                          output$XvartoXColMapping, l=l))
}

#' @title Calculate The Multiplicative Change in \eqn{\alpha}-diversity 
#' Statistic
#' @description Calculates the \eqn{\Delta_v\alpha_p(t)} statistic using the 
#' Hill diversity index from the \code{FunCZIDM} samples.
#' @param output (list) Output list from \code{FunCZIDM}
#' @param change (numeric) The amount change of the covariate.
#' @param l (numeric, default=0) Parameter l the diversity weight parameter. 
#' If 0, the Shannon diversity index is calculated, if 1, the count diversity 
#' index is calculated. Closer to 1 gives more weight to rare categories. 
#' @param covProfile (numeric vector, default=NULL) Covariate profile. If NULL, 
#' the baseline profile is used (ie all zeros).
#' @param forCovs (numeric vector, default=NULL) Indices of \code{covariates}
#' for which \eqn{\Delta_v\alpha_p(t)} to be calculated.
#' @return A list of 3-dimensional arrays of \eqn{\Delta_v\alpha_p(t)} samples
#' @export 
calcDeltaAlphaDiv <- function(output, change, l = 0, covProfile=NULL, 
                              forCovs=NULL) {
  checkList(output, "output")
  numCov <- length(output$colMapping)
  if (!is.null(covProfile)) {
    checkVector(covProfile, "covProfile", numericOnly=TRUE)
    checkLength(covProfile, numCov, "covProfile")
  } else { # covProfile is NULL
    covProfile <- c(1, rep(0, numCov)) # default to baseline profile
  }
  if (!is.null(forCovs)) {
    checkValuesBetween(forCovs, 1, numCov)
  } else {
    forCovs <- 2:numCov # Default to all covariates (except intercept)
  }
  if (l < 0 | l > 1)
    stop(paste("l must be in [0, 1)."))

  minVarCov <- min(output$varyingCov)
  maxVarCov <- max(output$varyingCov)
  testPoints <- seq(minVarCov,  maxVarCov, length.out = 250)
  basis <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                            knots=output$interiorKnots,
                            Boundary.knots=output$boundaryKnots)
  basis <- cbind(1, basis) # Add intercept
  numDf <- output$df + 1 # + 1 for intercept
  
  fit <- getFit(output$beta, basis, covProfile, output$XvartoXColMapping)
  sumExpFit <- getSumExpFit(fit)
  
  toReturn <- list()
  startTV <- match(2, output$XvartoXColMapping) # First instance of the first covariate
  numCat <- length(output$catNames) # number of categories
  i <- 1
  for (cov in forCovs) {
    while (cov > output$XvartoXColMapping[startTV]) {
      # If the cov is not current startTV, then the current startTV is not
      # not a covariate that is wanted to calculate deltaRA for.
      startTV <- startTV + 1
    }

    # if the next instance in the mapping is not the next covariate, then
    # the covariate does not have varying coefficients. So startTV == endTV.
    if (cov == output$XvartoXColMapping[startTV + 1]) {
      endTV <- startTV + numDf - 1
      betaVCFits <- getBetaVC(output$beta, c(startTV, endTV), basis)
    } else {
      endTV <- startTV
      betaVCFits <- output$beta[startTV, , , drop=FALSE]
      betaVCFits <- array(rep(betaVCFits, each = 100), 
                        dim = c(100, dim(betaVCFits)[2], dim(betaVCFits)[3]))
    }

    toReturn[[output$colMapping[[cov]]]] <- getDeltaHillDiversity(betaVCFits, 
                                                                  fit, 
                                                                  sumExpFit,
                                                                  l=l,
                                                               change=change[i])
    i <- i + 1
    startTV <- endTV + 1 # Move to next covariate
  }

  return(toReturn)
}

#' @title Combine \code{FunCZIDM} Outputs
#' @description Combines multiple FunCZIDM output files into a single output 
#' list.
#' @param outputFiles (Character vector) Character vector of file paths to the 
#' saved outputs.
#' @param saveToFile (logical, default=TRUE) Whether to save output to file.
#' @param fileName (character, default="combinedOutputs.rds") File name
#'  for output.
#' @param compress (character, default="xz") If \code{saveToFile} is TRUE, the 
#' compression method for saving the output file. Options are "xz", "gzip", or 
#' "bzip2".
#' @return If \code{saveToFile} is TRUE, returns NULL after saving the combined
#' output to file. If \code{saveToFile} is FALSE, returns the combined output 
#' list.
#' @export 
combineOutputs <- function(outputFiles, saveToFile = TRUE, 
                           fileName = "combinedOutputs.rds",
                           compress = "xz") {
  combinedOutput <- list()
  numOutputs <- length(outputFiles)
  for (f in outputFiles) {
    if (!file.exists(f)){
      stop(paste(f, "does not exist"))
    }
  }
  output <- readRDS(outputFiles[1]) # Read first output to initialize 
  numSamples <- dim(output$beta)[3] # Number of samples in the output

  # these keys are the same for all outputs
  keys <- c("interiorKnots", "boundaryKnots", "varyingCov", "colMapping",
            "XvartoXColMapping", "basisFunc", "df", "catNames",
            "centerScaleList")
  for (key in keys) {
    combinedOutput[[key]] <- output[[key]]
  }

  # getting the keys who values are different for each output
  outputKeys <- setdiff(names(output), keys)

  # Initialize combinedOutput with the first output
  for (key in outputKeys) {
    if (grepl("prop|mean|prop", key, ignore.case=TRUE)) {
      # if prop or mean are in the key, then we must combine by averaging
      combinedOutput[[key]] <- output[[key]]/numOutputs
    } else {
      # otherwise, we have samples of a parameter, so we can just combine
      dim <- dim(output[[key]])
      dim[length(dim)] <- numSamples*numOutputs # Combine samples across outputs
      combinedOutput[[key]] <- array(0, dim = dim) # Initialize with zeros

      # interleave the samples from each output
      sOrig = 1
      for (s in seq(1, numSamples*numOutputs, by = numOutputs)) {
        if (length(dim(output[[key]])) == 3) {
          combinedOutput[[key]][,,s] <- output[[key]][,,sOrig]
        } else {
          combinedOutput[[key]][,s] <- output[[key]][,sOrig]
        }
        sOrig <- sOrig + 1
      }
    }
  }

  # Combine the outputs
  fileNum <- 2 # Start from the second file
  for (file in outputFiles[-1]) { # Skip the first file as it is already read
    output <- readRDS(file)
    for (key in outputKeys) {
      if (grepl("prop|mean|prop", key, ignore.case=TRUE)) {
        # if prop or mean are in the key, then we must combine by averaging
        combinedOutput[[key]] <- combinedOutput[[key]] 
                                 + (output[[key]]/numOutputs)
      } else {
        # otherwise, we have samples of a parameter, so we can just combine
        
        # interleave the samples from each output
        sOrig = 1
        for (s in seq(fileNum, numSamples*numOutputs, by = numOutputs)) {
          if (length(dim(output[[key]])) == 3) {
            combinedOutput[[key]][,,s] <- output[[key]][,,sOrig]
          } else {
            combinedOutput[[key]][,s] <- output[[key]][,sOrig]
          }
          sOrig <- sOrig + 1
        }
      }
    }
    fileNum <- fileNum + 1 # Increment the file number
  }

  if (saveToFile) {
    saveRDS(combinedOutput, file = fileName, compress = compress)
    return(NULL) # Return NULL if saved to file
  }
  return(combinedOutput)   
}

#' @title Generate Data
#' @description Generates synthetic data for testing the model.
#' @param n (integer) Number of individuals.
#' @param c (integer) Number of categories.
#' @param p (integer) Number of functional covariates.
#' @param totalCountRange (numeric vector, default=c(100, 200)) Range of total 
#' counts.
#' @param numActive (integer, default=4) Number of categories with active 
#' covariates. Alway 1:numActive categories are active.
#' @return A list containing the generated data: counts, covariates, ids, 
#' timepoints and the true parameters.
#' @export 
generateData <- function(n, c, p, totalCountRange = c(100, 200),
                         numActive = 4) {
  ### Generate predictor matrices ###
  timePoints <- c()
  numObsPerID <- numeric(n)
  cumNumObsPerID <- numeric(n)
  sumNumObsPerID <- 0
  for (i in 1:n) {
    numObsPerID[i] <- sample(3:10, 1)
    cumNumObsPerID[i] <- sumNumObsPerID
    sumNumObsPerID <- sumNumObsPerID + numObsPerID[i]
    timePoints <- c(timePoints, sort(runif(numObsPerID[i], 0, 10)))
  }
  nTotal <- sum(numObsPerID)
  X <- matrix(rnorm(nTotal*p), nTotal, p)
  X <- scale(X)
  X <- cbind(1, X) # Add intercept
  colnames(X) <- c("intercept", paste0("cov", 1:p))
  personID <- rep(1:n, numObsPerID)
  Y <- makeBlockMatrix(X, c(1), personID)

  ### Set true parameter values ###
  f1 <- function(t) { (-.2*(t-5)^2 + 5)/7 }
  f2 <- function(t) { 1/(1.75 + exp(-1.25*(t - 5))) }
  f3 <- function(t) { .07*t }
  f4 <- function(t) { .5 + 0*t }
  f5 <- function(t) { 0*t }
  funcList <- list(f1, f2, f3, f4, f5)

  # Generate time-varying coefficients and true response w/ intercept or r
  tvMat <- matrix(0, nTotal, c)
  funcMat <- matrix(runif(p*c, .5, 5.5), p, c)
  funcMat[,(numActive+1):c] <- 5
  for (cov in 2:(p+1)) {
    for (i in 1:c) {
      f <- funcList[[round(funcMat[cov-1,i])]]
      plusMinus <- funcMat[cov-1,i] - round(funcMat[cov-1,i])
      tvMat[,i] <- tvMat[,i] + 
                   ifelse(plusMinus >= 0, 1, -1)*f(timePoints)*X[,cov]
    }
  }
  
  # Generate functional intercepts and r
  betaIntercepts <- matrix(NA, length(timePoints), c)
  betaTrue <- matrix(runif(c, .5, 4.5), 1, c)
  for (i in 1:c) {
    f <- funcList[[round(betaTrue[1,i])]]
    plusMinus <- betaTrue[1,i] - round(betaTrue[1,i])
    betaIntercepts[,i] <- ifelse(plusMinus >= 0, 1, -1)*f(timePoints)
  }
  rTrue <- matrix(runif(n*c, -.05, .05), n, c)

  # Generate gamma values
  gamma <- exp(betaIntercepts + Y%*%rTrue + tvMat)

  # Generate zero inflation
  numLowZeros <- floor(.30*c)
  numMidZeros <- floor(.50*c)
  numHighZeros <- c - numLowZeros - numMidZeros
  zeroProps <- c(runif(numLowZeros, 0, .15), runif(numMidZeros, .15, .75), 
                 runif(numHighZeros, .75, .90))
  zeroProps <- sample(zeroProps, c, replace=FALSE)
  oneZeroMultiplier <- matrix(1, nTotal, c)
  for (cat in 1:c){
    zeroInf <- rbinom(n, 1, 1-zeroProps[cat])
    personWithOne <- which(zeroInf == 1)
    while (length(personWithOne) < 5) {
      zeroInf <- rbinom(n, 1, 1-zeroProps[cat])
      personWithOne <- which(zeroInf == 1)
    }
    personWithZero <- which(zeroInf == 0)
    oneZeroMultiplier[,cat] <- ifelse(personID %in% personWithZero, 0, 1)
  }

  # check if any person has all zeros
  allZeroIndicators = rowSums(oneZeroMultiplier) == 0
  personWithAllZero <- unique(personID[which(allZeroIndicators)])
  for (person in personWithAllZero) {
    personIdx = which(personID == person)
    while (sum(rowSums(oneZeroMultiplier[personIdx,])) == 0) {
      for (cat in 1:c){
        zeroInf <- rbinom(1, 1, 1-zeroProps[cat])
        oneZeroMultiplier[personIdx, cat] <- rep(ifelse(zeroInf == 0, 0, 1),
                                                 length(personIdx))
      }
    }
  }

  # apply the zero inflation multiplier
  gamma <- gamma*oneZeroMultiplier

  # Generate RA values from a Dirichlet distribution
  RA <- t(apply(gamma, 1, function(a) MCMCpack::rdirichlet(1, a)))

  # Generate N (total counts) from a Uni
  N <- round(runif(nTotal, totalCountRange[1], totalCountRange[2]))

  # Combine N and RA into a data frame for easier iteration
  parms <- data.frame(N=N)
  parms$RA <- split(RA, row(RA))

  # Generate counts from a Multinomial distribution
  getCounts <- function(N, RA) {
    rmultinom(1, N, RA)
  }

  # Apply getCounts to each row of parms
  counts <- t(mapply(getCounts, parms$N, parms$RA, SIMPLIFY=TRUE))
  counts <- as.matrix(counts)

  covariates <- data.frame(X[, -1]) # Remove intercept column
  colnames(covariates) <- colnames(X)[-1] # Set column names to covariates
  return(list(counts = counts, ids = personID, covariates = covariates,
              timePoints = timePoints, betaIntercepts = betaTrue, rTrue = rTrue, 
              betaFuncMat = funcMat, trueZI = sum(gamma == 0)/length(gamma),
              trueRA=RA, funcList=funcList)) 
}