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
#' \item{"kappaShape"}: 100
#' \item{"kappaRate"}: 900
#' }
#' @param proposalVars (list, default=NULL) List with proposal 
#' sd for \eqn{\beta_{jpd}^*} and \eqn{r_{jp}}. If NULL, the default proposal 
#' sd are used, which are:
#' \itemize{
#' \item{"beta proposal sd"}: .1
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
#' output list. The output list contains:
#' \itemize{
#' \item{"beta"}: The \eqn{\beta_{jpd}^*} samples.
#' \item{"betaAcceptProp"}: The acceptance proportion for the 
#' \eqn{\beta_{jpd}^*} proposals after burn-in.
#' \item{"rAcceptProp"}: The acceptance proportion for the \eqn{r_{jp}} 
#' proposals after burn-in.
#' \item{"rMeans"}: The mean of the \eqn{r_{jp}} samples.
#' \item{"etaMeanPropZeros"}: The average level of zero-inflation indicated.
#' \item{"interiorKnots"}: The interior knots of the basis functions.
#' \item{"boundaryKnots"}: The boundary knots of the basis functions.
#' \item{"df"}: The degrees of freedom used for the basis functions.
#' \item{"basisFunc"}: The basis function used.
#' \item{"varyingCov"}: The values of the covariate that the varying 
#' coefficients are a function of.
#' \item{"colMapping"}: A list mapping the covariate indices to the column names
#' of the design matrix.
#' \item{"XvartoXColMapping"}: A vector mapping the columns of the design matrix
#' with the basis design matrix.
#' \item{"catNames"}: The category names from the counts.
#' \item{"centerScaleList"}: A list containing the center and scale of each 
#' continuous covariate. If a covariate is a factor, then "no c/s" is stored.
#' }
#' @export 
FunCZIDM <- function(counts, covariates, ids, varyingCov, 
                     iter = 10000, burnIn = 5000, thin = 1, adjustFreq = 250,
                     ZIGrouped = TRUE, returnBurnIn = FALSE, 
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
                            PROPOSAL_CAP = 0, # no capping the proposals
                            ZI_GROUPED = ZIGrouped,
                            ADJ_PROPOSALS = adjustProp,
                            RETURN_BURN_IN = returnBurnIn,
                            CAP_PROPOSALS = FALSE, # no capping the proposals
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
  rownames(betaVCFits) <- testPoints 

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
    covProfile <- getCenterScaledCovProfile(covProfile, output$centerScaleList)
  } else { # covProfile is NULL
    covProfile <- c(1, rep(0, numCov - 1)) # default to baseline profile
  }

  minVarCov <- min(output$varyingCov)
  maxVarCov <- max(output$varyingCov)
  testPoints <- seq(minVarCov, maxVarCov, length.out = 250)
  basis <- output$basisFunc(testPoints, df=output$df, intercept=FALSE,
                           knots=output$interiorKnots,
                           Boundary.knots=output$boundaryKnots)
  basis <- cbind(1, basis) # Add intercept
  numDf <- output$df + 1 # + 1 for intercept

  ret <- getRAFits(output$beta, basis, covProfile, 
                   output$XvartoXColMapping)
  dimnames(ret) <- list(testPoints, output$catNames, NULL)
  return(ret)
}

#' @title Calculate the Multipicative Change in Relateive Abundance
#' @description This function calculates the \eqn{\Delta_vRA_{jp}(t)} statistic
#' from the samples of \code{FunCZIDM}. 
#' @param output (list) Output list from \code{FunCZIDM}
#' @param change (numeric vector) The change of the covariate.
#' @param forCovs (integer vector) Indices of \code{covariates}
#' for which \eqn{\Delta_vRA_{jp}(t)} to be calculated.
#' @param covProfile (numeric vector, default=NULL) Covariate profile.
#' If NULL, the baseline profile is used (ie all zeros).
#' @param forCats (integer vector, default=NULL) Indices of categories  
#' for which \eqn{\Delta_vRA_{jp}(t)} is calculated. If NULL, then all 
#' categories are used.
#' @return A list of 3-dimensional array \eqn{\Delta_vRA_{jp}(t)} samples
#' @export 
calcDeltaRA <- function(output, change, forCovs, covProfile=NULL, 
                        forCats=NULL) {
  checkList(output, "output")
  numCov <- length(output$colMapping) 
  numCat <- length(output$catNames)
  if (!is.null(covProfile)) {
    checkVector(covProfile, "covProfile", numericOnly=TRUE)
    checkLength(covProfile, numCov, "covProfile")
    covProfile <- getCenterScaledCovProfile(covProfile, output$centerScaleList)
  } else { # covProfile is NULL
    covProfile <- c(1, rep(0, numCov - 1)) # default to baseline profile
  }
  
  checkValuesBetween(forCovs, 1, numCov)
  if (!is.null(forCats)) {
    checkValuesBetween(forCats, 1, numCat)
  } else {
    forCats <- 1:numCat # Default to all categories
  }
  i = 1
  for (cov in forCovs) {
    change[i] <- getCenterScaledChange(change[i], output$centerScaleList, cov,
                                       covProfile)
    i <- i + 1
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
      dimnames(toReturn[[covariteName]]) <- list(testPoints, output$catNames,
                                                 NULL)
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
        betaVCFits <- array(rep(betaVCFits, each = 250), 
                           dim = c(250, dim(betaVCFits)[2], dim(betaVCFits)[3]))
      }
      covariateName <- output$colMapping[[cov]]

      # combine the deltaRA for the specified categories as an array
      dRA <- abind::abind(lapply(forCats, function(cat)
                                     getDeltaRA(cat - 1, betaVCFits, fit,
                                                sumExpFit, change=change[i])),
                          along=3)
      toReturn[[covariateName]] <- aperm(dRA, c(1, 3, 2))
      dimnames(toReturn[[covariateName]]) <- list(testPoints, 
                                                  output$catNames[forCats], 
                                                  NULL)
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
#' @return A matrix of \eqn{\alpha_p(t)} samples
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

  ret <- getHillDiversity(output$beta, basis, covProfile, 
                          output$XvartoXColMapping, l=l)
  rownames(ret) <- testPoints
  return(ret)
}

#' @title Calculate the Multiplicative Change in \eqn{\alpha}-diversity 
#' Statistic
#' @description Calculates the \eqn{\Delta_v\alpha_p(t)} statistic using the 
#' Hill diversity index from the \code{FunCZIDM} samples.
#' @param output (list) Output list from \code{FunCZIDM}
#' @param change (numeric vector) The amount change of the covariate.
#' @param forCovs (numeric vector) Indices of \code{covariates}
#' for which \eqn{\Delta_v\alpha_p(t)} to be calculated.
#' @param l (numeric, default=0) Parameter l the diversity weight parameter. 
#' If 0, the Shannon diversity index is calculated, if 1, the count diversity 
#' index is calculated. Closer to 1 gives more weight to rare categories. 
#' @param covProfile (numeric vector, default=NULL) Covariate profile. If NULL, 
#' the baseline profile is used (ie all zeros).
#' @return A list of matrix of \eqn{\Delta_v\alpha_p(t)} samples
#' @export 
calcDeltaAlphaDiv <- function(output, change, forCovs, l = 0, covProfile=NULL) {
  checkList(output, "output")
  numCov <- length(output$colMapping)
  if (!is.null(covProfile)) {
    checkVector(covProfile, "covProfile", numericOnly=TRUE)
    checkLength(covProfile, numCov, "covProfile")
    covProfile <- getCenterScaledCovProfile(covProfile, output$centerScaleList)
  } else { # covProfile is NULL
    covProfile <- c(1, rep(0, numCov - 1)) # default to baseline profile
  }
  checkValuesBetween(forCovs, 1, numCov, "forCovs")
  i = 1
  for (cov in forCovs) {
    change[i] <- getCenterScaledChange(change[i], output$centerScaleList, cov,
                                       covProfile)
    i <- i + 1
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
  startTV <- match(2, output$XvartoXColMapping) # First instance of the first 
                                                # covariate
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
      betaVCFits <- array(rep(betaVCFits, each = 250), 
                          dim = c(250, dim(betaVCFits)[2], dim(betaVCFits)[3]))
    }

    toReturn[[output$colMapping[[cov]]]] <- getDeltaHillDiversity(betaVCFits, 
                                                                  fit, 
                                                                  sumExpFit,
                                                                  l=l,
                                                               change=change[i])
    rownames(toReturn[[output$colMapping[[cov]]]]) <- testPoints
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
    if (grepl("prop|mean", key, ignore.case=TRUE)) {
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
      if (grepl("prop|mean", key, ignore.case=TRUE)) {
        # if prop or mean are in the key, then we must combine by averaging
        combinedOutput[[key]] <- combinedOutput[[key]] +
                                 (output[[key]]/numOutputs)
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

#' @title Get Proportion of 0 in 95% CI of \eqn{\beta_{jp}(t)}
#' @description Calculates the proportion of \eqn{\beta_{jp}(t)} functions that 
#' have the zero function in their 95% credible interval.
#' @param output (list) Output list from \code{FunCZIDM}.
#' @param varCovRange (numeric vector, default=NULL) Range of the varying 
#' covariate to evaluate the \eqn{\beta_{jp}(t)} functions. By default, the
#' range is the min and max of the varying covariate in the data.
#' @return A matrix with columns "category", "covariate", and "Proportion 0 in 
#' 95 CI"
#' @export 
getPropNullInCI <- function(output, 
                            varCovRange = NULL) {
  XvartoXColMapping = output$XvartoXColMapping
  colMapping = output$colMapping
  varyingCov = output$varyingCov
  catNames = output$catNames
  interiorKnots = output$interiorKnots 
  boundaryKnots = output$boundaryKnots
  basisFunc = output$basisFunc 
  df = output$df

  # getting the points to evaluate
  minVarCov <- ifelse(is.null(varCovRange), min(varyingCov), varCovRange[1])
  maxVarCov <- ifelse(is.null(varCovRange), max(varyingCov), varCovRange[2])
  testPoints <- seq(minVarCov,  maxVarCov, length.out = 250)

  # get the basis functions for the test points
  basis <- basisFunc(testPoints, df=df, intercept=FALSE,
                     knots=interiorKnots,
                     Boundary.knots=boundaryKnots)
  basis <- cbind(1, basis)
  numDf <- df + 1
  
  # get the zero functions for checking the null in CI
  zeroFunc <- testPoints*0
  
  startTV <- match(2, XvartoXColMapping) # First instance of the first covariate
  numCov <- length(colMapping) # number of covariates
  numCat <- length(catNames) # number of categories
  returnMat <- matrix(NA, nrow=numCat*(numCov-1), ncol=3)
  colnames(returnMat) <- c("category", "covariate", "Proportion 0 in 95 CI")
  for (cov in 2:numCov) {
    # if the next instance in the mapping is not the next covariate, then
    # the covariate does not have varying coefficients. So startTV == endTV.
    if (cov == XvartoXColMapping[startTV + 1]) {
      endTV <- startTV + numDf - 1
      betaVCFits <- getBetaVC(output$beta, c(startTV, endTV), basis)
    } else {
      endTV <- startTV
      betaVCFits <- output$beta[startTV, , , drop=FALSE]
      betaVCFits <- array(rep(betaVCFits, each = 250), 
                          dim = c(250, dim(betaVCFits)[2], dim(betaVCFits)[3]))
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
      returnMat[(cov-2)*numCat + catToCheck, ] <- c(category, covariate, 
                                                    propZeroInCI[catToCheck])
    }
    startTV <- endTV + 1 
  }

  return(returnMat)
}

#' @title Get Min/Max Multiplicative Change in Alpha Diversity
#' @description Finds the covariate profiles corresponding to the minimum and 
#' maximum mean multiplicative change in alpha diversity 
#' (\eqn{\Delta_v\alpha_p(t)}) for a given covariate, within intervals of the 
#' varying covariate.
#' @param output (list) Output list from \code{FunCZIDM}.
#' @param changeFrom (numeric) Value to change from for the covariate.
#' @param changeTo (numeric) Value to change to for the covariate.
#' @param forCov (integer) Index of the covariate for which to calculate the change.
#' @param covariates (data.frame) Covariate data.
#' @param varyingCov (numeric vector) Values of the varying covariate.
#' @param varyingCovBreaks (numeric vector) Breaks for intervals of the varying covariate.
#' @param l (numeric, default=0) Diversity parameter for the Hill diversity index.
#' @param reqSig (logical, default=FALSE) If TRUE, only consider intervals where the 95% CI does not include 1.
#' @return A list with elements \code{maxCov} and \code{minCov}, each a matrix of covariate profiles corresponding to the maximum and minimum mean multiplicative change in alpha diversity for each interval.
#' @export
getMinMaxDeltaAlphaDiv <- function(output, changeFrom, changeTo, forCov, 
                                   covariates, varyingCov, varyingCovBreaks,
                                   l = 0, reqSig = FALSE) {
  checkList(output, "output")
  numCov <- length(output$colMapping)
  checkDataFrame(covariates, "covariates")
  X <- model.matrix(~., data=covariates) 
  for (row in 1:nrow(X)) {
    covProfiles <- getCenterScaledCovProfile(X, output$centerScaleList)
  }
  covProfiles[, forCov] <- changeFrom

  checkValuesBetween(forCov, 1, numCov, "forCov")
  changeTo <- getCenterScaledChange(changeTo, output$centerScaleList, forCov,
                                       covProfiles[1,])

  if (l < 0 | l > 1)
    stop(paste("l must be in [0, 1)."))

  basis <- output$basisFunc(varyingCov, df=output$df, intercept=FALSE,
                            knots=output$interiorKnots,
                            Boundary.knots=output$boundaryKnots)
  basis <- cbind(1, basis) # Add intercept
  numDf <- output$df + 1 # + 1 for intercept

  fit <- getFitOfData2(output$beta, basis, covProfiles, output$XvartoXColMapping)
  sumExpFit <- getSumExpFit(fit)
  
  startTV <- match(forCov, output$XvartoXColMapping) # First instance of the cov 

  while (forCov > output$XvartoXColMapping[startTV]) {
    # If the cov is not current startTV, then the current startTV is not
    # not a covariate that is wanted to calculate deltaRA for.
    startTV <- startTV + 1
  }

  # if the next instance in the mapping is not the next covariate, then
  # the covariate does not have varying coefficients. So startTV == endTV.
  if (forCov == output$XvartoXColMapping[startTV + 1]) {
    endTV <- startTV + numDf - 1
    betaVCFits <- getBetaVC(output$beta, c(startTV, endTV), basis)
  } else {
    endTV <- startTV
    betaVCFits <- output$beta[startTV, , , drop=FALSE]
    betaVCFits <- array(rep(betaVCFits, each = 250), 
                        dim = c(250, dim(betaVCFits)[2], dim(betaVCFits)[3]))
  }

  deltaDivs <- getDeltaHillDivMeanAndCI(betaVCFits, fit, sumExpFit, l=l, 
                                        change=changeTo)
  
  minMat <- matrix(covariates[1,], nrow=length(varyingCovBreaks)-1, 
                   ncol=ncol(covariates), byrow=T)
  maxMat <- matrix(covariates[1,], nrow=length(varyingCovBreaks)-1,
                   ncol=ncol(covariates), byrow=T)
  for (i in 2:length(varyingCovBreaks)) {
    low <- varyingCovBreaks[i - 1]
    high <- varyingCovBreaks[i]
    inVaryingCovBreaks <- which((varyingCov >= low)&(varyingCov < high))
    vals <- deltaDivs[inVaryingCovBreaks, ]
    if (reqSig) {
      minIdx <- which((vals[,2] == min(vals[,2]))&
                      ((vals[,1] > 1)|(vals[,2] < 1)))
      maxIdx <- which((vals[,2] == max(vals[,2]))&
                      ((vals[,1] > 1)|(vals[,2] < 1)))
      if (length(minIdx) != 0){
        minIdx <- inVaryingCovBreaks[minIdx] 
        minMat[i - 1,] <- covariates[minIdx,]
      } else {
        minMat[i - 1,] <- NA
      }
      if (length(maxIdx) != 0){
        maxIdx <- inVaryingCovBreaks[maxIdx]
        maxMat[i - 1,] <- covariates[maxIdx,]
      } else {
        maxMat[i - 1,] <- NA
      }
    } else {
      minIdx <- which.min(vals[,2])
      maxIdx <- which.max(vals[,2])
      minIdx <- inVaryingCovBreaks[minIdx]
      minMat[i - 1,] <- covariates[minIdx,]
      maxIdx <- inVaryingCovBreaks[maxIdx]
      maxMat[i - 1,] <- covariates[maxIdx,]
    }
  }

  minMat[, forCov - 1] <- NA
  maxMat[, forCov - 1] <- NA
  return(list("maxCov" = maxMat, "minCov" = minMat))
}

#' @title Get Min/Max Multiplicative Change in Relative Abundance
#' @description Finds the covariate profiles corresponding to the minimum and maximum mean multiplicative change in relative abundance (\eqn{\Delta_vRA_{jp}(t)}) for a given covariate and category, within intervals of the varying covariate.
#' @param output (list) Output list from \code{FunCZIDM}.
#' @param changeFrom (numeric) Value to change from.
#' @param changeTo (numeric) Value to change to.
#' @param forCov (integer) Index of the covariate for which to calculate.
#' @param forCat (integer) Index of the category for which to calculate.
#' @param covariates (data.frame) Covariate data.
#' @param varyingCov (numeric vector) Values of the varying covariate.
#' @param varyingCovBreaks (numeric vector) Breaks for intervals of the varying covariate.
#' @param reqSig (logical, default=FALSE) If TRUE, only consider intervals where the 95% CI does not include 1.
#' @return A list with elements "maxCov" and "minCov", each a matrix of covariate profiles.
#' @export
getMinMaxDeltaRA <- function(output, changeFrom, changeTo, forCov, forCat,
                             covariates, varyingCov, varyingCovBreaks,
                             reqSig = FALSE) {
  checkList(output, "output")
  numCov <- length(output$colMapping)
  checkDataFrame(covariates, "covariates")
  X <- model.matrix(~., data=covariates)
  for (row in 1:nrow(X)) {
    covProfiles <- getCenterScaledCovProfile(X, output$centerScaleList)
  }
  covProfiles[, forCov] <- changeFrom

  checkValuesBetween(forCov, 1, numCov, "forCov")
  checkValuesBetween(forCat, 1, length(output$catNames), "forCat")
  change <- getCenterScaledChange(changeTo, output$centerScaleList, forCov, covProfiles[1,])

  basis <- output$basisFunc(varyingCov, df=output$df, intercept=FALSE,
                            knots=output$interiorKnots,
                            Boundary.knots=output$boundaryKnots)
  basis <- cbind(1, basis)
  numDf <- output$df + 1

  fit <- getFitOfData2(output$beta, basis, covProfiles, output$XvartoXColMapping)
  sumExpFit <- getSumExpFit(fit)

  startTV <- match(forCov, output$XvartoXColMapping)
  while (forCov > output$XvartoXColMapping[startTV]) {
    startTV <- startTV + 1
  }

  if (forCov == output$XvartoXColMapping[startTV + 1]) {
    endTV <- startTV + numDf - 1
    betaVCFits <- getBetaVC(output$beta, c(startTV, endTV), basis)
  } else {
    endTV <- startTV
    betaVCFits <- output$beta[startTV, , , drop=FALSE]
    betaVCFits <- array(rep(betaVCFits, each = nrow(basis)),
                        dim = c(nrow(basis), dim(betaVCFits)[2], dim(betaVCFits)[3]))
  }

  # Get mean and CI for the specified category
  deltaRA <- getDeltaRAMeanAndCI(forCat - 1, betaVCFits, fit, sumExpFit, change=change)

  minMat <- matrix(covariates[1,], nrow=length(varyingCovBreaks)-1, 
                   ncol=ncol(covariates), byrow=T)
  maxMat <- matrix(covariates[1,], nrow=length(varyingCovBreaks)-1,
                   ncol=ncol(covariates), byrow=T)
  for (i in 2:length(varyingCovBreaks)) {
    low <- varyingCovBreaks[i - 1]
    high <- varyingCovBreaks[i]
    inVaryingCovBreaks <- which((varyingCov >= low) & (varyingCov < high))
    vals <- deltaRA[inVaryingCovBreaks, , drop=FALSE]
    if (nrow(vals) == 0) next
    if (reqSig) {
      minIdx <- which((vals[,2] == min(vals[,2])) &
                      ((vals[,1] > 1) | (vals[,3] < 1)))
      maxIdx <- which((vals[,2] == max(vals[,2])) &
                      ((vals[,1] > 1) | (vals[,3] < 1)))
      if (length(minIdx) != 0){
        minIdx <- inVaryingCovBreaks[minIdx]
        minMat[i - 1,] <- X[minIdx,]
      }
      if (length(maxIdx) != 0){
        maxIdx <- inVaryingCovBreaks[maxIdx]
        maxMat[i - 1,] <- X[maxIdx,]
      }
    } else {
      minIdx <- which.min(vals[,2])
      maxIdx <- which.max(vals[,2])
      minIdx <- inVaryingCovBreaks[minIdx]
      minMat[i - 1,] <- X[minIdx,]
      maxIdx <- inVaryingCovBreaks[maxIdx]
      maxMat[i - 1,] <- X[maxIdx,]
    }
  }

  minMat[, forCov - 1] <- NA
  maxMat[, forCov - 1] <- NA
  return(list("maxCov" = maxMat, "minCov" = minMat))
}