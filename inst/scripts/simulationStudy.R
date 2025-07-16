options(error = function() { traceback(); quit(status = 1) })

# Load necessary library
library(optparse)
library(FunCZIDM)

###############################################################################
#                              optparse options
###############################################################################

optionsList = list(
  make_option(c("-n", "--numUnique"), type="integer", default=50,
              help="Number of unique persons"),
  make_option(c("-a", "--numMeas"), type="integer", default=1, # TODO : delete
              help="Number of measurements per person"),
  make_option(c("-c", "--numCat"), type="integer", default=5,
              help="Number of categories"),
  make_option(c("-p", "--numFE"), type="integer", default=10,
              help="Number of functional covariates"),
  make_option(c("-l", "--totalCountLow"), type="integer", default=1500,#150,
              help="Low end of range of total counts"),
  make_option(c("-h", "--totalCountHigh"), type="integer", default=2000,
              help="High end of range of total counts"),
  make_option(c("-i", "--numIterations"), type="integer", default=10,
              help="Number of iterations"),
  make_option(c("--seed"), type="integer", default=66142,
              help="Seed for random number generator"),
  make_option(c("-f", "--pngFilePath"), type="character", default="plots/",
              help="File path to save plots"),
  make_option(c("-b", "--burnIn"), type="integer", default=0,
              help="Number of burn in iterations"),
  make_option(c("--propSigma"), type="numeric", default=.3,
              help="Proposal standard deviation for beta"),
  make_option(c("--toReturn"), type="character", default=c("tau", "lambda", "RA"), 
              help="Return values from function"),
  make_option(c("--saveTraceplots"), type="logical", default=TRUE,
              help="Save traceplots"),
  make_option(c("--VCStats"), type="logical", default=TRUE,
              help="Get statistics for varying coefficient function"),
  make_option(c("--outputToFile"), type="logical", default=TRUE,
              help="Output results to file"),
  make_option(c("--noZI"), type="logical", default=FALSE,
              help="Run model without zero inflation"),
  make_option(c("--noRE"), type="logical", default=FALSE,
              help="Run model without random effects"),
  make_option(c("--noZIRE"), type="logical", default=FALSE,
              help="Run model without zero inflation or random effects"),
  make_option(c("--thin"), type="integer", default=1,
              help="Thinning parameter")
)

optParser = OptionParser(option_list=optionsList, add_help_option=FALSE)
opt = parse_args(optParser)

source(system.file("scripts", "simHelper.R", package = "FunCZIDM"),
       chdir = TRUE)

###############################################################################
#                                   Set up
###############################################################################
# Set Seed
set.seed(opt$seed)

# Set dimensions
n <- opt$numUnique # Number of unique observations
c <- opt$numCat # Number of categories
p <- opt$numFE # Number of fixed effects
totalCount <- c(opt$totalCountLow, opt$totalCountHigh) # Range of total counts

# Set running parameters
iter <- opt$numIterations # Number of iterations

# File path to save results
dir.create(sprintf("%d", opt$seed), showWarnings = FALSE)
pngFilePath <- sprintf("%d/%s", opt$seed, opt$pngFilePath)

# Create Folders if not existing
dir.create(pngFilePath)

# Directing output to a file if specified
if (opt$outputToFile) {
  sinkFile <- sprintf("%d/output%d.csv", opt$seed, opt$seed)
  evalFile <- sprintf("%d/eval%d.csv", opt$seed, opt$seed)
  cat("metric,value\n", file=evalFile)
  sink(sinkFile)
}

# Print parameters
cat("------------------------------ Parameters ------------------------------------\n")
cat("Number of unique observations: ", n, "\n", sep = "")
cat("Number of categories: ", c, "\n", sep = "")
cat("Number of fixed effects: ", p, "\n", sep = "")
cat("Range of total counts: ", totalCount[1], ",", totalCount[2], "\n", sep = "")
cat("----------------------------- Running parameters -----------------------------\n")

###############################################################################
#                               Generate data
###############################################################################
generatedData <- generateData(n, c, p, totalCountRange=totalCount)
X <- generatedData[["covariates"]]
counts <- generatedData[["counts"]]
personID <- generatedData[["ids"]]
timePoints <- generatedData[["timePoints"]]
betaTrue <- generatedData[["betaIntercepts"]]
rTrue <- generatedData[["rTrue"]]
funcMat <- generatedData[["betaFuncMat"]]
RA <- generatedData[["trueRA"]]
nTotal <- length(timePoints)
funcList <- generatedData[["funcList"]]

if (!opt$outputToFile) {
  cat("Proportion zero inflated:", generatedData[["trueZI"]], "\n")
} else {
  cat("True proportion zero inflated", generatedData[["trueZI"]], 
      sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
}

if (!opt$outputToFile) {
  cat("Proportion of zeros in counts:", sum(counts == 0)/length(counts), "\n")
} else {
  cat("Proportion of zeros in counts", sum(counts == 0)/length(counts), 
      sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
}

###############################################################################
#                                 Run model
###############################################################################

#tvOnly <- TRUE
TVCols <- 1:(ncol(X)+1) # Include intercept 
#linearComponent <- TRUE
#tvList <- addTimeVaryingCoeff(X, timePoints, TVCols, df=numDf, tvOnly=tvOnly,
#                              linearComponent=linearComponent, withNames=TRUE)
#Xtv <- tvList$X
#basis <- tvList$basis
#
#randEffCols <- c(1)
#if (opt$noZI) {
#  sourceCpp(sprintf("%ssrc/ZIDMwRE_noZI.cpp", opt$codeLocation))
#} else if (opt$noRE) {
#  sourceCpp(sprintf("%ssrc/ZIDMwRE_noRE.cpp", opt$codeLocation))
#} else if (opt$noZIRE) {
#  sourceCpp(sprintf("%ssrc/ZIDMwRE_noZIRE.cpp", opt$codeLocation))
#} else {
#  sourceCpp(sprintf("%ssrc/ZIDMwRE.cpp", opt$codeLocation))
#}
#
## start timer
#output <- zidmWre(iter, counts, Xtv, personID, numObsPerID, n, burnIn=opt$burnIn, 
#                        betaPriorSigma=opt$priorSigma, betaPropSigma=opt$propSigma,
#                        basis=basis, TVCols=TVCols, toReturn=opt$toReturn, 
#                        ZIGrouped=opt$ZIGrouped, adjustProposals=opt$adjustProp,
#                        interceptIncluded=TRUE, thin=opt$thin)
print("Starting model run...")
startTime <- Sys.time()
output <- FunCZIDM(counts, X, personID, timePoints, iter=opt$numIterations, 
                   burnIn=opt$burnIn, thin=opt$thin, toReturn=opt$toReturn,
                   proposalVars=list("beta proposal sd"=opt$propSigma),
                   saveToFile=FALSE, saveNullInBetaCI = FALSE)
endTime <- Sys.time()
totalTime <- as.numeric(difftime(endTime, startTime, units = "secs"))
print("Model run finished.")

numDf <- output$df
tvList <- list(basisFunc = output$basis, interiorKnots = output$interiorKnots,
               boundaryKnots = output$boundaryKnots)
###############################################################################
#                               Get results
###############################################################################

if (!opt$outputToFile) {
  cat("Total time: ", totalTime, " seconds", "\n", sep = "")
} else {
  cat("Total time", totalTime, sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
}

if (!opt$outputToFile) {
  # Get Acceptance rates
  cat("------------------------------ Acceptance Rates ------------------------------------\n")
  cat("2.5%, 50%, 97.5% acceptance rate quantiles for beta:", 
      quantile(output$betaAcceptProp, c(.025, .5, .975)), "\n", sep = " ")
  cat("2.5%, 50%, 97.5% acceptance rate quantiles for r:", 
      quantile(output$rAcceptProp, c(.025, .5, .975)), "\n", sep = " ")

  # Get True RA distribution
  cat("------------------------------ True Prob ------------------------------------\n")
  cat("2.5%, 97.5% quantiles:", quantile(RA, c(.025, .975)), "\n", sep = " ")
  cat("Median: ", median(RA), "\n", sep = "")
  cat("Mean: ", mean(RA), "\n", sep = "")
} else {
  cat("2.5% acceptance rate quantile for beta", quantile(output$betaAcceptProp, .025), 
      sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
  cat("50% acceptance rate quantile for beta", quantile(output$betaAcceptProp, .5),
      sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
  cat("97.5% acceptance rate quantile for beta", quantile(output$betaAcceptProp, .975),
      sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
  cat("2.5% acceptance rate quantile for r", quantile(output$rAcceptProp, .025), 
      sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
  cat("50% acceptance rate quantile for r", quantile(output$rAcceptProp, .5),
      sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
  cat("97.5% acceptance rate quantile for r", quantile(output$rAcceptProp, .975),
      sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)

  if (!opt$noZI || !opt$noZIRE) {
    cat("Proportion of zero indicators", output$etaMeanPropZeros, sep = ",", 
        file=evalFile, append=TRUE)
  }
  cat("\n", file=evalFile, append=TRUE)
  cat("true prob 2.5% quantile", quantile(RA, .025), sep = ",", 
      file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
  cat("true prob 97.5% quantile", quantile(RA, .975), sep = ",", 
      file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
  cat("true prob median", median(RA), sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
  cat("true prob mean", mean(RA), sep = ",", file=evalFile, append=TRUE)
  cat("\n", file=evalFile, append=TRUE)
}


# Get RA catpure rate
if ("RA" %in% names(output)) {
  cat("------------------------------ Probability ------------------------------------\n")
  getCaptureRate(output$RA, RA, checkZeros = TRUE)
  getCaptureRate(output$RA, RA, .97, checkZeros = TRUE)
  getCaptureRate(output$RA, RA, .99, checkZeros = TRUE)
  getMAE(output$RA, RA)
  getAitchison(output$RA, RA)
}

if ("eta" %in% names(output)) {
  cat("------------------------------ Eta ------------------------------------\n")
  meanProportionZeros <- mean(output$eta == 0)
  cat("Mean proportion of zeros: ", meanProportionZeros, "\n", sep = "")
}

# Saving Traceplots
if (opt$saveTraceplots) {
  # check if betaTrue is the same size as output$beta
  if (dim(output$beta)[1] == dim(betaTrue)[1] && dim(output$beta)[2] == dim(betaTrue)[2]) {
    saveTraceplotGrid(output$beta, paste0(pngFilePath,"betaTraceplotGrid.pdf"),
                      trueValues = betaTrue)
  } else {
    if (c > 30) {
      catRandSample <- c(1:4, sample(5:c, 26))
    } else {
      catRandSample <- 1:c
    }
    saveTraceplotGrid(output$beta[1:5,catRandSample,], 
                      paste0(pngFilePath,"betaTraceplotGridIntercept.pdf"))
    saveTraceplotGrid(output$beta[6:10,catRandSample,], 
                      paste0(pngFilePath,"betaTraceplotGridCoeff.pdf"))
    if ("r" %in% names(output)) {
      saveTraceplotGrid(output$r[,catRandSample,], paste0(pngFilePath,"rTraceplotGrid.pdf"))
    }
  }
}
if ("RA" %in% names(output)) {
  obsRandSample <- sample(1:nTotal, 10)
  if (c > 15) {
    catRandSample <- sample(1:c, 15)
  } else {
    catRandSample <- 1:c
  }
  saveTraceplotGrid(output$RA[obsRandSample, catRandSample,], 
                    paste0(pngFilePath,"probTraceplotGrid.pdf"),
                    trueValues = RA[obsRandSample, catRandSample])
}

if (opt$VCStats) {
  useBs <- TRUE
  interceptCols <- 1:(numDf + 1)
  intercepts <- output$beta[interceptCols,, , drop=FALSE]
  minTimepoint <- min(timePoints)
  maxTimepoint <- max(timePoints)
  testPoints <- seq(minTimepoint, maxTimepoint, length.out = 100)
  bFunc <- tvList$basisFunc
  basisTP <- bFunc(testPoints, df=numDf, intercept=FALSE,
                   knots=tvList$interiorKnots,
                   Boundary.knots=tvList$boundaryKnots)
  basisTP <- cbind(1, basisTP)
  numDf <- numDf + 1

  intercepts <- FunCZIDM:::getBetaVC(intercepts, c(1, numDf), basisTP)
  betaTrueFunc <- matrix(0, length(testPoints), c)
  for (i in 1:c) {
    f <- funcList[[round(betaTrue[1,i])]]
    plusMinus <- betaTrue[1,i] - round(betaTrue[1,i])
    betaTrueFunc[,i] <- ifelse(plusMinus >= 0, 1, -1)*f(testPoints)
  }
  betaTrue <- betaTrueFunc
  interceptTrueRA <- exp(betaTrue)/rowSums(exp(betaTrue))

  numTVCov <- length(TVCols)
  startTV <- 1
  covariates <- c(1, rep(0, numTVCov - 1))
  fit <- FunCZIDM:::getFit(output$beta, basisTP, covariates)
  sumExpFit <- FunCZIDM:::getSumExpFit(fit)
  baselineRA <- FunCZIDM:::getRAFitsPrecalc(fit, sumExpFit)
  outputFile <- sprintf("%d/RAStats%d.csv", opt$seed, opt$seed)
  cat("category,covariate,mean RMSE,sd RMSE,RA 95 CI,zero prop,B(t) 95 CI,accept rate\n",
      file=outputFile)
  
  for (tvc in 1:numTVCov) {
    betaVCFits <- FunCZIDM:::getBetaVC(output$beta, c(startTV, startTV+(numDf-1)), basisTP)
    covariate <- ifelse(startTV == 1, "intercept", paste0("cov", tvc-1))
    png(paste0(pngFilePath, gsub("[<>%]", "", gsub(" ", "_", paste0(covariate, "CIs.png")))),
        width=1500, height=4800)
    par(mfrow=c(10, 2), mar=c(4,4,2,2), oma=c(1,1,1,1))# first 10 beta functions
    
    if (opt$saveTraceplots) {
      RAtracePlots <- matrix(0, floor((iter - opt$burnIn)/opt$thin), 10)
      BetaTracePolts <- matrix(0, floor((iter - opt$burnIn)/opt$thin), 10)
    }
    for (catToCheck in 1:c) {
      betaVCFit = betaVCFits[,catToCheck,]
      if (covariate != "intercept") {
        trueRA <- FunCZIDM:::getTrueVCRAFromMat(betaTrue, testPoints, catToCheck,
                                     funcMat[tvc-1,])
        RA <- FunCZIDM:::getDeltaRA(catToCheck - 1, betaVCFits, fit, sumExpFit) 

        # get 2.5th and 97.5th percentiles of beta function 
        betaFuncQuant <- apply(betaVCFit, 1,
                               quantile, probs=c(0.025, 0.975))
        f <- funcList[[round(funcMat[tvc-1,catToCheck])]]
        plusMinus <- funcMat[tvc-1,catToCheck] - round(funcMat[tvc-1,catToCheck])
        trueBetaFunc <- ifelse(plusMinus >= 0, 1, -1)*f(testPoints)
        inCIBt <- betaFuncQuant[1,] < trueBetaFunc & 
                  trueBetaFunc < betaFuncQuant[2,]
      } else {
        RA <- baselineRA[,catToCheck,]
        trueRA <-  interceptTrueRA[, catToCheck]

        betaFuncQuant <- apply(intercepts[,catToCheck,], 1, quantile, probs=c(0.025, 0.975))
        inCIBt <- betaFuncQuant[1,] < betaTrue[,catToCheck] & 
                   betaTrue[,catToCheck] < betaFuncQuant[2,]
      }
      # get 2.5th and 97.5th percentiles of each row in RA
      RAQuant <- apply(RA, 1, quantile, probs=c(0.025, 0.975))
      # check how many of the trueRA values are within the 95% CI
      inCI <- RAQuant[1,] < trueRA & trueRA < RAQuant[2,]
      # get average RMSE
      RADiff <- apply(RA, 2, function(a) sqrt(mean((a-trueRA)^2)))
      # proportion of zero counts in counts
      zeroProp <- sum(counts[,catToCheck] == 0)/length(counts[,catToCheck])
      
      cat(catToCheck, covariate, mean(RADiff), sd(RADiff), mean(inCI), zeroProp, mean(inCIBt),
          output$betaAcceptProp[startTV:(startTV+numDf-1),catToCheck], sep = ",", file=outputFile, append=TRUE)
      cat("\n", file=outputFile, append=TRUE)
      
      if (catToCheck <= 10) {  
        # plot the 95% CI
        x <- testPoints
        yMax <- max(c(RAQuant[2,], trueRA))
        yMin <- min(c(RAQuant[1,], trueRA))
        if (yMax - yMin < 100) {
          yTicks <- round(seq(yMin, yMax, by=.1), 3)
        } else {
          yTicks <- c(NA)
        }
        
        RAMean <- apply(RA, 1, mean)
        if (covariate != "intercept") {
          ylab <- "Multipicative Change in Relative Abundance"
        } else {
          ylab <- "Relative Abundance"
        }
        plot(x, trueRA, type="l", col="red", xlab="Time", ylab=ylab,
             ylim=c(yMin, yMax), yaxt="n")
        if (!is.na(yTicks[1])) {
          axis(2, at=yTicks)
        }
        lines(x, RAMean, col="black")
        lines(x, RAQuant[1,], col="blue")
        lines(x, RAQuant[2,], col="blue")
        if (opt$saveTraceplots) {
          RAtracePlots[,catToCheck] <- RA[50, ]
        }

        # plot the beta function 95% CI
        if (covariate != "intercept") {
          yMax <- max(c(betaFuncQuant[2,], trueBetaFunc))
          yMin <- min(c(betaFuncQuant[1,], trueBetaFunc))
          if (yMax - yMin < 100) {
            yTicks <- round(seq(yMin, yMax, by=.1), 3)
          } else {
            yTicks <- c(NA)
          }
          betaFuncMean <- apply(betaVCFit, 1, mean)
          plot(x, trueBetaFunc, type="l", col="red", xlab="Time",
               ylab="Beta", ylim=c(yMin, yMax), yaxt="n")
          if (!is.na(yTicks[1])) {
            axis(2, at=yTicks)
          }
          lines(x, betaFuncMean, col="black")
          lines(x, betaFuncQuant[1,], col="blue")
          lines(x, betaFuncQuant[2,], col="blue")
          if (opt$saveTraceplots) {
            BetaTracePolts[,catToCheck] <- betaVCFit[50, ]
          }
        } else {
          yMax <- max(c(betaFuncQuant[2,], betaTrue[,catToCheck]))
          yMin <- min(c(betaFuncQuant[1,], betaTrue[,catToCheck]))
          if (yMax - yMin < 100) {
            yTicks <- round(seq(yMin, yMax, by=.1), 3)
          } else {
            yTicks <- c(NA)
          }
          betaFuncMean <- apply(intercepts[,catToCheck,], 1, mean)
          y <- betaTrue[,catToCheck]
          plot(x, y, type="l", col="red", xlab="Time",
               ylab="Beta", ylim=c(yMin, yMax), yaxt="n")
          if (!is.na(yTicks[1])) {
            axis(2, at=yTicks)
          }
          lines(x, betaFuncMean, col="black")
          lines(x, betaFuncQuant[1,], col="blue")
          lines(x, betaFuncQuant[2,], col="blue")
          if (opt$saveTraceplots) {
            BetaTracePolts[,catToCheck] <- intercepts[50, catToCheck,]
          }
        }
      }
    }
    dev.off()
    # Plot traceplots of delta RA(t) and B(t) at t=5
    if (opt$saveTraceplots) {
      png(paste0(pngFilePath, gsub("[<>%]", "", gsub(" ", "_", paste0(covariate, "Traceplots.png")))),
          width=1500, height=4800)
      par(mfrow=c(10, 2), mar=c(4,4,2,2), oma=c(1,1,1,1))# first 10 beta functions
      x <- 1:dim(RAtracePlots)[1]
      yMinRA <- min(RAtracePlots)
      yMaxRA <- max(RAtracePlots)
      yMinB <- min(BetaTracePolts)
      yMaxB <- max(BetaTracePolts)
      for (cat in 1:10) {
        plot(x, RAtracePlots[,cat], type="l", col="black", xlab="",
             ylab="Delta RA", ylim=c(yMinRA, yMaxRA))
        plot(x, BetaTracePolts[,cat], type="l", col="black", xlab="",
             ylab="beta", ylim=c(yMinB, yMaxB))
      }
    }
    dev.off()
    startTV <- startTV + numDf
  }
}
