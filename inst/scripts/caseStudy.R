# Author: Brody Erlandson
#
# This file contains the code to run the FunC-ZIDM regression model for 
# the infant data case study.
library(FunCZIDM)
library(optparse)

###############################################################################
#                              commandline options
###############################################################################

optionsList = list(
  make_option(c("-i", "--numIterations"), type="integer", default=100,
              help="Number of iterations"),
  make_option(c("--seed"), type="integer", default=87808,
              help="Seed for random number generator"),
  make_option(c("-b", "--burnIn"), type="integer", default=0,
              help="Number of burn in iterations"),
  make_option(c("--toReturn"), type="character", default=NULL, 
              help="Return values from function"),
  make_option(c("--outputFolder"), type="character", default="infantApp",
              help="Output folder"),
  make_option(c("--thin"), type="integer", default=5,
              help="Thinning"),
  make_option(c("--df"), type="integer", default=4,
              help="number degrees of freedom"),
  make_option(c("--kappaShape"), type="numeric", default=100.0,
              help="Shape parameter for beta prior on fixed effects"),
  make_option(c("--kappaRate"), type="numeric", default=900.0,
              help="Rate parameter for beta prior on fixed effects"),
  make_option(c("--propSigma"), type="numeric", default=.1,
              help="Proposal standard deviation for beta")
)
optParser = OptionParser(option_list=optionsList, add_help_option=FALSE)
opt = parse_args(optParser)

###############################################################################
#                              Running model
###############################################################################

# Set seed for repoducability
set.seed(opt$seed)

# loading the data
data(infantData)

# Convert the categorical variables to factors
infantCovariates$gender <- as.factor(infantCovariates$gender)
infantCovariates$mode.of.birth <- as.factor(infantCovariates$mode.of.birth)
infantCovariates$Room.category <- as.factor(infantCovariates$Room.category)
infantCovariates$milk <- as.factor(infantCovariates$milk)
infantCovariates$milk <- relevel(infantCovariates$milk, ref="< 10%")
infantCovariates$Period.of.study <- as.factor(infantCovariates$Period.of.study)
infantCovariates$Period.of.study <- relevel(infantCovariates$Period.of.study,
                                            ref="before")

# Check count and covariate data are aligned
if (!all(infantCounts$Subject == infantCovariates$Subject) &&
    !all(infantCounts$Day.of.life.sample.obtained 
         == infantCovariates$Day.of.life.sample.obtained)) {
  stop("The subject IDs and time points in the counts and covariates data do not
        match. Please check the data.")
}

# getting the index of the subject id's and the time-varying variable
idIdx <- which(colnames(infantCovariates) == "Subject")
timeIdx <- which(colnames(infantCovariates) == "Day.of.life.sample.obtained")
# extracting the subject ids and tv variable
ids <- infantCovariates[, idIdx]
time <- infantCovariates[, timeIdx]
infantCovariates <- infantCovariates[, -c(idIdx, timeIdx), drop = FALSE]

idIdx <- which(colnames(infantCounts) == "Subject")
timeIdx <- which(colnames(infantCounts) == "Day.of.life.sample.obtained")
infantCounts <- infantCounts[, -c(idIdx, timeIdx), drop = FALSE]
infantCounts <- as.matrix(infantCounts)

# setting file names and paths
dir.create(opt$outputFolder, showWarnings = FALSE)
outputFile <- sprintf("%s/output%d.rds", opt$outputFolder, opt$seed)

FunCZIDM(infantCounts, 
         infantCovariates, 
         ids,
         time,
         iter=opt$numIterations,
         burnIn=opt$burnIn,
         thin=opt$thin,
         returnBurnIn = TRUE,
         toReturn=opt$toReturn,
         priors=list("first beta sd"=1, "a"=3, "b"=9, "alpha"=.01, 
                     "beta"=10, "kappaShape"=opt$kappaShape,
                     "kappaRate"=opt$kappaRate),
         proposalVars=list("beta proposal sd"=opt$propSigma),
         covWithoutVC=c(which(grepl("Period", colnames(infantCovariates)))),
         df=opt$df,
         saveToFile=TRUE,
         fileName=outputFile,
         saveNullInBetaCI=FALSE) # will get null proportions after combining