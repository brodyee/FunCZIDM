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
              help="number degrees of freedom")
)
optParser = OptionParser(option_list=optionsList, add_help_option=FALSE)
opt = parse_args(optParser)

###############################################################################
#                              Running model
###############################################################################

# Set seed for repoducability
set.seed(opt$seed)

# loading the data
data(infantCounts)
data(infantCovariates)

# getting the index of the subject id's and the time-varying variable
idIdx <- which(colnames(infantCovariates) == "Subject")
timeIdx <- which(colnames(infantCovariates) == "Day.of.life.sample.obtained")
ids <- infantCovariates[, idIdx]
time <- infantCovariates[, timeIdx]
infantCovariates <- infantCovariates[, -c(idIdx, timeIdx), drop = FALSE]

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
         toReturn=opt$toReturn,
         priors=list("first beta sd"=1, "a"=3, "b"=9, "alpha"=.01, 
                     "beta"=10),
         saveToFile=TRUE,
         fileName=outputFile,
         saveNullInBetaCI=FALSE) # will get null proportions after combining

