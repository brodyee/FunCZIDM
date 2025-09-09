## -----------------------------------------------------------------------------
library(FunCZIDM)

data(infantData)

## -----------------------------------------------------------------------------
# Check count and covariate data are aligned
(all(infantCounts$Subject == infantCovariates$Subject) ||
 all(infantCounts$Day.of.life.sample.obtained 
     == infantCovariates$Day.of.life.sample.obtained))

## -----------------------------------------------------------------------------
# Convert the categorical variables to factors
infantCovariates$gender <- as.factor(infantCovariates$gender)
infantCovariates$mode.of.birth <- as.factor(infantCovariates$mode.of.birth)
infantCovariates$Room.category <- as.factor(infantCovariates$Room.category)
infantCovariates$milk <- as.factor(infantCovariates$milk)
infantCovariates$milk <- relevel(infantCovariates$milk, ref="< 10%")
infantCovariates$Period.of.study <- as.factor(infantCovariates$Period.of.study)
infantCovariates$Period.of.study <- relevel(infantCovariates$Period.of.study,
                                            ref="before")

# Check the types of the covariates
str(infantCovariates)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
str(infantCovariates)

## -----------------------------------------------------------------------------
if (is.matrix(infantCounts)) {
  print("Counts data is a matrix.")
} 

colnames(infantCounts)

## ----eval=FALSE---------------------------------------------------------------
# # This code is not evaluated.
# FunCZIDM(infantCounts,
#          infantCovariates,
#          ids,
#          time,
#          iter=85000,
#          burnIn=75000,
#          thin=10,
#          toReturn=c("RA"),
#          priors=list("first beta sd"=1, "a"=3, "b"=3, "alpha"=.01,
#                      "beta"=10, "kappaShape"=100, "kappaRate"=900),
#          saveToFile=TRUE,
#          saveNullInBetaCI=FALSE)

## ----eval=FALSE---------------------------------------------------------------
# # This code is not evaluated.
# outputFiles <- c("output1.rds", "output2.rds", "output3.rds", "output4.rds")
# combinedOutput <- combineOutputs(outputFiles, saveToFile=FALSE)

## -----------------------------------------------------------------------------
combinedOutput <- readRDS(system.file("extdata", "caseStudySamples.rds",
                          package = "FunCZIDM"))
names(combinedOutput)

## -----------------------------------------------------------------------------
combinedOutput$colMapping

## -----------------------------------------------------------------------------
# always put 1 for the intercept
covProfile <- c(1, 27, 0, 0, 0, 0, 0, 0, .3)
# calculate the relative abundance
RA <- calcRA(combinedOutput, covProfile = covProfile)
dim(RA)

## -----------------------------------------------------------------------------
# calculate the change in relative abundance
change <- c(2, 1, 1) # two week change in gestational age and a 1 for indicating 
                  # the change from <10% to >50% breast milk in the diet
forCovs <- c(2, 7, 8) # GA is 2nd idx and >50% BM in diet is 6th idx
forCats <- match(c("Bacilli", "Clostridia", "Gammaproteobacteria"),
                 combinedOutput$catNames) # gets the idx of these taxa
deltaRA <- calcDeltaRA(combinedOutput, change, covProfile = covProfile,
                       forCovs = forCovs, forCats = forCats)
names(deltaRA)
colnames(deltaRA$Gestational.age.at.birth...weeks)
dim(deltaRA$Gestational.age.at.birth...weeks)

## -----------------------------------------------------------------------------
l <- 0.75 # Hill's diversity parameter
alphaDiv <- calcAlphaDiv(combinedOutput, l=l, covProfile=covProfile)
dim(alphaDiv)

## -----------------------------------------------------------------------------
# calculate the change in relative abundance
change <- c(2, 1) # two week change in gestational age and a 1 for indicating 
                  # the change from <10% to >50% breast milk in the diet
forCovs <- c(2, 6) # GA is 2nd idx and >50% BM in diet is 6th idx

deltaAlphaDiv <- calcDeltaAlphaDiv(combinedOutput, change, l=l, 
                                   covProfile=covProfile, forCovs=forCovs)
names(deltaAlphaDiv)
dim(deltaAlphaDiv$Gestational.age.at.birth...weeks)

## -----------------------------------------------------------------------------
betas <- getBetaFunctions(combinedOutput)
dim(betas)

## -----------------------------------------------------------------------------
betas <- getBetaFunctions(combinedOutput, cov=6)
dim(betas)

## -----------------------------------------------------------------------------
n <- 50 # Number of individuals
c <- 10 # Number of categories
p <- 10 # Number of functional covariates
dataList <- generateData(n, c, p)

## -----------------------------------------------------------------------------
cov <- 1 # the first covariate, not the intercept.
cat <- 1
f <- dataList$funcList[[round(dataList$betaFuncMat[cov, cat])]]
plusMinus <- dataList$betaFuncMat[cov, cat] - round(dataList$betaFuncMat[cov, cat])
betaFunc <- ifelse(plusMinus >= 0, 1, -1)*f(dataList$timePoints)

