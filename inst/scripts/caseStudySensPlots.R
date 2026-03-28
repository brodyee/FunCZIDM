# Author: Brody Erlandson
#
# Build comparative sensitivity plots for the case study prior scenarios.

library(reticulate)
library(FunCZIDM)

helpersFile <- system.file("app", "helpers.R", package = "FunCZIDM")
if (!nzchar(helpersFile)) {
  helpersFile <- "inst/app/helpers.R"
}
source(helpersFile)

plotModule <- system.file("scripts", "caseStudySensPlotFunctions.py",
                          package = "FunCZIDM")
if (!nzchar(plotModule)) {
  plotModule <- "inst/scripts/caseStudySensPlotFunctions.py"
}

pyPackages <- c("numpy", "pandas", "matplotlib", "seaborn")
missing <- pyPackages
tryCatch(
  expr = {
    reticulate::use_condaenv("shinyPy", required = TRUE)
    missing <- pyPackages[!vapply(pyPackages, reticulate::py_module_available,
                                   logical(1))]
    if (length(missing)) {
      stop("Missing Python modules in env shinyPy: ",
           paste(missing, collapse = ", "))
    }
  },
  error = function(e) {
    message("shinyPy conda env not available. Installing required Python packages in active environment.")
    message("Configured Python:")
    py_config()
    if (length(missing)) {
      py_require(missing)
    }
  }
)
reticulate::source_python(plotModule)

data(infantData, package = "FunCZIDM")
counts <- infantCounts
idIdx <- which(colnames(counts) == "Subject")
timeIdx <- which(colnames(counts) == "Day.of.life.sample.obtained")
counts <- counts[, -c(idIdx, timeIdx), drop = FALSE]
counts <- as.matrix(counts)
rowTotals <- rowSums(counts)
rowTotals[rowTotals == 0] <- 1
trueRA <- counts/rowTotals

gaCovariate <- "Gestational.age.at.birth...weeks"
gaExposureLabels <- c("GA +1 week", "GA -1 week")
milkCovariates <- c("milk> 50%", "milk10–50%")
milkExposureLabels <- c("Breast milk > 50% vs < 10%",
                        "Breast milk 10-50% vs < 10%")
aitchisonEps <- 1e-8

summarizeScenario <- function(output, scenario, taxa, exposureSettings, covProfile) {
  xRange <- c(min(output$varyingCov), 50) # max time 5 for inference
  scenarioRows <- list()
  for (idx in seq_len(nrow(exposureSettings))) {
    covName <- exposureSettings$covariate[idx]
    covLabel <- exposureSettings$label[idx]
    covChange <- exposureSettings$change[idx]
    deltaData <- deltaRAData(output, covName, taxa, covProfile, covChange,
                             xRange, betas = FALSE)
    for (tax in taxa) {
      if (is.null(deltaData[[tax]])) {
        stop("Missing delta RA data for ", tax, " in scenario ", scenario)
      }
      raFrame <- deltaData[[tax]][["RA"]]
      scenarioRows[[length(scenarioRows) + 1]] <- data.frame(
        scenario = scenario,
        taxa = tax,
        exposure = covLabel,
        time = raFrame$testpoints,
        mean = raFrame$mean,
        lower = raFrame$lowerCI,
        upper = raFrame$upperCI,
        stringsAsFactors = FALSE
      )
    }
    rm(deltaData)
    gc()
  }
  do.call(rbind, scenarioRows)
}

scenarioDirs <- c(
  "CS_Orig",
  "CS_unifTh",
  "CS_midTh",
  "CS_weakInfoKappa",
  "CS_highShrKappa",
  "CS_lowShrKappa",
  "CS_highShrPhi",
  "CS_lowShrPhi",
  "CS_highGlobParam",
  "CS_lowGlobParam",
  "CS_lowIntVar",
  "CS_highIntVar",
  "CS_vLowIntVar",
  "CS_vHighIntVar"
)

scenarioFiles <- file.path(scenarioDirs, "plotCaseStudySamples.rds")
names(scenarioFiles) <- scenarioDirs

missing <- scenarioFiles[!file.exists(scenarioFiles)]
if (length(missing)) {
  stop("Missing scenario outputs: ", paste(missing, collapse = ", "))
}

if (!"CS_Orig" %in% names(scenarioFiles)) {
  stop("CS_Orig scenario is required for comparisons")
}

origOutput <- readRDS(scenarioFiles["CS_Orig"])
gaIdx <- which(unlist(origOutput$colMapping) == gaCovariate)
milk1050Idx <- which(unlist(origOutput$colMapping) == milkCovariates[1])
milk50Idx <- which(unlist(origOutput$colMapping) == milkCovariates[2])
minVarCov <- min(origOutput$varyingCov)
maxVarCov <- max(origOutput$varyingCov)
testPoints <- seq(minVarCov,  maxVarCov, length.out = 250)
lastIdx <- max(which(testPoints < 50))
origGaBeta <- apply(FunCZIDM::getBetaFunctions(origOutput, cov = gaIdx), c(1,2), mean)[1:lastIdx,]
origMilk1050Beta <- apply(FunCZIDM::getBetaFunctions(origOutput, cov = milk1050Idx), c(1,2), mean)[1:lastIdx,]
origMilk50Beta <- apply(FunCZIDM::getBetaFunctions(origOutput, cov = milk50Idx), c(1,2), mean)[1:lastIdx,]

betaDiffSummary <- data.frame(
  scenario = names(scenarioFiles),
  meanAbsBetaDiffGA = NA,
  meanAbsBetaDiffMilk1050 = NA,
  meanAbsBetaDiffMilk50 = NA,
  stringsAsFactors = FALSE
)
aitchisonSummary <- data.frame(
  scenario = names(scenarioFiles),
  meanAitchisonDistance = NA,
  stringsAsFactors = FALSE
)

message("Generating sensitivity plot from ", length(scenarioFiles),
        " scenarios ...")

taxa <- c("Bacilli", "Clostridia", "Gammaproteobacteria")

gaEntries <- vector("list", length(scenarioFiles))
milkEntries <- vector("list", length(scenarioFiles))
for (i in seq_along(scenarioFiles)) {
  scenario <- names(scenarioFiles)[i]
  message("Processing ", scenario, " ...")
  if (scenario == "CS_Orig") {
    output <- origOutput
  } else {
    output <- readRDS(scenarioFiles[i])
    gaIdx <- which(unlist(output$colMapping) == gaCovariate)
  }
  covProfile <- c(1, rep(0, length(output$colMapping) - 1))
  covProfile <- getCenterScaledCovProfile(covProfile, output$centerScaleList)
  gaChanges <- c(
    "GA +1 week" = getCenterScaledChange(1, output$centerScaleList, gaIdx,
                                          covProfile),
    "GA -1 week" = getCenterScaledChange(-1, output$centerScaleList, gaIdx,
                                          covProfile)
  )
  gaSettings <- data.frame(
    covariate = rep(gaCovariate, length(gaChanges)),
    label = names(gaChanges),
    change = as.numeric(gaChanges),
    stringsAsFactors = FALSE
  )
  milkChanges <- vapply(milkCovariates, function(covName) {
    covIdx <- which(unlist(output$colMapping) == covName)
    if (!length(covIdx)) {
      stop("Unable to locate ", covName, " covariate in ", scenario)
    }
    getCenterScaledChange(1, output$centerScaleList, covIdx, covProfile)
  }, numeric(1))
  milkSettings <- data.frame(
    covariate = milkCovariates,
    label = milkExposureLabels,
    change = as.numeric(milkChanges),
    stringsAsFactors = FALSE
  )
  gaEntries[[i]] <- summarizeScenario(output, scenario, taxa, gaSettings,
                                      covProfile)
  milkEntries[[i]] <- summarizeScenario(output, scenario, taxa, milkSettings,
                                        covProfile)
  scenarioGaBeta <- apply(FunCZIDM::getBetaFunctions(output, cov = gaIdx), c(1,2), mean)[1:lastIdx,]
  scenarioMilk1050Beta <- apply(FunCZIDM::getBetaFunctions(output, cov = milk1050Idx), c(1,2), mean)[1:lastIdx,]
  scenarioMilk50Beta <- apply(FunCZIDM::getBetaFunctions(output, cov = milk50Idx), c(1,2), mean)[1:lastIdx,]
  betaDiffSummary$meanAbsBetaDiffGA[i] <- mean(abs(scenarioGaBeta - origGaBeta))
  betaDiffSummary$meanAbsBetaDiffMilk1050[i] <- mean(abs(scenarioMilk1050Beta - origMilk1050Beta))
  betaDiffSummary$meanAbsBetaDiffMilk50[i] <- mean(abs(scenarioMilk50Beta - origMilk50Beta))
  scenarioRAMean <- apply(output$RA, c(1, 2), mean)
  distances <- vapply(seq_len(nrow(trueRA)), function(rowIdx) {
    FunCZIDM:::aitchison_distance(trueRA[rowIdx, ] + aitchisonEps,
                                  scenarioRAMean[rowIdx, ] + aitchisonEps)
  }, numeric(1))
  aitchisonSummary$meanAitchisonDistance[i] <- mean(distances)
  if (scenario != "CS_Orig") {
    rm(output)
  }
  gc()
}

summaryDfGA <- do.call(rbind, gaEntries)
summaryDfMilk <- do.call(rbind, milkEntries)

labels <- c(
  "Default",
  "Uniform ZI Proportion",
  "Weaker ZI Proportion",
  "Weakly Informative Shrinkage",
  "Higher Shrinkage",
  "Lower Shrinkage",
  "Higher Individual Effect Shrinkage",
  "Lower Individual Effect Shrinkage",
  "Higher Global Shrinkage Parameter",
  "Lower Global Shrinkage Parameter",
  "Lower Intercept Variance",
  "Higher Intercept Variance",
  "Very Low Intercept Variance",
  "Very High Intercept Variance"
)
summaryDfGA$label <- NA
summaryDfMilk$label <- NA
for (s in 1:length(scenarioDirs)) {
  summaryDfGA$label[summaryDfGA$scenario == scenarioDirs[s]] <- labels[s]
  summaryDfMilk$label[summaryDfMilk$scenario == scenarioDirs[s]] <- labels[s]
}

plotPathGA <- "caseStudySens_deltaRA_GA.png"
createCaseStudySensPlot(summaryDfGA, plotPathGA,
                        width = 16, height = 10, dpi = 150)
message("Saved GA plot to ", normalizePath(plotPathGA, mustWork = FALSE))

plotPathMilk <- "caseStudySens_deltaRA_milk.png"
createCaseStudySensPlot(summaryDfMilk, plotPathMilk,
                        width = 16, height = 10, dpi = 150)
message("Saved breast milk plot to ", normalizePath(plotPathMilk, mustWork = FALSE))

cat("\nMean absolute beta-function difference vs CS_Orig (GA covariate):\n")
print(betaDiffSummary[,c(1,2)])

cat("\nMean absolute beta-function difference vs CS_Orig (Milk 10-50% covariate):\n")
print(betaDiffSummary[,c(1,3)])

cat("\nMean absolute beta-function difference vs CS_Orig (Milk 50% covariate):\n")
print(betaDiffSummary[,c(1,4)]) 

cat("\nMean Aitchison distance between observed RA and scenario mean RA:\n")
print(aitchisonSummary)
