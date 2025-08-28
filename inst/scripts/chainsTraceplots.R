# Author: Brody Erlandson
#
# This file contains the code to make traceplots of the chains

library(FunCZIDM)

files <- c("caseStudy/output18342.rds", "caseStudy/output33261.rds", 
           "caseStudy/output67235.rds", "caseStudy/output87808.rds")
taxa <- c("Bacilli", "Clostridia", "Gammaproteobacteria")

covProfile <- c(1, 27, 0, 0, 0, 0, 0, 0, .3)
change <- c(2, 1, 1) 
forCovs <- c(2, 7, 8) 
covNames <- c("GA", "milk50", "milk1050")
output <- readRDS(files[1])
forCats <- match(taxa, output$catNames)
deltaRA <- list()
deltaRA[[1]] <- calcDeltaRA(output, change, covProfile = covProfile, 
                       forCovs = forCovs, forCats = forCats)
rm(output)
gc()

for (f in 2:length(files)) {
  output <- readRDS(files[f])
  deltaRA[[f]] <- calcDeltaRA(output, change, covProfile = covProfile, 
                       forCovs = forCovs, forCats = forCats)
  rm(output)
  gc()
}

for (c in 1:length(forCovs)) {
  # make a plot for each covariate with all taxa in one plot
  png(paste0("traceplot_", covNames[c], ".png"), width=800, height=600)
  par(mfrow=c(length(taxa), 1), mar=c(4, 4, 2, 1))
  for (t in 1:length(taxa)) {
    plot(deltaRA[[1]][[c]][109, t, ], type="l", col="black", 
         xlab="Iteration", ylab="Change in RA sample",
         main=paste("Traceplot for", taxa[t], "and", covNames[c]))
    for (f in 2:length(files)) {
      lines(deltaRA[[f]][[c]][109, t, ], col=f)
    }
  }
  dev.off()
} 