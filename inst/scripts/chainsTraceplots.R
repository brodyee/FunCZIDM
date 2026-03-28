# Author: Brody Erlandson
#
# This file contains the code to make traceplots of the chains

library(FunCZIDM)
library(optparse)

optionsList = list(
  make_option(c("--outputFolder"), type="character", default=".",
              help="Output folder")
)
optParser = OptionParser(option_list=optionsList, add_help_option=FALSE)
opt = parse_args(optParser)

files <- c(paste0(opt$outputFolder, "/output18342.rds"), 
           paste0(opt$outputFolder, "/output33261.rds"), 
           paste0(opt$outputFolder, "/output67235.rds"), 
           paste0(opt$outputFolder, "/output87808.rds"))
taxa <- c("Bacilli", "Clostridia", "Gammaproteobacteria")

covProfile <- c(1, 27, 0, 0, 0, 0, 0, 0, .3)
change <- c(2, 1, 1) 
forCovs <- c(2, 7, 8) 
covNames <- c("Gestational Age", "Diet of > 50\% Breast Milk", "Diet of 10-50\% Breast Milk")
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
  png(paste0(opt$outputFolder, "/traceplot_", covNames[c], ".png"), units="in", width=8, height=6, res = 300)
  par(mfrow=c(length(taxa), 1), mar=c(4, 4, 2, 1))
  minYs <- numeric(length(taxa))
  maxYs <- numeric(length(taxa))
  for (t in 1:length(taxa)) {
    minYs[t] <- min(deltaRA[[1]][[c]][109, t, ])
    maxYs[t] <- max(deltaRA[[1]][[c]][109, t, ])
    for (f in 2:length(files)) {
      minYs[t] <- min(minYs[t], min(deltaRA[[f]][[c]][109, t, ]))
      maxYs[t] <- max(maxYs[t], max(deltaRA[[f]][[c]][109, t, ]))
    }
  }
  for (t in 1:length(taxa)) {
    plot(deltaRA[[1]][[c]][109, t, ], type="l", col="black", 
         xlab="Iteration", ylab="Change in RA sample",
         main=paste("Traceplot for", taxa[t], "and", covNames[c]),
         ylim=c(minYs[t], maxYs[t]))
    for (f in 2:length(files)) {
      lines(deltaRA[[f]][[c]][109, t, ], col=f)
    }
  }
  dev.off()
} 