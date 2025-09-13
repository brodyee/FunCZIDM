# Author: Brody Erlandson
#
# This file contains the code to combine the outputs from multiple chains of
# caseStudy.R

library(FunCZIDM)

files <- c("output18342.rds", "output33261.rds", "output67235.rds", 
           "output87808.rds")
for (file in files) {
  output <- readRDS(paste0("caseStudy/", file))
  output$beta <- output$beta[,, 7501:8500] # only keeping samples after burn-in
  saveRDS(output, file=file, compress="xz")
}
rm(output)
gc()

combinedFile <- "caseStudySamples.rds"
combineOutputs(files, saveToFile=TRUE, fileName=combinedFile, compress="xz")