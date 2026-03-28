# Author: Brody Erlandson
#
# This file contains the code to combine the outputs from multiple chains of
# caseStudy.R

library(FunCZIDM)
library(optparse)

optionsList = list(
  make_option(c("--outputFolder"), type="character", default="infantApp",
              help="Output folder")
)
optParser = OptionParser(option_list=optionsList, add_help_option=FALSE)
opt = parse_args(optParser)

files <- c("output18342.rds", "output33261.rds", "output67235.rds", 
           "output87808.rds")
files <- paste0(opt$outputFolder, "/", files)
for (file in files) {
  output <- readRDS(file)
  output$beta <- output$beta[,, seq(4501, 8500, by = 4)] # only keeping samples after burn-in
  output$RA <- output$RA[,, seq(4501, 8500, by = 4)] # only keeping samples after burn-in
  saveRDS(output, file=file, compress="xz")
}
rm(output)
gc()

combinedFile <- paste0(opt$outputFolder, "/plotCaseStudySamples.rds")
combineOutputs(files, saveToFile=TRUE, fileName=combinedFile, compress="xz")