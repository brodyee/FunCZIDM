###############################################################################
#                   Helper functions for plot and results 
############################################################################### 

getTraceplot <- function(sample, trueValue=NA) {
  plot(sample, type="l", xlab="", ylab="")

  if (!is.na(trueValue)) {
    abline(h=trueValue, col="red")
  }
}

saveTraceplotGrid <- function(samples, filename, width=1000, height=1000, 
                              trueValues=NA) {
  if (length(dim(samples)) == 3) {
    totalSamples <- dim(samples)[1]*dim(samples)[2]
    dims <- dim(samples)[1:2]
    doubleLoop <- TRUE
  } else if (length(dim(samples)) == 2) {
    totalSamples <- dim(samples)[1]
    dims <- dim(samples)[1]
    doubleLoop <- FALSE
  }
  NUM_COL <- 2
  NUM_ROW <- ceiling(totalSamples/NUM_COL)
  height <- height*ceiling((NUM_ROW/12))
  
  filename <- paste0(substr(filename, 1, nchar(filename)-3), "png")
  png(filename, width=width, height=height)

  par(mfrow=c(NUM_ROW, NUM_COL), mar=c(0, 2, 0, 0), oma=c(0,0,0,0), 
      mgp=c(0, 0, -.5), tck=-.01)
  if (doubleLoop) {
    for (j in 1:dims[2]) {
      for (i in 1:dims[1]) {
        if (length(trueValues) > 1) {
          getTraceplot(samples[i, j, ], trueValues[i, j])
        } else {
          getTraceplot(samples[i, j, ])
        }
      }
    }
  } else {
    for (i in 1:dims) {
      if (length(trueValues) > 1) {
        getTraceplot(samples[i, ], trueValues[i])
      } else {
        getTraceplot(samples[i, ])
      }
    }
  }

  dev.off()
}

getCaptureRate <- function(sample, trueValues, credLevel=.95, checkZeros=FALSE) {
  NUM_ROW <- dim(trueValues)[1]
  NUM_COL <- dim(trueValues)[2]
  captured <- 0
  credProbs <- c((1-credLevel)/2, 1-(1-credLevel)/2)

  nearZero <- 0
  for (i in 1:NUM_ROW) {
    for (j in 1:NUM_COL) {
      credIterval <- quantile(sample[i, j, ], probs = credProbs)
      if (credIterval[1] <= trueValues[i, j] 
          & trueValues[i, j] <= credIterval[2]) {
        captured <- captured + 1
      } else if (checkZeros) {
         if (trueValues[i, j] < .001) {
           nearZero <- nearZero + 1
           if (round(credIterval[1], 2) == 0) {
             captured <- captured + 1
           }
         }
      } 
    }
  }

  cat("------------------- ", credLevel*100,
      "% Credible Interval Capture Rate ------------------------\n",
       sep = "")
  cat("Capture rate: ", captured/(NUM_ROW*NUM_COL), "\n", sep = "")
  if (checkZeros) {
    cat("Number of values near zero that weren't captured: ", nearZero, "\n", sep = "")
    cat("Proprtion of total values near zero that weren't captured: ", 
        (nearZero)/(NUM_ROW*NUM_COL), "\n", sep = "")
  }

}

getMAE <- function(sample, trueValues) {
  NUM_ROW <- dim(trueValues)[1]
  NUM_COL <- dim(trueValues)[2]
  mae <- 0

  for (i in 1:NUM_ROW) {
    for (j in 1:NUM_COL) {
      mae <- mae + abs(trueValues[i, j] - mean(sample[i, j, ]))
    }
  }

  cat("-------------------------- Mean Absolute Error --------------------------------\n")
  cat("MAE: ", mae/(NUM_ROW*NUM_COL), "\n", sep = "")
}

getAitchison <- function(sample, trueValues) {
  NUM_ROW <- dim(trueValues)[1]
  NUM_COL <- dim(trueValues)[2]
  mad <- 0
  eps <- .0000001

  for (i in 1:NUM_ROW) {
    mad <- mad + FunCZIDM:::aitchison_distance(trueValues[i,] + eps, rowMeans(sample[i,,]) + eps)
  }

  cat("------------------------- Mean Aitchison Distance -----------------------------\n")
  cat("MAD: ", mad/NUM_ROW, "\n", sep = "")
}