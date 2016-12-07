Data.Preprocess <- function(X,epsilon,ratio) {
  l <- dim(X)[1]
  tsize <- floor(log(l,base = 2))
  if (tsize > 10) {
    tsize <- 10
  }
  th <- X[,2]
  #selecting thresholds
  thresholds <- GetThresholds(th,epsilon * ratio,tsize)
  thresholds <- sort(thresholds, decreasing = TRUE)

  X.df <- as.data.frame(X)
  X <- as.matrix(X.df[rev(order(X.df$V2)),])

  acc <- NULL
  tot <- NULL

  num <- 0
  r <- 1
  curth <- thresholds[r]
  curacc <- 0
  curtot <- 0

  if (l > 0) {
    for (i in 1:l) {
      flag <- TRUE
      while (flag) {
        if (X[i,2] >= curth) {
          curacc <- curacc +  X[i,1]
          curtot <- curtot + 1
          flag <- FALSE
        } else {
          acc <- c(acc,curacc)
          tot <- c(tot,curtot)
          r <- r + 1
          if (r <= length(thresholds)) {
            curth <- thresholds[r]
          } else {
            curth <- 0
          }
          curacc <- 0
          next #do we really need this
        }
      }
    }
  }
  acc <- c(acc,curacc)
  tot <- c(tot, curtot)
  return(list(acc = acc, tot = tot))
}

# Compute noisy ROC curve
ComputeROC <- function(accNoise,accNoiseF) {
  size <- length(accNoise) # bad code. Rewrite later, after unit test
  NoiseTPR <- accNoise / accNoise[size]
  NoiseFPR <- accNoiseF / accNoiseF[size]
  NoiseROC <- cbind(NoiseTPR,NoiseFPR)
  return(NoiseROC)
}

Consistency <- function(NoiseROC) {
  TPR <- NoiseROC[,1]
  FPR <- NoiseROC[,2]
  size <- length(TPR)

  #monotonic consistency
  accT <- c(0,cumsum(TPR))
  accF <- c(0,cumsum(FPR))

  MT <- matrix(0,size,size) #only fill upper triangualr matrix
  MF <- matrix(0,size,size)
  for ( i in 1:size) {
    for (j in i:size) {
      MT[i,j] <- (accT[j+1] - accT[i]) / (j-i+1)
      MF[i,j] <- (accF[j+1] - accF[i]) / (j-i+1)
    }
  }

  ans <- matrix(0,size,2) #include both T and F. T in col 1, F in col 2

  for (k in 1:size) {
    curminT <- 10.0
    curminF <- 10.0
    for (j in k:size) {
      curT <- max(MT[1:j,j])
      curF <- max(MF[1:j,j])
      if (curT < curminT) {
        curminT <- curT
      }
      if (curF < curminF) {
        curminF <- curF
      }
    }
    if (curminT < 0) {
      curminT <- 0
    }
    if (curminT > 1) {
      curminT <- 1
    }
    if (curminF < 0) {
      curminF <- 0
    }
    if (curminF > 1) {
      curminF <- 1
    }
    ans[k,1] <- curminT
    ans[k,2] <- curminF
  }
  return(ans)
}

#compute AUC
Area <- function(NoiseROC) {
  NoiseA <- 0.0
  size <- dim(NoiseROC)[1]
  for (i in 1:(size-1)) {
    NoiseA <- NoiseA + (NoiseROC[i+1,2] - NoiseROC[i,2])*(NoiseROC[i+1,1] + NoiseROC[i,1])*0.5
  }
  return(NoiseA)
}

#' Compute differentially private ROC curve
#'
#' @param X n by 2 matrix, with row_i = (L_i,P_i) for the label and prediction of record_i
#' @param epsilon privacy budget
#' @param ratio budget epsilon*ratio spent on selecting thresholds, default is 0.2
#'
#' @return A list containing nAUC for noisy AUC and nROC for noisy ROC curve of (TPR,FPR) points
#' @export
#' @examples
#' library(VerificationMeasures)
#' library(Matrix)
#' library(smoothmest)
#' options(expressions = 10000)
#' example.file <- system.file("extdata", "ratios_twitter.csv", package = "VerificationMeasures")
#' data <- read.csv(example.file,header = FALSE)
#' X <- as.matrix(data[,c(1,2)]) #only use the first two columns?
#' epsilon <- 1
#' result <- PriROC(X,epsilon)
#' FPR <- result$nROC[,2]
#' TPR <- result$nROC[,1]
#' #plot(FPR, TPR, type = 'l', xlab = "False Positive Rate", ylab = "True Positive Rate")

PriROC <- function(X,epsilon,ratio = 0.2) {

  #compute "acc" and "tot" based on all thresholds

  dp <- Data.Preprocess(X,epsilon,ratio)
  acc <- dp$acc
  tot <- dp$tot

  #acc is TP, accf is FP
  accf <- tot - acc - c(0,tot[-length(tot)])

  # using DAWA to perturb TPs and FPs
  Noise <- Dawa(acc,epsilon*0.5*(1.0-ratio))
  NoiseF <- Dawa(accf,epsilon*0.5*(1.0-ratio))

  accNoise <- cumsum(Noise)
  accNoiseF <- cumsum(NoiseF)


  #compute ROC curve
  NoiseROC <- ComputeROC(accNoise,accNoiseF)
  NoiseROC <- rbind(c(0.0,0.0),NoiseROC)


  #make consistency
  nROC <- Consistency(NoiseROC)
  nAUC <- Area(nROC)

  return(list(nAUC = nAUC, nROC = nROC))
}
