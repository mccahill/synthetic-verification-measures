
L1partition <- function(x, epsilon,ratio = 0.5, gethist = FALSE, approx = FALSE) {
  n <- length(x)
  if (approx) {
    hist <- l1partition_approx(x,epsilon,ratio,sample.int(500000,1))
  } else {
    hist <- l1partition(x,epsilon,ratio,sample.int(500000,1))
  }
  rb <- n
  if (gethist) {
    bucks <- NULL
    for (lb.index in 2:length(hist)) { # use while loop instead?
      bucks <- rbind(c(hist[lb.index],rb-1),bucks)
      rb <- hist[lb.index]
      if (hist[lb.index] == 0) {
        break;
      }
    }
    return(bucks) #return zero-based indexed in R
  } else {
    hatx <- matrix(0,nrow = n, ncol = 1)
    for (lb.index in 2:length(hist)) {
      lb <- hist[lb.index]
      hatx[(lb + 1):rb] <- max(0,sum(x[(lb+1):rb]) + rdoublex(1, 0.0, 1.0 / (epsilon * (1.0-ratio)) )) / (rb - lb)
      rb <- lb
      if (lb== 0) {
        break;
      }
    }
    return(hatx)
  }
}
