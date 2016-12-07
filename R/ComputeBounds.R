#D is a l by p matrix, with each column is a variable
#return a row vector for each column in D

ComputeBounds <- function(D,eps,unit,alpha = 0.95) {

  l <- dim(D)[1]
  p <- dim(D)[2]
  Epsilon <- eps / 2

  D.abs <- abs(D)
  output <- matrix(0,1,p)
  for (i in 1:p) {
    x <- unit
    th <-  rdoublex(1, l * alpha, 2.0 / Epsilon)
    while (TRUE) {
      cnt <- sum(D.abs[,i] < x)
      ncnt <- rdoublex(1, cnt, 4.0 / Epsilon)
      if (ncnt >= th) {
        break #hate this.
      } else {
        x <- x * 2
      }
    }
    output[i] <- x
  }
  return(output)
}
