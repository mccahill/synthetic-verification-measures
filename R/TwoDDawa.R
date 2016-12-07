hilbert <- function(N) {
  if (N == 2) {
    return(list(x = c(0,0,1,1), y = c(0,1,1,0)))
  } else {
    r <- hilbert(N/2)
    x <- c(r$y, r$x, N/2 + r$x, N -1 - r$y)
    y <- c(r$x, N/2 + r$y, N/2 + r$y, N/2 -1 - r$x)
    return(list(x = x, y = y))
  }
}

Dawa2D <- function(X, epsilon, ratio = 0.5) {
  n1 <- dim(X)[1]
  n2 <- dim(X)[2]
  d = 2 ^ ceiling(log(max(n1,n2),base = 2))
  tmp <- matrix(0,d,d)
  tmp[1:n1,1:n2] <- X
  X <- tmp
  coords <- hilbert(d)

  #hilbert coordinates are zero-based, so need to be addjusted to one-based
  x1d <- unlist(mapply(function(r,c, M=M) M[r,c], coords$x + 1, coords$y + 1, MoreArgs = list(M = X)))

  Partition <- L1partition(x1d,epsilon,ratio,gethist = TRUE, approx = FALSE)
  nans <- matrix(0,1,length(x1d))
  for (i in 1:dim(Partition)[1]) {
    a <- Partition[i,1]
    b <- Partition[i,2]
    tot <- sum(x1d[(a+1):(b+1)])
    ntot <- rdoublex(1, tot, 1.0 / (epsilon*(1-ratio)))
    nans[(a+1):(b+1)] <- ntot / (b-a+1)
  }
  nans[nans<0] <-0
  nans <- as.integer(nans)
  hatx2d <- matrix(0,d,d)
  Index <- coords$y * d + coords$x + 1 #convert 2d indeices to 1d column-major, one-based index
  hatx2d[Index] <- nans
  hatx2d <- hatx2d[1:n1, 1:n2]
  return(hatx2d)
}


