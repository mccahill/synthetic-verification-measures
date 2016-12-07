ComputeSS <- function(X,epsilon) {
  X <- sort(X)
  beta <- epsilon * 0.5
  l <- length(X)
  if (l <= 2) {
    return (X[2]-X[1])
  }
  m <- l %/% 2
  n <- l - m - 2
  SS <- 0
  for (k in 0:n) {
    cur <- 0
    for (t in 0:(k+1)) {
      temp <- (X[m+t+1] - X[m + t - k])
      if (temp > cur) {
        cur <- temp
      }
    }
    SS <- max(SS,cur / exp(k*beta))
  }
  return(SS)
}

GetNoiseSS <- function(X,epsilon) {
  u <- 0
  while (u == 0) {
    u <- runif(1)
  }
  u <- u - 0.5
  z <- tan( pi * u)
  noise <- ComputeSS(X,epsilon) / (epsilon / 8.0) * z
  return(noise)
}
