v.norm2 <- function(x) {
  return(sqrt(sum(x^2)))
}

greedyHierByLv <- function(fullQ, n, offset, depth = 0, branch = 2) {
  #print(paste(n,offset,depth))
  if (n == 1) {
    result <- list(serr = v.norm2(fullQ[,offset + 1])^2,
                      sinv = matrix(1.0,1,1),
                      sdist = matrix(1.0,1,1),
                      sq = matrix(c(offset,offset),1,2) #offset, sq are zero based
    )
    return(result)
  }

  Q <- fullQ[,(offset+1):(offset + n)]
  if (all(apply(Q,1,min) == apply(Q,1,max))) {
    mat <- matrix(1.0/n/n,n,n)
    result <- list(serr = v.norm2(Q[,1])^2,
                        sinv = mat,
                        sdist = matrix(1.0,1,1),
                        sq = matrix(c(offset,offset+n-1),1,2) #sq is zero-based
    )
    return(result)
  }

  if (n <= branch) {
    bound <- cbind(0:(n-1),1:n)
  } else {
    rem <- n %% branch
    step <- (n-rem) / branch
    swi <- (branch - rem) * step
    if (n-1 >= swi) {
      sep <- c(seq(0,swi-1,step),seq(swi,n-1, step + 1),n)
    } else {
      sep <- c(seq(0,swi-1,step),n)
    }
    bound <- cbind(sep[1:(length(sep)-1)],sep[2:length(sep)])
  }

  result <- list()
  for (c in 1:dim(bound)[1]) {
    result[[c]] <- greedyHierByLv(fullQ, bound[c,2] - bound[c,1], offset + bound[c,1], depth +1, branch)
  }

  invAuList <- lapply(result,function(x) colSums(x$sinv))
  ncol <- length(invAuList)
  invAu <- unlist(invAuList)
  k <- sum(invAu)

  m1 <- 0
  for (i in 1:ncol) {
    temp <- Q[,(bound[i,1]+1):bound[i,2]] #a fix for matrix collapsing to vector
    dim(temp) <- c(dim(Q)[1],length(invAuList[[i]]))
    m1  <- m1 + v.norm2(temp %*% invAuList[[i]]) ^2
  }
  m <- v.norm2(Q %*% invAu) ^2

  sumerr <- sum(sapply(result,function(x) x$serr))

  granu <- 100
  decay <- 1.0 / ( branch^(depth / 2.0))
  err1 <- seq(granu,1) ^ 2
  err2 <- seq(0,granu-1) ^ 2 *decay
  toterr <- 1.0/err1 * (sumerr - ((m-m1)*decay+m1) * err2 / (err1+err2*k))
  err = min(toterr) * granu^2

  #use zero based-index. To be consistent with python
  perc <- 1.0 - (min(which(toterr == min(toterr)))-1.) / granu
  #print(perc)
  DIAG <- as.matrix(bdiag(lapply(result, function(x) x$sinv)))
  dim(invAu) <- c(n,1)

  inv <- (1.0/perc)^2 * (DIAG - (1.0-perc)^2 / ( perc^2 + k * (1-perc)^2 ) * invAu %*% t(invAu))
  dist <- c(1.0-perc,perc*unlist(sapply(result, function(x) x$sdist)))
  sq <- rbind(c(offset, offset+n-1),do.call(rbind, lapply(result, function(x) x$sq)))

  result <- list(serr = err,
                 sinv = inv,
                 sdist = dist,
                 sq = sq #offset, sq are zero based
  )
  return(result)
}

dawa <- function(Q,x,epsilon, ratio = 0.25, branch = 2) {
  #save(Q,x,epsilon,ratio,branch, file = "testdawa.RData")
  n <- length(x)
  x = x*1.0
  if (ratio > 0) {
    hist <- L1partition(x,epsilon,ratio,gethist = TRUE, approx = TRUE)
  } else {
    hist <- cbind(0:(n-1), 0:(n-1))
  }

  n2 <- dim(hist)[1]
  cum <- seq(0,dim(Q)[1] - 1, n)
  QtQ <- matrix(0,n2,n2)
  for (i in 1:length(cum)) {
    c0 <- cum[i]
    nrow <- min(dim(Q)[1],n)
    Q0mat <- matrix(0,nrow,n)
    for (c in 1: nrow) {
      wt <- Q[c + c0,1]
      lb <- Q[c + c0,2]
      rb <- Q[c + c0,3]
      Q0mat[c,(lb+1):(rb+1)] = wt #The input interval is [lb,rb] in Python or C
    }

    Qmat <- matrix(0,nrow,n2)
    for (c in 1: n2) {
      lb <- hist[c,1]
      rb <- hist[c,2]
      if (lb == rb) {
        Qmat[,c] <- Q0mat[,(lb+1):(rb+1)] #R column matrix will colapse to to row vector
      } else {
        Qmat[,c] <- apply(Q0mat[,(lb+1):(rb+1)],1,mean)
      }
    }

    QtQ = QtQ + t(Qmat) %*% Qmat
  }
  result <- greedyHierByLv(QtQ,n2,0,depth = 0, branch = branch)
  result$sinv <- result$sinv / ((1.0 - ratio)^2)

  qmat <- matrix(0,length(result$sdist),n2)
  y2 <- matrix(0,n2,1)
  count <-1
  for (c in 1:length(result$sdist)) {
    if (result$sdist[c] > 0) {
      lb <- result$sq[c,1]
      rb <- result$sq[c,2]
      qmat[count,(lb+1):(rb+1)] <- result$sdist[c] * (1.0 - ratio)
      y2[count] <- sum(x[(hist[lb+1,1]+1):(hist[rb+1,2]+1)])*result$sdist[c] * (1.0 - ratio)
      count <- count+ 1
    }
  }
  qmat <- qmat[1:count-1,]
  y2 <- y2 + rdoublex(length(y2), 0.0, 1.0 / epsilon)
  estv <- result$sinv %*% (t(qmat) %*% y2)

  estx <- matrix(0,1,n)
  for (c in 1:n2) {
    estx[(hist[c,1]+1):(hist[c,2]+1)] <- estv[c] / (hist[c,2] - hist[c,1] + 1 )
  }
  return(estx)
}

Dawa <- function(X,epsilon) {
  l <- length(X)
  Q <- matrix(0,l,3)
  Q[,1] <- 1
  Q[,2] <- 0
  Q[,3] <- seq(0,l-1)
  A <- dawa(Q,X,epsilon,0.2)
  return(A)
}


