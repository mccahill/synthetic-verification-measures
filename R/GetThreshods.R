.divide <- function(hist,e,l,r) {
  hist <- sort(hist)

  med <- hist[length(hist) %/%2 + 1]
  newhist <- c(hist,l,r)
  noise.med <- med + GetNoiseSS(newhist,e)
  left <- NULL
  right <- NULL
  if (length(hist) > 0) {
    for (x in 1: length(hist)) {
      if (hist[x] < noise.med) {
        left <- c(left, hist[x])
      } else {
        if (hist[x] > noise.med) {
          right <- c(right,hist[x])
        }
      }
    }
  }
  ans <- -1
  if (l <= noise.med && noise.med <= r) {
    ans <- noise.med
  }
  return(list(left = left, right = right, ans = ans))
}

.preprocess <- function(hist,e,level,levelmax,lleft,rright) {
  if (level >= levelmax) {
    return(list(p1 = hist, p2 = NULL))
  }
  Allset <- NULL
  if (length(hist) == 0) {
    ans = -1
  } else {
    result <- .divide(hist,e,lleft, rright)
    left <- result$left
    right <- result$right
    ans <- result$ans
  }
  if (ans == -1) {
    q <- (lleft + rright) * 0.5
    Allset <- c(Allset,q)
    left <- NULL
    right <- NULL
    if (length(hist) > 0) {
      for (x in 1:length(hist)) {
        if (hist[x] < q) {
          left <- c(left,hist[x])
        } else {
          if (hist[x] > q) {
            right <- c(right,hist[x])
          }
        }
      }
    }
    lresult <- .preprocess(left,e,level+1,levelmax, lleft,q)
    rresult <- .preprocess(right,e,level+1,levelmax,q,rright)
  } else {
    Allset <- c(Allset,ans)
    if (length(left) >0) {
      lresult<- .preprocess(left,e,level+1,levelmax,lleft,ans)
    } else {
      lresult <- list(p1 =NULL, p2 = NULL)
    }
    if (length(right) >0) {
      rresult <- .preprocess(right,e,level+1, levelmax, ans,rright)
    } else {
      rresult <- list(p1 = NULL, p2 = NULL)
    }
  }

  Allset <- c(Allset,lresult$p2,rresult$p2)
  t<- c(lresult$p1,rresult$p1)
  return(list(p1 = t, p2 = Allset))
}

GetThresholds <- function(X,epsilon,levelmax) {
  eps <- epsilon / levelmax
  result <- .preprocess(X,eps,0,levelmax,0,1)
  Allset <- sort(result$p2)
  return(Allset)
}
