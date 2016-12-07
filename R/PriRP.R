#sample plots based on noisy cells
GenerateNoiseResiduals <- function(X,bounds,grids) {
  if (length(bounds) == 2) {
    x1 <- -bounds[1]
    y1 <- -bounds[2]
    x2 <- bounds[1]
    y2 <- bounds[2]
  } else {
    x1 <- bounds[1]
    x2 <- bounds[2]
    y1 <- bounds[3]
    y2 <- bounds[4]
  }

  result <- matrix(0,sum(as.integer(X)),2)
  print(dim(result))
  count <- 0
  x <- ((x2-x1) / grids)
  y <- ((y2-y1) / grids)
  for (i in 1:grids) {
    for (j in 1:grids) {
      N <- round(X[i,j])
      if (N > 0) {
        xmin <- x1 + x * (i-1)
        xmax <- xmin + x
        ymin <- y1 + y * (j-1)
        ymax <- ymin + y

        result[(count+1):(count+N),1] <- runif(N,xmin,xmax)
        result[(count+1):(count+N),2] <- runif(N,ymin,ymax)
        count <- count + N

      }
    }
  }
  return(result)
}

#' Compute private regression residuals
#'
#' @param D n by 2 matrix, with row i = (y_i, r_i) where y_i and r_i are the predicted value and residual of record i
#' @param epsilon privacy budget
#' @param unit unit - unit of the data (e.g. for "age", unit is 1; for "salary", unit is 1000 or 10000)
#' @param ratio budget epsilon*ratio spent on computing bounds, default is 0.3
#' @param bounds a vector of length 4 (x1, x2, y1, y2) for user defined bounds. Default to NULL and the bounds will be computed automatically
#'
#' @return
#' nplots a list of noisy plots
#' @export
#' @examples
#' library(VerificationMeasures)
#' library(Matrix)
#' library(smoothmest)
#' options(expressions = 10000)
#' example.file <- system.file("extdata", "CPS_log.txt", package = "VerificationMeasures")
#' data <- read.table(example.file,header = FALSE,sep = " ")
#' D <- as.matrix(data)
#' epsilon <- 1.0
#' unit <- 1.0
#' #ratio <- 0.3 #default
#' nplots = PriRP(D,epsilon,unit)
#' #print(nplots)
#' #plot(nplots)
PriRP <- function(D, epsilon, unit, ratio = 0.3, bounds = NULL) {
  if (is.null(bounds)) {
    bounds <- ComputeBounds(D,epsilon * ratio, unit)
    x1 <- -bounds[1]
    y1 <- -bounds[2]
    x2 <- bounds[1]
    y2 <- bounds[2]
  } else {
    x1 <- bounds[1]
    x2 <- bounds[2]
    y1 <- bounds[3]
    y2 <- bounds[4]
  }
  l = dim(D)[1]
  tot <- sum(D[,1] <= x2 & D[,1] >= x1 & D[,2] <= y2 & D[,2] >= y1)
  if (tot < 10) {
    warning("less than 10 data points in the user provided bounds.")
  }
  #compute the number of cells
  totgrid <- floor(0.95 * tot * epsilon / 10.0)
  grids = 1
  while ((grids ^2) < totgrid) {
    grids <- grids + 1
  }

  if (grids < 4) {
    grids <- 4
  }

  x <- ((x2-x1) / grids)
  y <- ((y2-y1) / grids)
  X <- matrix(0,grids, grids)
  for (i in 1:l) {
    if (D[i,1] <= x2 && D[i,1] >= x1 && D[i,2] <= y2 && D[i,2] >= y1) {
      xcoord <- trunc((D[i,1] - x1) / x)
      ycoord <- trunc((D[i,2] - y1) / y)

      if (xcoord < 0) {xcoord <- 0} #simply index logic later
      if (xcoord >= grids) {xcoord <- grids -1}
      if (ycoord < 0) {ycoord <- 0}
      if (ycoord >= grids) {ycoord <- grids -1}

      xcoord <- xcoord + 1
      ycoord <- ycoord + 1
      X[xcoord,ycoord] <- X[xcoord, ycoord] + 1
    }
  }

  # perturb the plots
  if (is.null(bounds)) {
    noiseX <- Dawa2D(X,epsilon*(1-ratio))
  }  else {
    noiseX <- Dawa2D(X,epsilon)
  }

  # sample plots
  nplots <- GenerateNoiseResiduals(noiseX,bounds,grids)
  return(nplots)

}
