% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PriRP.R
\name{PriRP}
\alias{PriRP}
\title{Compute private regression residuals}
\usage{
PriRP(D, epsilon, unit, ratio = 0.3, bounds = NULL)
}
\arguments{
\item{D}{n by 2 matrix, with row i = (y_i, r_i) where y_i and r_i are the predicted value and residual of record i}

\item{epsilon}{privacy budget}

\item{unit}{unit - unit of the data (e.g. for "age", unit is 1; for "salary", unit is 1000 or 10000)}

\item{ratio}{budget epsilon*ratio spent on computing bounds, default is 0.3}

\item{bounds}{a vector of length 4 (x1, x2, y1, y2) for user defined bounds. Default to NULL and the bounds will be computed automatically}
}
\value{
nplots a list of noisy plots
}
\description{
Compute private regression residuals
}
\examples{
library(VerificationMeasures)
library(Matrix)
library(smoothmest)
options(expressions = 10000)
example.file <- system.file("extdata", "CPS_log.txt", package = "VerificationMeasures")
data <- read.table(example.file,header = FALSE,sep = " ")
D <- as.matrix(data)
epsilon <- 1.0
unit <- 1.0
#ratio <- 0.3 #default
nplots = PriRP(D,epsilon,unit)
#print(nplots)
#plot(nplots)
}

