rm(list = ls())
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
plot(nplots)

nplots = PriRP(D,epsilon,unit,bounds = c(-16,16,-2,2))
plot(nplots)
