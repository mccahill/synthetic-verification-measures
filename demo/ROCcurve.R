library(VerificationMeasures)

library(Matrix)
library(smoothmest)
options(expressions = 10000)
example.file <- system.file("extdata", "ratios_twitter.csv", package = "VerificationMeasures")
data <- read.csv(example.file,header = FALSE)
X <- as.matrix(data[,c(1,2)]) #only use the first two columns?

epsilon <- 1
result <- PriROC(X,epsilon)
FPR <- result$nROC[,2]
TPR <- result$nROC[,1]
plot(FPR, TPR, type = 'l', xlab = "False Positive Rate", ylab = "True Positive Rate")
