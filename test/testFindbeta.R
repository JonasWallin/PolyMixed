####
# tries locating t-values
# D:2019-03-28
####

rm(list=ls())
library(latentMixed)
library(testthat)
###
# mixed model with intercept
#'   Y       \sim N(X \gamma, \sigma_Y^2 I_{n \times n})
#'   \gamma  \sim N(u, \tau^2 I_{p \times p})
###
set.seed(12)
sigma <- 0.01
tau <- 0.2
beta0 <- 0.1
u   <- 0.2
n   <- 300
p1  <- 11
p2  <- 2*n
markers1 <- c(ceiling(p1/2), p1 - ceiling(p1/2))
markers2 <- c(ceiling(p2/2), p2 - ceiling(p2/2))
X1  <- matrix(rnorm(n*p1), nrow = n, ncol=p1)
X2  <- matrix(rnorm(n*p2), nrow = n, ncol=p2)

beta1 <- u + tau  *rnorm(n = p1)
beta2 <- u + tau  *rnorm(n = p2)

Y1 <-  beta0 + X1%*%beta1 + sigma * rnorm(n = n)
Y2 <-  beta0 + X2%*%beta2 + sigma * rnorm(n = n)


MixGeneObj1 <- SetupGeneMix('Y ~ 1', data = data.frame(Y=Y1), X=X1)
MixGeneObj1b <- SetupGeneMix('Y ~ 1', data = data.frame(Y=Y1), X=X1, meanMuOff=T)
MixGeneObj1c <- SetupGeneMix('Y ~ -1', data = data.frame(Y=Y1), X=X1, meanMuOff=T)
MixGeneObj2 <- SetupGeneMix('Y ~ 1', data = data.frame(Y=Y2), X=X2)
MixGeneObj1b <- mixedModel(MixGeneObj1b)
MixGeneObj1c <- mixedModel(MixGeneObj1c)
MixGeneObj1 <- mixedModel(MixGeneObj1)
MixGeneObj2 <- mixedModel(MixGeneObj2)

MixGeneObj1b <- mixedModelForwardBackward(MixGeneObj1b, markers1)
MixGeneObj1c <- mixedModelForwardBackward(MixGeneObj1c, markers1)
MixGeneObj2 <- mixedModelForwardBackward(MixGeneObj2, markers2)
