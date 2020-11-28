####
#  testing simple mixed estimation model
#
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
set.seed(1)
sigma <- 0.01
tau <- 0.2
beta0 <- 0.1
u   <- 0.2
n   <- 300
p1  <- 11
p2  <- 2*n
X1  <- matrix(rnorm(n*p1), nrow = n, ncol=p1)
X2  <- matrix(rnorm(n*p2), nrow = n, ncol=p2)

beta1 <- u + tau  *rnorm(n = p1)
beta2 <- u + tau  *rnorm(n = p2)

Y1 <-  beta0 + X1%*%beta1 + sigma * rnorm(n = n)
Y2 <-  beta0 + X2%*%beta2 + sigma * rnorm(n = n)


MixGeneObj1 <- SetupGeneMix('Y ~ 1', data = data.frame(Y=Y1), X=X1)
MixGeneObj2 <- SetupGeneMix('Y ~ 1', data = data.frame(Y=Y2), X=X2)
MixGeneObj1 <- estimateGeneMix(MixGeneObj1)
MixGeneObj2 <- estimateGeneMix(MixGeneObj2)

context("mixed model estimation")
test_that("testing mixed model with large and small p model", {

expect_equal(as.vector(MixGeneObj1$beta), as.vector(c(beta0,u)), tolerance=1e-1)
expect_equal(as.vector(MixGeneObj2$beta), as.vector(c(beta0,u)), tolerance=1e-1)
expect_equal(as.vector(c(MixGeneObj1$sigma,MixGeneObj1$tau)), as.vector(c(sigma,tau)), tolerance=1e-1)
expect_equal(as.vector(c(MixGeneObj2$sigma,MixGeneObj2$tau)), as.vector(c(sigma,tau)), tolerance=1e-1)
})
