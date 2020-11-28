#' testing various methods ability to find true coeffients
#'
#'
#'
#' D:2019-07-24
#'
#'
library(testthat)
library(latentMixed)
set.seed(1)
betas <- c(1, 60)
simulatedata <-simulation_GeneData(1,
                                   k = 10,
                                   n = 200,
                                   L = 10,
                                   tau = 0.0,
                                   mu0 =0.00,
                                   sigma = 0.1,
                                   betas =betas)
MixGeneObj <- SetupGeneMix('Y ~ 1', data = data.frame(Y=simulatedata$Y), X=simulatedata$X)
MixGeneObj1 <- mixedModelForwardBackward(MixGeneObj, simulatedata$markers)
MixGeneObj2 <- GeneForward(MixGeneObj, restricted = F)
MixGeneObj3 <- GeneBackward(MixGeneObj2, restricted = F)
test_that("testing simple coeff to find", {
  expect_equal(sort(MixGeneObj1$find), betas)
})

test_that("testing simple coeff to find", {
  expect_equal(sort(MixGeneObj3$find), betas)
})
