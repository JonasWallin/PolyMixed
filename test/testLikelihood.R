###
#' testing the likelihood functions
#'
#'
#'
###

library(testthat)
library(latentMixed)

rawlikleihood <- function(theta, X, Xk, Y, likObj){
  tau     <- exp(theta[2])
  sigma_e <- exp(theta[1])
  Sigma <- tau^2 * (X%*%t(X))
  diag(Sigma) <- diag(Sigma) + sigma_e^2
  Xm <- likObj$Xu
  H <- t(Xk)%*%solve(Sigma, Xk)
  H <- t(Xk)%*%solve(Sigma, Xk)

  beta_hat =solve(H, (t(Xk)%*%(solve(Sigma, Y))))
  rm     <- Y - Xk%*%beta_hat
  loglik <- -0.5 * determinant(Sigma)$modulus[1]- 0.5 *(t(rm)%*%solve(Sigma, rm))
  loglik_res <- loglik + 0.5 * determinant(H)$modulus[1]
  BIC        <- loglik - 0.5 * determinant(H)$modulus[1]
return(list(loglik = loglik,
            loglik_res = loglik_res,
            BIC        = BIC))
}
n <- 100
p <-  101
sigma_e <- 0.1
tau     <- 0.2
X  <- matrix(rnorm(n*p), nrow = n, ncol = p)
Xk <- cbind(rep(1, n), as.matrix(rnorm(n), nrow=n))
Y <-  Xk%*%c(-2,3)  + sigma_e * rnorm(n) + X%*%rnorm(p, sd = tau)
MixGeneObj <- SetupGeneMix('Y ~ 1 + cov ',
                           data = data.frame(Y=Y, cov = Xk[,2]),
                           X=X,
                           meanMuOff  = T)
theta <- log(c(sigma_e, tau))
res <- rawlikleihood(theta, X, Xk, Y, MixGeneObj$LikObj)
loglik     <- res$loglik
loglik_res <- res$loglik_res
BIC        <- res$BIC
test_that("testing likelihood", {
  loglik_model <- loglikSimpleR(theta,
                           MixGeneObj$LikObj,
                           restricted = F)
  expect_equal(loglik_model, loglik)
})
test_that("testing likelihood restricted", {
  loglik_model <- loglikSimpleR(theta,
                                MixGeneObj$LikObj,
                                restricted = T)
  expect_equal(loglik_model, loglik_res)
})
test_that("testing likelihood BIC", {
  loglik_model <- loglikSimpleR(theta,
                                MixGeneObj$LikObj,
                                BIC = T)
  expect_equal(loglik_model, BIC)
})

library(numDeriv)
lambda <- function(theta){
  #'
  #' here we have parametrization theta = (\sigma^2 \tau^2) as parameters
  res <- rawlikleihood(theta, X, Xk, Y, MixGeneObj$LikObj)
  return(res$BIC)
}

H <- -hessian(lambda, theta)

test_that("testing likelihood BIC", {
  BICfull  <- BICmixed.R(theta,
                         MixGeneObj$LikObj)
  expect_equal(BICfull, -2*res$BIC+determinant(H)$modulus[1])
})

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
dupl <- findDuplicate(simulatedata$X)
MixGeneObj <- mixedModel(MixGeneObj,
                         dupl = dupl,
                         restricted=F)
Lik1 <- MixGeneObj$Lik
MixGeneObj_new <- mixedModel(MixGeneObj, find = 1, dupl = dupl, restricted=F)
Lik_old <- MixGeneObj_new$Lik
rawlikleihood(MixGeneObj_new$theta0, simulatedata$X, simulatedata$X[,1], simulatedata$Y, MixGeneObj$LikObj)

loglikSimpleR(MixGeneObj_new$theta0,
              MixGeneObj$LikObj,
              BIC = F)
