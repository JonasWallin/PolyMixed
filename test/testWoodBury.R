
library(testthat)
library(latentMixed)
set.seed(1)
m <- 3
n <- 5
D <- 1+runif(n)
U <- matrix(rnorm(n=n*m), nrow= n, ncol = m)

Res <- solve(diag(D) + U%*%t(U))

Res_cpp <- WoodburySpecial(D,U)
context("Rcpp coding")
test_that("testing woodbury matrix inverse", {
  expect_equal(Res, Res_cpp)
})


test_that("testing cholesky update", {
  R <- CholUpdate(sqrt(diag(D)), as.matrix(U[,1]))
  expect_equal(R, chol(diag(D) + U[,1]%*%t(U[,1])))
})

A <- matrix(rnorm(n * n), nrow=n, ncol=n)
A <- A%*%t(A) + diag(n)
Ra <- chol(A)
test_that("testing cholesky update 2", {
  R <- CholUpdate(Ra, as.matrix(U[,1]))
  expect_equal(R, chol(A + U[,1]%*%t(U[,1])))
})
test_that("testing cholesky update multi", {
  R <- CholUpdate(Ra, U)
  expect_equal(R, chol(A + U%*%t(U)))
})

