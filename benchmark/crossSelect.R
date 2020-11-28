##
#'
#'
#'
#'
#'
#'
#'

n  <- 4
N  <- 5
set.seed(1)
X1 <- matrix(rnorm(n*N),N,n)
X2 <- matrix(rnorm(n*N),N,n)
Z1 <- matrix(as.integer(runif(n*N)<0.5),N,n)
Z2 <- matrix(as.integer(runif(n*N)<0.5),N,n)
mixing_population(X1,
                  X2,
                  Z1,
                  Z2,
                  2,
                  1)
