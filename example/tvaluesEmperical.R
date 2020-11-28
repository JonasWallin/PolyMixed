##
# Emperical invstigation of t-values under H_0
#
# Conern the formula is invalid using L rather it should use 1-p...
# D: 2018-12-23
##

rm(list=ls())
set.seed(1)
source('simulateDesignMatrix.R')
library(latentMixed)
sim <- 100
model  = 'constant'
n <- 200
mu <- 0.004
sigma <- 1
tau   <- 0.01
d <- 1
m <- 150
chr <- 10
nc  <- NULL
if(model=='constant')
  nc = chr
simulate.X = T
markers <- rep(m, chr)   # markers on each chromosome
alpha = 0.05
C <- length(markers)
L <- sum(markers)
v <- function(y) (2 / y * (pnorm(y / 2) - 0.5)) / (y / 2 * pnorm(y / 2) + dnorm(y / 2))
crit.vec    <- matrix(0, nrow=sim, ncol = 2)
tvalues.vec <- matrix(0, nrow=sim, ncol = 2)
beta.vec    <- matrix(0, nrow=sim, ncol = 2)
select.vec  <- matrix(0, nrow=sim, ncol = 2)
for(j in 1:sim){
  cat('j = ',j,'\n')
  if(simulate.X ==T | j==1){
    X <- simulateDesignMatrix(chr, m, n, d)
  }
  y <- X %*%(mu +  tau *rnorm(dim(X)[2])) + sigma * rnorm(dim(X)[1])
  tvalues <-  gettavlues(y,
                         X,
                         model = model,
                         nc  = nc)
  tvalues_no_est <-  gettavlues(y,
                         X,
                         model = 'simple',
                         theta = c(sigma, tau))
  Tvalues <- matrix(tvalues, nrow=m)
  X_t=Tvalues[1:(m-1),]
  Y_t=Tvalues[2:m,]
  tvalues.vec[j, 1] <- max(abs(tvalues))
  tvalues.vec[j, 2] <- max(abs(tvalues_no_est))
  a <- sum(X_t*Y_t)/sum(X_t*X_t)
  cor.t <- calculateCorrelations(Tvalues, markers= markers)
  x <- 1:max(sum(cor.t > 0.3), 2)
  cor.t <- cor.t[x]
  (beta <- -coef(lm(log(cor.t) ~ x - 1)))
  beta_est <- -log(a)
  beta.vec[j,] <- c(beta, beta_est)
  #res <- Y_t -a*X_t
  Pr <- function(z) {
    1 - exp(-2*C*(1 - pnorm(z)) - 2 * beta * L * z * dnorm(z) * v(z* sqrt(2 * beta * d))) -
      alpha
  }
  crit <- uniroot(Pr, c(1, 10))$root
  Pr <- function(z) {
    1 - exp(-2*C*(1 - pnorm(z)) - 2 * beta * L * z * dnorm(z) * v(z* sqrt(2 * beta_est * d))) -
      alpha
  }
  crit2 <- uniroot(Pr, c(1, 10))$root
  crit.vec[j, 1] <- crit
  crit.vec[j, 2] <- crit2
  cat('t_vals = ',round(tvalues.vec[j,],1),'\n')
  cat('crit   = ',round(crit.vec[j,],1),'\n')
  select.vec[j,] <- tvalues.vec[j,]>crit.vec[j,]
  cat('select   = ',round(select.vec[j,],1),'\n')
}
