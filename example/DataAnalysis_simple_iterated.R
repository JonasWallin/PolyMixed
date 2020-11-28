##
# Analysing the simulated data using the simpler model:
#'   Y = X_k\beta + X v + \epsilon
#'   epsilon \sim N(0, \sigma_Y^2 I)
#'   v_j     \sim N(\mu_0, \tau^2)
# D: 2018-11-29
##
rm(list=ls())
graphics.off()
library(latentMixed)
set.seed(3)
signal_nr <- 1
betas <- NULL #c(50, 250, 750)
nc <- 10 #only used for plotting
###
#
# loading the signal
# and setting up data
###
data <- simulation_GeneData(type  = signal_nr,
                            k     = nc,
                            n     = 500,
                            L     = 100,
                            betas = NULL)
X <- data$X
Y <- data$Y
betas_full  <- iterativeForwardSelect(Y,
                                      X,
                                      model='constant',
                                      nc = nc,
                                      kmax = 20,
                                      p.val = 0.998)
betas_simple <- iterativeForwardSelect(Y,
                                       X,
                                       model='simple',
                                       kmax =20,
                                       p.val = 0.998)
#betas_full_L  <- iterativeForwardSelect(Y, X, model='constant', nc = 10, kmax = 10,likelihood_method=T)
#betas_simple_L <- iterativeForwardSelect(Y, X, model='simple', kmax =10, likelihood_method=T)
cat('true     = ',which(data$beta!=0),'\n')
cat('simple   = ',which(betas_simple!=0),'\n')
#cat('simple_L = ',which(betas_simple_L!=0),'\n')
cat('full     = ',which(betas_full!=0),'\n')
#cat('full_L   = ',which(betas_full_L!=0),'\n')

