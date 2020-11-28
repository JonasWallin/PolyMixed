##
# Analysing the simulated data using the simpler model:
#'   Y = X_k\beta + X v + \epsilon
#'   epsilon \sim N(0, \sigma_Y^2 I)
#'   v_j     \sim N(\mu_0, \tau^2)
# D: 2018-11-23
##
rm(list=ls())
graphics.off()
library(latentMixed)
save.fig = F
datafolder <- 'data/data_sets/'
estimate_beta <- TRUE #if false use beta the beta given
signal_nr <- 3
no_signal = FALSE # data with no signal i.e. no beta
nc <- 10 #only used for plotting
###
#
# loading the signal
# and setting up data
###
if(no_signal==F){
  load(paste(datafolder, 'signal_', signal_nr,'.Rdata', sep=''))
}else{
  load(paste(datafolder, 'data/data_sets/no_signal_', signal_nr, sep=''))
}


if(no_signal==FALSE){
  Xk <- X[,beta != 0] #only use the non-zero elements
  if(estimate_beta==FALSE){
    Y = Y- X%*%beta
    Xk <- NULL
  }
}else{
  Xk <- NULL
  estimate_beta=FALSE
}
likObj <- loglikSetup(Y, Xk, X)

theta0 <- c(log(sd(Y)), -1)
# setting up the optim function
lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
res <- optim(theta0, lambda)
res <- exp(res$par)
names(res) <- c('sigma_Y','sigma_tau')
print(round(res,4))

#colleting the estimate of \mu_0, \beta
beta_mu <- getBetaMu0(log(res), likObj)
Sigma_est <- getCov(log(res), likObj)
cat('beta = ',round(beta_mu$beta,2),'\n')
cat('mu = ',round(beta_mu$mu,5),'\n')
if(save.fig){
  pdf(paste('signal_simple_model',signal_nr,'.pdf',sep=''))
}else{
  x11()
}
# colleting posterior distribution of v|Y
mus <- getMU(log(res), likObj)
sd_v <-sqrt(diag(Sigma_est$Sigma_vY))
ylim = c(beta_mu$mu + min(mus$v_hat-2*sd_v),
         beta_mu$mu + max(mus$v_hat+2*sd_v))
plot(c(mus$v_hat),
     ylim = ylim,
     type='n',
     ylab='v|y',
     xlab='')
nl <- dim(X)[2]/nc
abline(h=beta_mu$mu,lwd=2)
for(i in 1:nc){
  index = 1:nl + nl*(i-1)
  lines(index,beta_mu$mu +  mus$v_hat[index],col='blue')
  abline(v=nl*i,lty=2, lw=0.8, col = 'gray')
  lines(index,beta_mu$mu + mus$v_hat[index] + 2*sd_v[index],col='blue',lty=2)
  lines(index,beta_mu$mu + mus$v_hat[index] - 2*sd_v[index],col='blue',lty=2)

}

if(save.fig)
  dev.off()



