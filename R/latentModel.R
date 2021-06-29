##
# R code for fitting simplest latent polienic model
#
# D: 2018-11-23
##
library(svMisc)
library(RcppEigen)
library(numDeriv)
#' setup for computing the loglikelihood for the model
#' Assumes equally spaced covariates
#'   Y       \sim N(X_k\beta + X \gamma, \sigma_Y^2 I_{n \times n})
#'   \gamma  \sim N(u, \tau^2 I_{p \times p})
#'   u       = [u_1,\ldots,u_{n_c}]
#'   u_i     \sim N(\mu_0, \sigma^2_u 1 %*% 1 ^T)
#'  @param Y     (n x 1) the data
#'  @param X     (n x p) the gene matrixs
#'  @param Xk    (n x k) the covariates (for \beta) if no covariates set NULL
#'  @param nc    (int - matrix)   number of chromosones, if null simpler model
#'
#'  @param SVDX  (list) output of SVD(X)
#'
loglikSetup <- function(Y, X, Xk = NULL, nc=NULL, SVDX = NULL){

  likObj <- list()
  likObj$Y   <- Y
  likObj$X   <- X
  likObj$Xk  <- Xk

  if(is.null(SVDX)){
    likObj$SVDX   <- svd(X, nu= nrow(X))
    d_ <- nrow(X) - length(likObj$SVDX$d)
    if(d_ > 0)
      likObj$SVDX$d <- c(likObj$SVDX$d, rep(0, d_))
  }else{
    likObj$SVDX   <- SVDX
    d_ <- nrow(X) - length(likObj$SVDX$d)
    if(d_ > 0)
      likObj$SVDX$d <- c(likObj$SVDX$d, rep(0, d_))
  }
  p <- dim(X)[2]
  if(is.null(nc)==F){
    SigmaHalf <- kronecker(diag(nc), rep(1, p/nc))
    H <- (likObj$SVDX$d*t(likObj$SVDX$v))%*%SigmaHalf
    likObj$SigmaHalf <- SigmaHalf
    likObj$Sigma_mu <-(SigmaHalf)%*%t(SigmaHalf)
    likObj$H <- H
    likObj$HHt <- H%*%t(H)
  }
  if(is.null(Xk) ==F){
    likObj$Xu   <- t(likObj$SVDX$u)%*%as.matrix(likObj$Xk)
  }else{
    likObj$Xu   <- NULL
  }
  likObj$Yu   <- t(likObj$SVDX$u)%*%Y
  likObj$UX   <- t(likObj$SVDX$u)%*%likObj$X
  return(likObj)
}

BICmixed.R <- function(theta, likObj, tau_fixed=F, Xu = NULL, find=NULL, C = NULL){
  #'
  #' computing the BIC for the model in loglikSetup function
  #' Not implimented for \sigma_U
  #'  @param theta (3 x 1) log(\sigma_Y), log(\tau), log(\sigma_U)
  #'  @param likobj (list) contaning:
  #'                 Yu     (n x 1) the data (modifed by U in the SVD(X))
  #'                 Xu     (n x k + 1) the covariates, for \beta (modifed by U in the SVD(X))
  #'                 H      (n x n_c) covariance matrix
  #'                 HHt    (n x n)
  #'                 SVDX   (list) singular value decomposition of X
  #'  @param tau_fixed - (bool) assuming that \tau,\sigma are fixed
  #'  @param Xu        - (k x 1) covariates supplied (optional)
  #'  @param find      - (m x 1) number of ocv
  #'  @param C         - (double) if non NULL it is a MBIC
  #'
  #'

  if(is.null(Xu))
    Xu  = likObj$Xu
  if(is.null(find)==F){
    if(is.null(Xu)){
      Xe = likObj$UX[,find]
    }else{
      Xu <- cbind(Xu, likObj$UX[,find])
    }
  }


  lik <- -2*loglikSimpleR(theta, likObj, Xu = Xu, BIC = T, restricted=F)
  if(is.null(C)==F){
    #mBIC2
    p = dim(likObj$X)[2]
    k = dim(Xu)[2]
    if(is.null(k))
      k = 0
    lik <- lik + 2 * k * log(p/C) - 2 * lgamma(k+1)

  }
  if(tau_fixed == T)
    return(lik)


  lambda <- function(theta){
    return(loglikSimpleR(theta, likObj, Xu = Xu, BIC = T, restricted=F))
  }
  J <- -numDeriv::hessian(lambda, theta)
  BIC = lik + determinant(J)$modulus[1]
  return(BIC)
}


loglikSimpleR <- function(theta, likObj, restricted = T, Xu = NULL, BIC = F){
  #'
  #' computing the loglikelihood for the model in loglikSetup function
  #'
  #'  @param theta (3 x 1) log(\sigma_Y), log(\tau), log(\sigma_U)
  #'  @param likobj (list) contaning:
  #'                 Yu     (n x 1) the data (modifed by U in the SVD(X))
  #'                 Xu     (n x k + 1) the covariates, for \beta (modifed by U in the SVD(X))
  #'                 H      (n x n_c) covariance matrix
  #'                 HHt    (n x n)
  #'                 SVDX   (list) singular value decomposition of X
  #'  @param restricted (bool) use restricted maximum likelihood
  #'  @param Xu   - (k x 1) covariates supplied (optional)
  #'  @param BIC  - (bool) caculate BIC adjusted
  #'
  #'
  if(BIC)
    restricted = F
  Sigma_U <- getSigmaUY(exp(theta), likObj)
  if(is.vector(Sigma_U)){
    lik = - 0.5 * sum(log(Sigma_U))
    SinvY <- likObj$Yu/Sigma_U
  }else{
    R <- chol(Sigma_U)
    lik <- - sum(log(diag(R)))
    SinvY <- backsolve(R,forwardsolve(R, likObj$Yu, upper.tri=T,transpose=T))
  }
  lik = lik - 0.5 * ( t(likObj$Yu)%*%SinvY)
  #if no covarietes

  if(is.null(Xu))
    Xu  = likObj$Xu
  if(is.null(Xu))
    return(lik)

  XtSinvY <- t(Xu)%*%SinvY
  if(is.vector(Sigma_U)){
    SinvX <- Xu / Sigma_U
  }else{
    SinvX <- backsolve(R,forwardsolve(R, Xu, upper.tri=T,transpose=T))
  }
  XtSinvX <- t( Xu)%*%SinvX
  XtSinvX <- 0.5 * (XtSinvX + t(XtSinvX))
  R2 <- chol(XtSinvX)

  if(restricted)
    lik = lik - sum(log(diag(R2)))

  if(BIC)
    lik = lik - sum(log(diag(R2)))
  XSiXi_XSiY <- backsolve(R2, forwardsolve(R2, XtSinvY, upper.tri=T,transpose=T))
  lik = lik + 0.5 * ( t(XtSinvY) %*% XSiXi_XSiY)
  return(lik)
}


#'
#' whiting selected vectors
#' @param GeneMixObj  -
#' @param theta       - (vec) parameters for variance model
#' @param Xout        - (bool) whitning entire gene matrix
#' @param Yout        - (bool) whitning output
#' @param XmOut       - (bool) whitning covariates
#' @param Xextra      - (n x p_0) extra covariates to whiten (should be in U coordinates (i.e X- SVD basis))
whiteData <- function(GeneMixObj, theta = NULL, XOut=T, YOut=T, XmOut = T, Xextra=NULL){

  if(is.null(theta ))
    theta <- GeneMixObj$theta0

  res <- list(X = NULL, Y = NULL, Xm = NULL)
  if (!GeneMixObj$tauOff) {
    Sigma_U <- getSigmaUY(exp(theta), GeneMixObj$LikObj)
    if(is.vector(Sigma_U)){
      if(XOut)
        res$X  <- GeneMixObj$LikObj$UX/sqrt(Sigma_U)
      if(YOut)
        res$Y  <- GeneMixObj$LikObj$Yu/sqrt(Sigma_U)
      if(XmOut)
        res$Xm <- GeneMixObj$LikObj$Xu/sqrt(Sigma_U)
      if(!is.null(Xextra))
        res$Xextra <-Xextra /sqrt(Sigma_U)
    }else{
      #needs speed up!
      #and wrong need speed up using inverses!!
      stop("error message")
      R  <- chol(Sigma_U)
      y  <- backsolve(R,forwardsolve(R, GeneMixObj$LikObj$Yu, upper.tri=T,transpose=T))
      X  <- backsolve(R,forwardsolve(R, GeneMixObj$LikObj$UX, upper.tri=T,transpose=T))
      Xm <- backsolve(R,forwardsolve(R,GeneMixObj$LikObj$Xu,  upper.tri=T,transpose=T))
    }
  }else{
    if(YOut)
      res$Y  <- GeneMixObj$LikObj$Yu
    if(XOut)
      res$X  <- GeneMixObj$LikObj$UX
    if(XmOut)
      res$Xm <- GeneMixObj$LikObj$Xu
    if(!is.null(Xextra))
      res$Xextra <- Xextra
  }
  return(res)
}
#'
#' gets the ML estimate of the covariates given coeffients
#'  @param theta (3 x 1) log(\sigma_Y), log(\tau), log(\sigma_U)
#'  @param likobj (list) contaning:
#'                 Yu     (n x 1) the data (modifed by U in the SVD(X))
#'                 Xu     (n x k + 1) the covariates, for \beta (modifed by U in the SVD(X))
#'                 H      (n x n_c) covariance matrix
#'                 HHt    (n x n)
#'                 SVDX   (list) singular value decomposition of X
#'  @param restricted (bool) use restricted maximum likelihood
#'  @param Xu   - (k x 1) covariates supplied (optional)
getBeta <- function(theta, likObj, Xu = NULL){


  if(is.null(Xu))
    Xu  = likObj$Xu
  if(is.null(Xu))
    return(matrix(nrow=0,ncol=0))

  Sigma_U <- getSigmaUY(exp(theta), likObj)
  if(is.vector(Sigma_U)){
    SinvY <- likObj$Yu/Sigma_U
  }else{
    R <- chol(Sigma_U)
    SinvY <- backsolve(R,forwardsolve(R, likObj$Yu, upper.tri=T,transpose=T))
  }
  XtSinvY <- t(Xu)%*%SinvY
  if(is.vector(Sigma_U)){
    SinvX <- Xu / Sigma_U
  }else{
    SinvX <- backsolve(R,forwardsolve(R, Xu, upper.tri=T,transpose=T))
  }
  XtSinvX <- t( Xu)%*%SinvX
  R2 <- chol(XtSinvX)
  beta <- backsolve(R2, forwardsolve(R2, XtSinvY, upper.tri=T,transpose=T))
  return(beta)
}

getSigmaBeta <- function(theta, likObj){
  #'
  #'  Returns the covariance matrix of Beta
  #'
  #'  @param theta (3 x 1) \sigma_Y, \tau, \sigma_U
  #'  @param likobj (list) contaning:
  #'                 Yu     (n x 1) the data (modifed by U in the SVD(X))
  #'                 Xu     (n x k + 1) the covariates, for \beta, \mu_0 (modifed by U in the SVD(X))
  #'                 H      (n x n_c) covariance matrix
  #'                 HHt    (n x n)
  #'                 SVDX   (list) singular value decomposition of X

  Sigma_U <- getSigmaUY(exp(theta), likObj)
  if(length(theta) == 3){
    CovBeta <- solve( t(likObj$Xu) %*% solve(Sigma_U, likObj$Xu))
  }else{
    CovBeta <- solve( t(likObj$Xu) %*% (likObj$Xu / Sigma_U))
  }
  return(CovBeta)
}
getSigmaUY  <- function(theta, likObj){
  #'
  #'  Returns the covariance matrix of t(U)%*%Y
  #'  if model is simple return the diagonal vector
  #'  @param theta (3 x 1) \sigma_Y, \tau, \sigma_U
  #'  @param likobj (list) contaning:
  #'                 Yu     (n x 1) the data (modifed by U in the SVD(X))
  #'                 Xu     (n x k + 1) the covariates, for \beta, \mu_0 (modifed by U in the SVD(X))
  #'                 H      (n x n_c) covariance matrix
  #'                 HHt    (n x n)
  #'                 SVDX   (list) singular value decomposition of X
  sigma_Y = theta[1]
  if(length(theta) == 1)
    D <- sigma_Y^2* diag(nrow=length(likObj$SVDX$d))
  if(length(theta)> 1){
    tau     = theta[2]
    D       <-  tau^2*likObj$SVDX$d^2 + sigma_Y^2
  }
  if(length(theta)>2){
    sigma_U = theta[3]
    Sigma_U <- sigma_U^2 * likObj$HHt
    diag(Sigma_U) <- diag(Sigma_U) + D
  }else{
    Sigma_U <- D
  }

  return(Sigma_U)
}

getMU<- function(theta, likObj){
  #'
  #'  Get the \mu,v estimate (i.e. kriging)
  #'
  #'
  #'  @param theta (3 x 1) log(\sigma_Y), log(\tau), log(\sigma_U)
  #'  @param likobj (list) contaning:
  #'                 Yu     (n x 1) the data (modifed by U in the SVD(X))
  #'                 Xu     (n x k + 1) the covariates, for \beta, \mu_0 (modifed by U in the SVD(X))
  #'                 H      (n x n_c) covariance matrix
  #'                 HHt    (n x n)
  #'                 SVDX   (list) singular value decomposition of X
  #'
  #'

  beta_mu <- getBetaMu0(theta,likObj)
  beta_tilde <- c(beta_mu$beta, beta_mu$mu)


  Sigma_U <- getSigmaUY(exp(theta), likObj)


  tau = exp(theta[2])
  Yu <- likObj$Yu - likObj$Xu%*%beta_tilde
  if(length(theta)>2){
    R        <- chol(Sigma_U)
    SinvYu   <- backsolve(R,forwardsolve(R, Yu, upper.tri=T,transpose=T))
    sigma_U = exp(theta[3])
    Sigma_mu <- sigma_U^2*likObj$Sigma_mu

    mu_hat   <- Sigma_mu%*%( t(likObj$UX)%*%SinvYu)
    Sigma_v  <- Sigma_mu
    diag(Sigma_v) <- diag(Sigma_v) + tau^2
    hats<-list(mu_hat = mu_hat)
  }else{
    SinvYu <- Yu/Sigma_U
    Sigma_v <- tau^2 * diag(dim(X)[2])
    hats <- list()
  }
  v_hat    <- Sigma_v%*%( t(likObj$UX)%*%SinvYu)
  hats$v_hat <- v_hat
  return(hats)
}

getCov <- function(theta, likObj){
  #'
  #'  Get the \mu,v covariances (not joint)
  #'
  #'
  #'  @param theta (3 x 1) log(\sigma_Y), log(\tau), log(\sigma_U)
  #'  @param likobj (list) contaning:
  #'                 Yu     (n x 1) the data (modifed by U in the SVD(X))
  #'                 Xu     (n x k + 1) the covariates, for \beta, \mu_0 (modifed by U in the SVD(X))
  #'                 H      (n x n_c) covariance matrix
  #'                 HHt    (n x n)
  #'                 SVDX   (list) singular value decomposition of X
  #'
  #'

  Sigma_U <- getSigmaUY(exp(theta), likObj)

  tau  = exp(theta[2])

  iSigma_U <- solve(Sigma_U)
  XtiSX <- t(likObj$UX)%*% (iSigma_U %*%likObj$UX)
  if(length(theta)>2){
    sigma_U <- exp(theta[3])
    Sigma_h <- sigma_U*likObj$SigmaHalf
    Sigma_muY <- XtiSX
    diag(Sigma_muY) <- diag(Sigma_muY) - 1
    Sigma_muY <- Sigma_h%*%(Sigma_muY %*% t(Sigma_h))
    Sigma_v  <- Sigma_mu
    diag(Sigma_v) <- diag(Sigma_v) + tau^2

    hats<-list(Sigma_muY = Sigma_muY)
  }else{

    Sigma_v <- tau^2 * diag(dim(X)[2])
    hats <- list()
  }
  Sigma_vY <- Sigma_v - Sigma_v%*%XtiSX%*%Sigma_v
  hats$Sigma_vY <- Sigma_vY
  return(hats)
}
getCov_old<- function(theta, likObj){
  #'
  #'  Get the \mu,v covariances (not joint)
  #'
  #'
  #'  @param theta (3 x 1) log(\sigma_Y), log(\tau), log(\sigma_U)
  #'  @param likobj (list) contaning:
  #'                 Yu     (n x 1) the data (modifed by U in the SVD(X))
  #'                 Xu     (n x k + 1) the covariates, for \beta, \mu_0 (modifed by U in the SVD(X))
  #'                 H      (n x n_c) covariance matrix
  #'                 HHt    (n x n)
  #'                 SVDX   (list) singular value decomposition of X
  #'
  #'

  Sigma_U <- getSigmaUY(exp(theta), likObj)

  tau  = exp(theta[2])


  if(length(theta)>2){
    iSigma_U <- solve(Sigma_U)
    XtiSX <- t(likObj$UX)%*% (iSigma_U %*%likObj$UX)
    sigma_U <- exp(theta[3])
    Sigma_mu <- sigma_U^2*likObj$Sigma_mu
    Sigma_muY <- Sigma_mu - Sigma_mu%*%XtiSX%*%Sigma_mu
    Sigma_v  <- Sigma_mu
    diag(Sigma_v) <- diag(Sigma_v) + tau^2

    hats<-list(Sigma_muY = Sigma_muY)
  }else{
    XtiSX <- t(likObj$UX)%*% (likObj$UX / Sigma_U)
    Sigma_v <- tau^2 * diag(dim(X)[2])
    hats <- list()
  }
  Sigma_vY <- Sigma_v - Sigma_v%*%XtiSX%*%Sigma_v
  hats$Sigma_vY <- Sigma_vY
  return(hats)
}


simulatetvalues<- function(Y,
                           X,
                           model = c('simple','constant'),
                           betas = NULL,
                           nc   = NULL,
                           theta = NULL,
                           SVDX = SVDX){

  model <- match.arg(model)
  p <- dim(X)[2]

  if(model=='simple'){
    theta0 <- c(log(sd(Y)), log(sd(Y))-2)
  }else{

    theta0 <- c(log(sd(Y)), log(sd(Y))-2, log(sd(Y))-2)
  }


  if(is.null(betas))
    betas <- rep(0, p)

  # the iteration procedure
  likObj <- loglikSetup(Y, Xk = X[,betas != 0], X, nc = nc, SVDX = SVDX)
  if(is.null(theta)){

    lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
    res <- optim(theta0, lambda)
    theta <- exp(res$par)
  }

  Sigma_U <- getSigmaUY(theta, likObj)
  if(is.vector(Sigma_U)==F){
    iSY <- solve(Sigma_U,likObj$Yu)
  }else{
    iSY <- likObj$Yu / Sigma_U
  }
  n_p <- dim(likObj$Xu)[2]
  if(is.vector(Sigma_U)){
    XtiSX <- t(likObj$Xu)%*%(likObj$Xu/ Sigma_U)
  }else{
    XtiSX <- t(likObj$Xu)%*%solve(Sigma_U, likObj$Xu)
  }
  CovBeta <- solve(XtiSX)
  beta_est <- CovBeta%*%( t(likObj$Xu)%*%iSY)

  if(model=='simple'){
    y_sim <- (X%*% rnorm(n=dim(X)[2],mean = beta_est[n_p], sd = theta[2])) + rnorm(n=dim(X)[1],sd=theta[1])
    likObj <- loglikSetup(y_sim, Xk = NULL, X, nc = nc, SVDX = SVDX)
    iSY <- likObj$Yu / Sigma_U
    tvalues = getTValuesDiagS(1:p,
                              likObj$Xu,
                              likObj$UX,
                              Sigma_U,
                              iSY)
  }
  return(tvalues)
}

gettvaluesCond  <- function(index,
                            Y,
                            X,
                            model = c('simple','constant'),
                            betas = NULL,
                            nc    = NULL,
                            theta = NULL){
  #'
  #' Caculates the t-values for X_i in :
  #'   Y       \sim N(X_i \beta_i + X \gamma, \sigma_Y^2 I_{n \times n})
  #'   \gamma  \sim N(u, \tau^2 I_{p \times p})
  #'   u       = [u_1,\ldots,u_{n_c}]
  #'   u_i     \sim N(\mu_0, \sigma^2_u 1 %*% 1 ^T)
  #'
  #' @param index - (p_1 x 1) vector of position where to calc t-values
  #' @param Y     - (n x 1) observations
  #' @param X     - (n x p) covariates
  #' @param model - (string)
  #'                 simple   - simpliefed model (sigma^2_u = 0)
  #'                 constant - the full model
  #' @param nc    - (int) number of chromosones
  #' @param theta - (k x 1) parameters of the model if NULL estimated from data
  #'
  #' @return t - (p x 1) vector of the t-values
  #'

  model <- match.arg(model)
  p <- dim(X)[2]

  if(model=='simple'){
    theta0 <- c(log(sd(Y)), log(sd(Y))-2)
  }else{

    theta0 <- c(log(sd(Y)), log(sd(Y))-2, log(sd(Y))-2)
  }


  if(is.null(betas))
    betas <- rep(0, p)

  # the iteration procedure
  likObj <- loglikSetup(Y, Xk = X[,betas != 0], X, nc = nc)
  if(is.null(theta)){

    lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
    res <- optim(theta0, lambda)
    theta <- exp(res$par)
  }

  Sigma_U <- getSigmaUY(theta, likObj)
  if(is.vector(Sigma_U)==F){
    iSY <- solve(Sigma_U,likObj$Yu)
  }else{
    iSY <- likObj$Yu / Sigma_U
  }
  n_p <- dim(likObj$Xu)[2]
  if(is.vector(Sigma_U)){
    XtiSX <- t(likObj$Xu)%*%(likObj$Xu/ Sigma_U)
  }else{
    XtiSX <- t(likObj$Xu)%*%solve(Sigma_U, likObj$Xu)
  }
  CovBeta <- solve(XtiSX)
  beta_est <- CovBeta%*%( t(likObj$Xu)%*%iSY)
  if(is.vector(Sigma_U)==F){
    iSY <-iSY - solve(Sigma_U, likObj$Xu[,-n_p,drop = F]%*%beta_est[-n_p,drop = F])
  }else{
    iSY <- iSY - (likObj$Xu[,-n_p,drop = F]%*%beta_est[-n_p,drop = F]) / Sigma_U
  }
  uX <- t(likObj$SVDX$u)%*%X

   tvalues <- getTValues_R(index,
                            likObj$Xu,
                            uX,
                            Sigma_U,
                            iSY)

  return(tvalues[index])
}
gettvaluesCond2  <- function(index,
                             Y,
                             X,
                             model = c('simple','constant'),
                             betas = NULL,
                             nc   = NULL,
                             theta = NULL,
                             X_rem = NULL,
                             SVDX  = NULL){
  #'
  #' Caculates the t-values for X_i in :
  #'   Y       \sim N(X_i \beta_i + X \gamma, \sigma_Y^2 I_{n \times n})
  #'   \gamma  \sim N(u, \tau^2 I_{p \times p})
  #'   u       = [u_1,\ldots,u_{n_c}]
  #'   u_i     \sim N(\mu_0, \sigma^2_u 1 %*% 1 ^T)
  #'
  #' @param index - (p_1 x 1) vector of position where to calc t-values
  #' @param Y     - (n x 1) observations
  #' @param X     - (n x p) covariates
  #' @param model - (string)
  #'                 simple   - simpliefed model (sigma^2_u = 0)
  #'                 constant - the full model
  #' @param nc    - (int) number of chromosones
  #' @param theta - (k x 1) parameters of the model if NULL estimated from data
  #' @param XSVD  - (list) singular value decomposition of X
  #'
  #' @return t - (p x 1) vector of the t-values
  #'

  model <- match.arg(model)
  p <- dim(X)[2]

  if(model=='simple'){
    theta0 <- c(log(sd(Y)), log(sd(Y))-2)
  }else{


    theta0 <- c(log(sd(Y)), log(sd(Y))-2, log(sd(Y))-2)
  }


  if(is.null(betas))
    betas <- rep(0, p)

  # the iteration procedure
  likObj <- loglikSetup(Y, Xk = X[,betas != 0], X, nc = nc, SVDX = SVDX)
  if(is.null(theta)){
    lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
    res <- optim(theta0, lambda)
    theta <- exp(res$par)
  }

  Sigma_U <- getSigmaUY(theta, likObj)
  if(is.vector(Sigma_U)==F){
    iSY <- solve(Sigma_U,likObj$Yu)
  }else{
    iSY <- likObj$Yu / Sigma_U
  }

  if(is.vector(Sigma_U)){
    tvalues = getTValuesDiagS(setdiff(index, X_rem),
                               likObj$Xu,
                               likObj$UX,
                               Sigma_U,
                               iSY)
  }else{
    tvalues = getTValues(setdiff(index, X_rem),
                          likObj$Xu,
                          likObj$UX,
                          Sigma_U,
                          iSY)
   }
    return(tvalues)
}
gettvaluesCond3 <- function(){

}
gettavlues  <- function(Y,
                        X,
                        model = c('simple','constant'),
                        betas = NULL,
                        nc   = NULL,
                        theta = NULL,
                        SVDX = NULL){
  #'
  #' Caculates the t-values for X_i in :
  #'   Y       \sim N(X_i \beta_i + X \gamma, \sigma_Y^2 I_{n \times n})
  #'   \gamma  \sim N(u, \tau^2 I_{p \times p})
  #'   u       = [u_1,\ldots,u_{n_c}]
  #'   u_i     \sim N(\mu_0, \sigma^2_u 1 %*% 1 ^T)
  #'
  #' @param Y     - (n x 1) observations
  #' @param X     - (n x p) covariates
  #' @param model - (string)
  #'                 simple   - simpliefed model (sigma^2_u = 0)
  #'                 constant - the full model
  #' @param nc    - (int) number of chromosones
  #' @param theta - (k x 1) parameters of the model if NULL estimated from data
  #'
  #' @return t - (p x 1) vector of the t-values
  #'

  model <- match.arg(model)
  p <- dim(X)[2]

  if(model=='simple'){
    theta0 <- c(log(sd(Y)), log(sd(Y))-2)
  }else{

    theta0 <- c(log(sd(Y)), log(sd(Y))-2, log(sd(Y))-2)
  }


  if(is.null(betas))
    betas <- rep(0, p)

  # the iteration procedure
  if(is.null(SVDX))
    SVDX = svd(X, nu= nrow(X))
  likObj <- loglikSetup(Y, Xk = X[,betas != 0], X, nc = nc, SVDX = SVDX)
  if(is.null(theta)){

    lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
    res <- optim(theta0, lambda)
    theta <- exp(res$par)
  }

  Sigma_U <- getSigmaUY(theta, likObj)
  likObj <- loglikSetup(Y, Xk = NULL, X, nc = nc, SVDX = SVDX)
  if(is.vector(Sigma_U)==F){
    iSY <- solve(Sigma_U,likObj$Yu)
  }else{
    iSY <- likObj$Yu / Sigma_U
  }

  uX <- t(likObj$SVDX$u)%*%X

  if(is.vector(Sigma_U)){
    tvalues = getTValuesDiagS(1:p,
                         likObj$Xu,
                         uX,
                         Sigma_U,
                         iSY)
  }else{
    tvalues = getTValues(1:p,
                         likObj$Xu,
                         uX,
                         Sigma_U,
                         iSY)
  }
  return(tvalues)
}



backStepIter  <- function(Y,
                          X,
                          model = c('simple','constant'),
                          betas,
                          critval = qnorm(1-(1-0.95)/2),
                          nc   = NULL,
                          return_full = FALSE,
                          SVDX = NULL
                          ){
  #'
  #' backwards removes \beta_i in model
  #'   Y       \sim N(X \beta + X \gamma, \sigma_Y^2 I_{n \times n})
  #'   \gamma  \sim N(u, \tau^2 I_{p \times p})
  #'   u       = [u_1,\ldots,u_{n_c}]
  #'   u_i     \sim N(\mu_0, \sigma^2_u 1 %*% 1 ^T)
  #'
  #' @param Y        - (n x 1) observations
  #' @param X        - (n x p) covariates
  #' @param model    - (string)
  #'                    simple   - simpliefed model (sigma^2_u = 0)
  #'                    constant - the full model
  #' @param nc       - (int) number of chromosones
  #' @param critval     - (double) value to minimum accepted level
  #' @param return_full - (bool) return either betas or betas and params
  #' @return betas      - (p x 1) vector of the beta values
  #'

  model <- match.arg(model)
  p <- dim(X)[2]

  if(model=='simple'){
    theta0 <- c(log(sd(Y)), log(sd(Y))-2)
  }else{

    theta0 <- c(log(sd(Y)), log(sd(Y))-2, log(sd(Y))-2)
  }
  if(is.null(SVDX))
    SVDX = svd(X, nu= nrow(X))
  likObj <- loglikSetup(Y, Xk = X[,betas != 0], X, nc = nc, SVDX = SVDX)
  lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
  res <- optim(theta0, lambda)
  while(sum(betas!=0) > 0){
    # the iteration procedure

    Sigma_U <- getSigmaUY(exp(res$par), likObj)
    if(is.vector(Sigma_U)){
      iSY <- likObj$Yu / Sigma_U
    }else{
      iSY <- solve(Sigma_U, likObj$Yu)
    }

    uX <- t(likObj$SVDX$u)%*%X
    tvalues <- rep(0, p)
    n_p <- dim(likObj$Xu)[2]
    if(is.vector(Sigma_U)){
      XtiSX <- t(likObj$Xu)%*%(likObj$Xu/ Sigma_U)
    }else{
      XtiSX <- t(likObj$Xu)%*%solve(Sigma_U, likObj$Xu)
    }
    CovBeta <- solve(XtiSX)
    beta_est <- CovBeta%*%( t(likObj$Xu)%*%iSY)
    tvalues <- abs(beta_est[1:(n_p-1)]/sqrt(diag(CovBeta)[1:(n_p-1)]))
    if(min(tvalues) > critval){
      betas[betas != 0] = beta_est[1:(n_p-1)]
      if(return_full){
        return(list(betas = betas,theta = exp(res$par)))
      }
      return(betas)
    }
    beta_ind <- which(betas!=0)
    betas[beta_ind[which(tvalues == min(tvalues))]] =0
    likObj <- loglikSetup(Y, Xk = X[,betas != 0], X, nc = nc, SVDX = XSVD)
    lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
    res <- optim(theta0, lambda)
  }
  if(return_full){
    return(list(betas = betas,theta = exp(res$par)))
  }
  return(betas)
}

backStepIter2  <- function(Y,
                          X,
                          model = c('simple','constant'),
                          betas,
                          critval = qnorm(1-(1-0.95)/2),
                          nc   = NULL,
                          SVDX = NULL){
  #'
  #' backwards removes \beta_i in model
  #'   Y       \sim N(X \beta + X \gamma, \sigma_Y^2 I_{n \times n})
  #'   \gamma  \sim N(u, \tau^2 I_{p \times p})
  #'   u       = [u_1,\ldots,u_{n_c}]
  #'   u_i     \sim N(\mu_0, \sigma^2_u 1 %*% 1 ^T)
  #'
  #' @param Y        - (n x 1) observations
  #' @param X        - (n x p) covariates
  #' @param model    - (string)
  #'                    simple   - simpliefed model (sigma^2_u = 0)
  #'                    constant - the full model
  #' @param nc       - (int) number of chromosones
  #' @param critval  - (double) value to minimum accepted level
  #'
  #' @return betas - (p x 1) vector of the beta values
  #'

  model <- match.arg(model)
  p <- dim(X)[2]

  if(model=='simple'){
    theta0 <- c(log(sd(Y)), log(sd(Y))-2)
  }else{
    theta0 <- c(log(sd(Y)), log(sd(Y))-2, log(sd(Y))-2)
  }

  if(is.null(SVDX))
    SVDX = svd(X, nu= nrow(X))
  while(sum(betas!=0) > 0){
    # the iteration procedure
    beta_ind <- which(betas!=0)
    tvalues <- rep(0, length(beta_ind))
    for(i in 1:length(beta_ind)){
      betas_i <- betas
      betas_i[beta_ind[i]] <- 0
      likObj <- loglikSetup(Y, Xk = X[,betas_i != 0], X, nc = nc, SVDX = SVDX)
      lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
      res <- optim(theta0, lambda)
      theta0 <- res$par
      Sigma_U <- getSigmaUY(exp(res$par), likObj)
      if(is.vector(Sigma_U)==F){
        iSY <- solve(Sigma_U,likObj$Yu)
      }else{
        iSY <- likObj$Yu  / Sigma_U
      }
      n_p <- dim(likObj$Xu)[2]
      if(is.vector(Sigma_U)){
        XtiSX <- t(likObj$Xu)%*%(likObj$Xu/ Sigma_U)
      }else{
        XtiSX <- t(likObj$Xu)%*%solve(Sigma_U, likObj$Xu)
      }
      CovBeta <- solve(XtiSX)
      beta_est <- CovBeta%*%( t(likObj$Xu)%*%iSY)
      if(is.vector(Sigma_U)==F){
        iSY <-iSY - solve(Sigma_U, likObj$Xu[,-n_p,drop = F]%*%beta_est[-n_p,drop = F])
      }else{
        iSY <- iSY - (likObj$Xu[,-n_p,drop = F]%*%beta_est[-n_p,drop = F]) / Sigma_U
      }
      uX <- t(likObj$SVDX$u)%*%X
      t_val_i <- getTValues_R(beta_ind[i],
                              likObj$Xu,
                              uX,
                              Sigma_U,
                              iSY)
      tvalues[i] = t_val_i[beta_ind[i]]
    }
    if(min(abs(tvalues)) > critval){
      likObj <- loglikSetup(Y, Xk = X[,betas != 0], X, nc = nc, SVDX = SVDX)
      lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
      res <- optim(theta0, lambda)
     Sigma_U <- getSigmaUY(exp(res$par), likObj)
      if(is.vector(Sigma_U)==F){
        iSY <- solve(Sigma_U,likObj$Yu)
      }else{
        iSY <- likObj$Yu  / Sigma_U
      }
      n_p <- dim(likObj$Xu)[2]
      if(is.vector(Sigma_U)){
        XtiSX <- t(likObj$Xu)%*%(likObj$Xu/ Sigma_U)
      }else{
        XtiSX <- t(likObj$Xu)%*%solve(Sigma_U, likObj$Xu)
      }
      CovBeta <- solve(XtiSX)
      beta_est <- CovBeta%*%( t(likObj$Xu)%*%iSY)
      betas[betas != 0] = beta_est[-n_p]

      return(betas)
    }
    betas[beta_ind[which(abs(tvalues) == min(abs(tvalues)))]] =0
  }
  return(betas)
}
