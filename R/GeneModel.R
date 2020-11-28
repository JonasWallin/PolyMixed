#mixed model

# fit various mixed model
# caculate t-values

#likObj
#: model-type
#: additional covariates
#: ignore tau
#: find - the Xs we are intressted in

library(numDeriv)
library(RcppEigen)

#'
#' Base class for gene model
#'   Y       \sim N(X_k\beta + X \gamma, \sigma_Y^2 I_{n \times n})
#'   \gamma  \sim N(u, \tau^2 I_{p \times p})
#'   u       = [u_1,\ldots,u_{n_c}]
#'   u_i     \sim N(\mu_0, \sigma^2_u 1 %*% 1 ^T)
#' @param formula     - (formula)   the covariates formula
#' @param data        - (data.frame) for covariates (X_k) these should not include X
#' @param X           - (n x p) the gene covariates
#' @param tauOff      - (bool) fix \tau = 0,
#' @param sigmaMuOff  - (bool) \sigma_u = 0, if sigmaMuOff=FALSE needs to include nc
#' @param meanMuOff   - (bool) \mu_0,= 0,
#' @param nc          - (int) number of chromosones, if null simpler model
#' @param SVDX        - (obj) singular value decomposition object of X
#' @param nPC         - (int) how many PC to use in fixed effect
#' @param testTau     - (bool) should we try to turn off tau?
#'
#' @export
SetupGeneMix<- function(formula,
                        data,
                        X,
                        tauOff     = FALSE,
                        sigmaMuOff = TRUE,
                        meanMuOff  = FALSE,
                        nc         = NULL,
                        SVDX       = NULL,
                        testTau    = FALSE,
                        nPC        = 0){

  mf <- model.frame(formula = formula, data = data)
  Y  <- as.matrix(model.extract(mf, "response"))
  Xk <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))
  covNames <- c()
  if(dim(Xk)[2]>0)
    covNames <- dimnames(Xk)[[2]]
  if(meanMuOff==FALSE){
    Xk <- as.matrix(cbind(Xk, rowSums(X)))
    covNames <- c(covNames,'mu_0')
  }

  if(is.null(SVDX))
    SVDX <- svd(X, nu=dim(X)[1])

  if(dim(SVDX$u)[2] != dim(X)[1])
    stop("the dim(SVD$u)[2] must equal dim(X)[1]")

  if(nPC > 0){
    Xk <- cbind(Xk, SVDX$d[1:nPC]*SVDX$u[,1:nPC] )
    covNames <- c(covNames,paste('PC',1:nPC,sep=""))
  }
  GeneMixObj        <- list(formula = formula,
                            data = data,
                            tauOff = tauOff,
                            sigmaMuOff = sigmaMuOff,
                            meanMuOff  = meanMuOff,
                            covNames   = covNames,
                            nPC        = nPC)
  GeneMixObj$testTau <- testTau
  GeneMixObj$LikObj <- loglikSetup(Y, X, Xk = Xk, nc=nc, SVDX = SVDX)
  class(GeneMixObj) <- "GeneMixObj"

  return(GeneMixObj)
}

#' estimate GeneMixObj
#' @param GeneMixObj - gene mixture object
#' @param theta0     - known intial parameters deafult null
estimateGeneMix <- function(GeneMixObj, theta0 = NULL, find =NULL, restricted=T){

  if(is.null(theta0)){
    if(is.null(GeneMixObj$theta0)){
      sdY <- sd(GeneMixObj$LikObj$Y)
      theta0 <- log(sdY)
      if(GeneMixObj$tauOff ==F)
        theta0 <- c(theta0, log(sdY) -2)
      if(GeneMixObj$sigmaMuOff == F)
        theta0 <- c(theta0, log(sdY) -2)
    }else{
      theta0 <- GeneMixObj$theta0
    }
  }
  Xu = GeneMixObj$LikObj$Xu
  if(!is.null(find))
    Xu <- cbind(Xu, GeneMixObj$LikObj$UX[,find])
  if(length(theta0) > 1){
    lambda <- function(theta){-loglikSimpleR(theta, GeneMixObj$LikObj, Xu = Xu ,restricted = restricted)}
    res <- optim(theta0, lambda)
    GeneMixObj$theta0 <- res$par
    GeneMixObj$RLik <- -res$value
    GeneMixObj$sigma <- exp(res$par[1])
  }else{
    linear <- lm('Yu ~ -1 + .', data = data.frame(Yu = MixGeneObj$LikObj$Yu, Xu))
    GeneMixObj$sigma <- summary(linear)$sigma
    GeneMixObj$beta <- as.vector(linear$coefficients)
  }


  if(GeneMixObj$tauOff ==F){
    GeneMixObj$tau <- exp(res$par[2])
    if(GeneMixObj$sigmaMuOff == F)
      GeneMixObj$sigmaMu <- exp(res$par[3])
  }else{
    GeneMixObj$tau <- 0
    if(GeneMixObj$sigmaMuOff == F)
      GeneMixObj$sigmaMu <- exp(res$par[2])
  }
  if(length(theta0) > 1)
    GeneMixObj$beta <- as.vector(getBeta(GeneMixObj$theta0, GeneMixObj$LikObj, Xu))
  if(length(GeneMixObj$beta)>0){
    if(length(find)==0){
      names(GeneMixObj$beta) <- GeneMixObj$covNames
    }else{
      names(GeneMixObj$beta) <- c(GeneMixObj$covNames,paste('X',find,sep=''))
    }
  }
  return(GeneMixObj)
}

#'
#' mBIC
#' takes a an object and a find list
#' and caculates the mBIC value
#'
#' @param GeneMixObj
#' @param find       - (k x 1) known non zero genes
#'
mBIC2_GeneMixObj <- function(GeneMixObj, find, const=  4){

  GeneMixObj <- estimateGeneMix(GeneMixObj, find = find, restricted = F)
  loglik     <- -loglikSimpleR(GeneMixObj$theta0,
                              GeneMixObj$LikObj,
                              restricted = F)
  k <- length(GeneMixObj$beta) + length(GeneMixObj$theta0)
  n <- length(GeneMixObj$LikObj$Y)
  p <- dim(GeneMixObj$LikObj$UX)[2]
  penalty <- 2 * lgamma(k+1)
  mbic = -2 * loglik + k * log(n) + 2 * k * log(p/const - 1)  - penalty
  return(mbic)
}
#'
#'
#' backward selection of find
#'
backward_GeneMixObj <-function(GeneMixObj, find, crit = mBIC2_GeneMixObj, silent = T){
  crit_v <- crit(GeneMixObj, find)
  stop = 0
  i <- -1
  while(stop==0){
    step=1
    if(silent==F)
      cat('find length = ',length(find),'rem = ',i,' critv = ',crit_v,'\n')
    for(i in 1:length(find)){
      find_i <- find[-i]
      crit_i <- crit(GeneMixObj, find_i)
      if(crit_i < crit_v){
        crit_v <- crit_i
        find <- find_i
        stop = 0
        break
      }
    }
  }
  return(find)
}
#'
#' estimate the gene model
#' @param GeneMixObj
#' @param find       - (k x 1) known non zero genes
#' @param dupl       - (c x 1) duplicated rows of X (found in GeneMixObj$LikObj$X)
#'
mixedModelOptim <- function(GeneMixObj, find = NULL, dupl = NULL, estPar = TRUE, restricted=T) {

  loglik <- NULL

  if (estPar & !GeneMixObj$tauOff) {
    GeneMixObj <- estimateGeneMix(GeneMixObj, find = find, restricted = restricted)
    L1 <- loglikSimpleR(GeneMixObj$theta0,
                        GeneMixObj$LikObj )
    # probably not correct below discuss
    L0 <- loglikSimpleR(c(GeneMixObj$theta0[1],
                          log(0)),
                        GeneMixObj$LikObj )
    if(GeneMixObj$testTau)
      GeneMixObj$tauOff <- L1 - L0 < qchisq(0.95, 1)/2 # LRT
    GeneMixObj$Lik <- loglikSimpleR(GeneMixObj$theta0,
                                    GeneMixObj$LikObj,
                                    restricted = T)
    GeneMixObj$BIC <- loglikSimpleR(GeneMixObj$theta0,
                                    GeneMixObj$LikObj,
                                    BIC = T)
  }

  GeneMixObj$find <- find

  return(GeneMixObj)
}



#'
#' estimate the gene model and tries to find non-zero coeffients
#' and also give t-values
#' @param GeneMixObj
#' @param find       - (k x 1) known non zero genes
#' @param dupl       - (c x 1) duplicated rows of X (found in GeneMixObj$LikObj$X)
#'
#' @export
mixedModel <- function(GeneMixObj, find = NULL, dupl = NULL, estPar = TRUE, restricted=T) {

  GeneMixObj <- mixedModelOptim(GeneMixObj, find, dupl, estPar, restricted)
  WhiteData <- whiteData(GeneMixObj)
  Xm <- WhiteData$Xm
  y  <- WhiteData$Y
  X  <- WhiteData$X
  Xm <- cbind(rep(0, dim(X)[1]), Xm)
  t <- rep(0, dim(X)[2])
  if (!is.null(find)) Xm <- cbind(Xm, X[, find])
  for (j in setdiff(1:dim(X)[2], dupl[, 1])) {
    Xm[, 1] <- X[, j]
    if (j %in% find) {
      fit <- RcppEigen::fastLmPure(unique(Xm, MARGIN = 2), y)
    } else {
      fit <- RcppEigen::fastLmPure(Xm, y)
    }
    t[j] <- fit$coefficients[1] / fit$se[1]
  }
  if(!is.null(dupl)){
    while(sum(abs(t[dupl[,1]]-t[dupl[,2]]))>0)
      t[dupl[,1]] <- t[dupl[,2]]
  }
  GeneMixObj$t <- t
  GeneMixObj$find <- find

  return(GeneMixObj)
}


#'
#'
#' A more general backward model for geneObject
#' @pram restricted - using restricted ML to estimate parameters
#' @param C         - mBIO2 parameter
GeneBackward <- function(GeneMixObj,
                        dupl = NULL,
                        restricted=T,
                        C = 4,
                        silent =F) {

  if(is.null(GeneMixObj$theta0))
    GeneMixObj <- mixedModelOptim(GeneMixObj,
                             find = GeneMixObj$find,
                             dupl = dupl,
                             restricted=restricted)
  BIC0  <- BICmixed.R(GeneMixObj$theta0,
                      GeneMixObj$LikObj,
                      find =  GeneMixObj$find,
                      C= C)

  while(TRUE){
    find_base <- GeneMixObj$find
    step <- FALSE
    if(silent==F)
      cat('*')
    for(i in 1:length(find_base)){ # liberal critical value
      find_i <- setdiff(find_base, find_base[i])
      GeneMixObj_new <- mixedModelOptim(GeneMixObj,
                                   find = find_i,
                                   dupl = dupl,
                                   restricted=restricted)
      GeneMixObj_new$BIC <- BICmixed.R(GeneMixObj_new$theta0,
                                       GeneMixObj_new$LikObj,
                                       find = find_i,
                                       C= C)


      if(GeneMixObj_new$BIC < BIC0){
        step=T
        GeneMixObj_selected <- GeneMixObj_new
        BIC0 <- GeneMixObj_new$BIC
        if(silent==F)
          cat('+')
      }
    }
    if(step==T){
      GeneMixObj <- GeneMixObj_selected
      step=F
      if(silent==F)
        cat(' ')
    }else{
      return(GeneMixObj)
    }
  }



}

#'
#'
#' A more general forward model for geneObject
#'
#'
GeneForward <- function(GeneMixObj,
                        candidates = NULL,
                        dupl = NULL,
                        C = 4,
                        restricted=T,
                        silent = F) {

  M <- ncol(GeneMixObj$LikObj$UX)
  find <- c()
  if(is.null(candidates))
    candidates = 1:M
  if(is.null(dupl))
    dupl <- findDuplicate(GeneMixObj$LikObj$X)
  GeneMixObj <- mixedModelOptim(GeneMixObj,
                           dupl = dupl,
                           restricted=restricted)
  if(GeneMixObj$tauOff==T)
    GeneMixObj$tau = 0
  # forward
  candidates <- setdiff(candidates, dupl[,2])
  n <- length(MixGeneObj$data$Y)
  GeneMixObj$BIC <- BICmixed.R(GeneMixObj$theta0,
                               GeneMixObj$LikObj,
                               C= C)
  for(i in 1:length(candidates)){ # liberal critical value
    if(silent==F){
      cat('*')
    }
    find <- c(GeneMixObj$find, candidates[i])
    GeneMixObj_new <- mixedModelOptim(GeneMixObj, find = find, dupl = dupl, restricted=restricted)
    GeneMixObj_new$BIC <- BICmixed.R(GeneMixObj_new$theta0,
                                     GeneMixObj_new$LikObj,
                                     find = find,
                                     C= C)
    if(GeneMixObj_new$BIC  < GeneMixObj$BIC ){
      GeneMixObj <- GeneMixObj_new
      cat(i,' ')
    }


  }
  return(GeneMixObj)
}



#'
#' fitting a mixed model with selection on X usign forward selection
#' @param GeneMixObj
#' @param markers    - (m x 1) length of the m chromosnes
#' @param dCM        - (double) distance in centiMorgan
#' @param dupl       - (c x 1) duplicated rows of X (found in GeneMixObj$LikObj$X)
mixedModelForwardBackward <- function(GeneMixObj,
                                      markers ,
                                      dCM =1,
                                      dupl = NULL,
                                      liberal = TRUE,
                                      qval = NULL,
                                      alpha = 0.05) {


  M <- ncol(GeneMixObj$LikObj$UX)
  find <- c()
  GeneMixObj <- mixedModel(GeneMixObj, dupl = dupl)
  if(GeneMixObj$tauOff==T)
    GeneMixObj$tau = 0
  t <- GeneMixObj$t
  maxv <- max(abs(t))
  if(is.null(qval))
    qval <- qnorm(1 - 0.05 / (M*2) * 5)
  # forward
  while (maxv > qval) { # liberal critical value
    find <- c(GeneMixObj$find, setdiff(which(abs(t) == maxv), dupl[,1]))
    GeneMixObj <- mixedModel(GeneMixObj, find = find, dupl = dupl)
    t <- GeneMixObj$t
    if(!is.null(dupl)){
      while(sum(abs(t[dupl[,1]]-t[dupl[,2]]))>0)
        t[dupl[,1]] <- t[dupl[,2]]
    }
    maxv <- max(abs(t[setdiff(1:M, c(find,dupl[,1]))]), na.rm = TRUE)
  }
  beta_critval = 0.02
  if(GeneMixObj$tauOff==FALSE){
    beta_critval = NULL
  }
  crit <- calculateCrit(t, markers, d = dCM, beta = beta_critval,alpha = alpha )

  find <- GeneMixObj$find
  GeneMixObj$find_forward <- find
  # backward step
  if (length(find) > 0) {
    wmin <- 0
    t.back <- abs(t[find])
    while (length(find) && min(t.back, na.rm = TRUE) < crit) {
      wmin <- which.min(t.back)
      find <- find[-wmin]
      GeneMixObj <- estimateGeneMix(GeneMixObj, find = find)
      WhiteData <- whiteData(GeneMixObj, XOut = F, Xextra = GeneMixObj$LikObj$UX[, find])
      Xm <- cbind(WhiteData$Xm,WhiteData$Xextra)
      fit <- RcppEigen::fastLmPure(Xm, WhiteData$Y)

      if (length(find)) t.back <- abs(fit$coef / fit$se)[(ncol(WhiteData$Xm)+1):ncol(Xm)]
    }
    if (wmin){
      GeneMixObj <- mixedModel(GeneMixObj, find = find)
    }else{
      GeneMixObj <- estimateGeneMix(GeneMixObj, find = find)
      GeneMixObj <- mixedModel(GeneMixObj, find = find)
    }
  }

  GeneMixObj$crit <- crit
  return(GeneMixObj)
}
