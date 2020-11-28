##
#
# outdated models
#
##


iterativeForwardSelectR<-function(Y,
                                  X,
                                  model = c('simple','constant'),
                                  kmax = NULL,
                                  nc   = NULL,
                                  X_rem  = NULL,
                                  critval = qnorm(1-(1-0.95)/2),
                                  likelihood_method = FALSE,
                                  silent=1){
  #'
  #' Iterative procedure for selecting covariates uder base model:
  #'   Y       \sim N(X_k\beta + X \gamma, \sigma_Y^2 I_{n \times n})
  #'   \gamma  \sim N(u, \tau^2 I_{p \times p})
  #'   u       = [u_1,\ldots,u_{n_c}]
  #'   u_i     \sim N(\mu_0, \sigma^2_u 1 %*% 1 ^T)
  #' In e
  #'
  #' @param Y     - (n x 1) observations
  #' @param X     - (n x p) covariates
  #' @param model - (string)
  #'                 simple   - simpliefed model (sigma^2_u = 0)
  #'                 constant - the full model
  #' @param Kmax  - (int) maximum number of iterations
  #' @param nc    - (int) number of chromosones
  #' @param X_rem - (p< x 1) number of columns X that should not be checked (like duplicates)
  #' @param critval  - (double) value to minimum accepted level
  #' @param likelihood_method - (bool) use highest likelihood to select coeff
  #' @param silent (int) 0 - no output
  #'                     1 - result of each loop
  #'                     2 - with stars for iteration
  #'
  #' @return beta - (p x 1) vector with the selected betas (and values)
  #'

  model <- match.arg(model)

  p <- dim(X)[2]
  if(is.null(kmax))
    kmax = min(p, 20)
  Check <- ArgumentCheck::newArgCheck()

  if (kmax  > p)
    ArgumentCheck::addError(
      msg = "kmax must be smaller then p",
      argcheck = Check
    )

  if(model=='simple'){
    theta0 <- c(log(sd(Y)), log(sd(Y))-2)
  }else{
    if(is.null(nc))
      ArgumentCheck::addError(
        msg = "if model not simple must set nc",
        argcheck = Check
      )

    theta0 <- c(log(sd(Y)), log(sd(Y))-2, log(sd(Y))-2)
  }
  ArgumentCheck::finishArgCheck(Check)
  beta <- rep(0, p)

  # the iteration procedure
  SVDX <- svd(X)
  likObj <- loglikSetup(Y, NULL, X, nc = nc, SVDX = SVDX)
  lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
  res <- optim(theta0, lambda)
  Sigma_U <- getSigmaUY(exp(res$par), likObj)
  if(is.vector(Sigma_U)==F){
    iSY <- solve(Sigma_U,likObj$Yu)
  }else{
    iSY <- likObj$Yu / Sigma_U
  }
  uX <- t(likObj$SVDX$u)%*%X
  t_values <- c()
  find     <- c()
  beta_found <- list()
  for(k in 1:kmax){
    if(likelihood_method == FALSE){
      Beta_vec <- getTValues_R(setdiff(which(beta==0), X_rem),
                               likObj$Xu,
                               uX,
                               Sigma_U,
                               iSY)
      Beta_vec <- abs(Beta_vec)
    }else{
      for(i in which(beta==0)){
        if(silent>1)
          cat('*')
        beta_i <- beta
        beta_i[i] <- 1
        likObj_i <- loglikSetup(Y, Xk = X[,beta_i == 1], X, nc = nc, SVDX = SVDX)
        lambda_i <- function(theta0){-loglikSimpleR(theta0, likObj_i)}
        res_i <- optim(res$par, lambda_i)
        Beta_vec[i] <- (-res_i$value)
      }
    }
    ind_beta <- which(Beta_vec == max(Beta_vec))[1]
    if(silent > 0){
      cat('selected in iteration ',k,' coeff_nr = ',ind_beta,',')
      cat(' tau =' ,format(exp(res$par[2]), scientific = T, digits=2),',')
    }
    beta[ind_beta] = 1
    likObj <- loglikSetup(Y, Xk = X[,beta == 1], X, nc = nc, SVDX = SVDX)
    lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
    res <- optim(res$par, lambda)
    betaCov <- getSigmaBeta(res$par, likObj)
    beta_mu <- getBetaMu0(res$par, likObj)
    nb <-  sum(beta == 1)
    # create marginal confidence interval for all betas
    sds <- sqrt( diag(as.matrix( betaCov[1:nb,1:nb])))
    Z_Beta <- beta_mu$beta/sds #+ qnorm(1-(1-p.val)/2)* cbind(sds,-sds)
    Z_small <- min(abs(Z_Beta))
    # one of the betas was not signficant!
    if(silent > 0)
      cat(' Z_min = ', round(Z_small,2),'\n')
    if( Z_small < critval){
      beta[ind_beta] = 0
      break
    }
    t_values <- c(t_values, Z_small)
    find     <- c(find, ind_beta)
    beta_found[[length(find)]] <- beta_mu$beta
  }
  likObj <- loglikSetup(Y, Xk = X[,beta == 1], X, nc = nc, SVDX = SVDX)
  lambda <- function(theta0){-loglikSimpleR(theta0, likObj)}
  res <- optim(res$par, lambda)
  beta_mu <- getBetaMu0(res$par, likObj)
  beta[beta == 1] <- beta_mu$beta
  return(list(betas = beta,
              t_values = t_values,
              find = find,
              beta_found = beta_found))
}

getTValues_R <- function(index, Xu, uX, Sigma_U, Y){

  p <- ncol(uX)
  ni <- length(index)
  tvalues <- rep(0, p)
  n_p <- ncol(Xu) + 1
  for(j in 1:ni){
    i = index[j]
    Xcov <- cbind(Xu, uX[,i])
    if(is.vector(Sigma_U)){
      XtiSX <- t(Xcov)%*%(Xcov / Sigma_U)
    }else{
      XtiSX <- t(Xcov)%*%solve(Sigma_U, Xcov)
    }

    CovBeta <- solve(XtiSX)
    beta_est <- CovBeta%*%( t(Xcov)%*%Y)
    tvalues[i] <- beta_est[n_p]/sqrt(CovBeta[n_p, n_p])
  }
  return(tvalues)
}
