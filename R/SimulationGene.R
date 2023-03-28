#
# simulation of different types
#
# D: 2019-11-05
#



#'
#' @param type (int) 1,2,3
#'                  1 - pure mixed effect
#'                  2 - varying effect over chromosones
#'                  3 - varying effect over each chromosone
#' @param k     (int) number of chromosone
#' @param n     (int) number of observations
#' @param L     (int) length of chromosne (i.e number locus)
#' @param tau   (double) beta variation
#' @param mu    (double) baseline beta
#' @param sigma (double) standard deviation of measurement error
#' @param betas (k*L) which betas are nonzero
#' @param betasVal (k*L or 1x1) which are the values of betas
#' @export
simulation_GeneData <- function(type,
                                k = 10,
                                n = 500,
                                L = 100,
                                tau = 0.02,
                                mu0 =0.005,
                                sigma = 1,
                                betas = NULL,
                                betasVal = NULL){

  # setting up X
  d <- 1 #odleglosc miedzy sasiednimi markerami w cM
  r <- 0.5*(1-exp(-0.02*d)) #p-stwo rekombinacji miedzy sasiednimi markerami
  X <- matrix(0, nrow=n, ncol= k * L)
  markers <- rep(L, k)
  for (i in 1:k)
  {
    P  <- rbinom(n,1,0.5);#genotypy w pierwszym markerze
    P  <- runif( n ) < 0.5
    R  <- matrix(runif( n * (L - 1) ) < r, nrow = n, ncol = L - 1)
    RP <- t(apply(cbind2(P,R), 1, cumsum))
    X[,seq(from=(i-1)*L+1, to =(i*L))]<-RP%%2;#finalna macierz genotypow
  }

  if(type == 1){
    mu = rep(0, L*k)
  }
  if(type == 2){
    c  <- rnorm(k)/10;
    mu <- rep(c,each = L)
  }
  if(type == 3){
    mu                 <- rep(0, L*k)
    c                  <- rnorm(5,0,1)/10;
    mu[1:(5*L)]        <- rep(c, each = L)
    mu[1:L]            <- 0.1 * sin(pi*seq(1:L)/100)
    mu[(L+1):(2*L)]    <- 0.1*sin(pi*seq(1:L)/100/2)
    mu[(2*L+1):(3*L)]  <- 0.1*sin(pi*seq(1:L)/100*2)

  }
  beta <- mu0 + mu  + tau * rnorm(L*k);
  Y<-X%*%beta + sigma * rnorm(n)


  beta  <- rep(0, L * k)
  if(is.null(betas)){
    index <- round(c(0.25,0.5,0.75) * (L*k))
  }else{
    index <- rep(0, L * k)
    index[betas] = 1
    index <- index ==1
  }
  if(is.null(betasVal)){
    beta[index ] = -0.8
  }else{
    beta[index] = betasVal
  }
  Y <- Y + X%*%beta
  return(list(Y = Y, X = X, beta = beta, k = k, mu = mu, markers = markers))
}


###
#' simulation of Xs through selection
#' input
#' @param beta         - (n x 1) coeffients
#' @param X            - (N x n x 2) original coeffients
#' @param ngen         - (int)  number of generations to simulate
#' @param p            - (N x 1) probability of generating offspring
#' @param lambda_cross - (double) intensty of number of crossing
#' @export
###
simulation_GeneDataSelectionR <-function(beta, X, ngen, p, lambda_cross){

  n <- dim(X)[2]
  N <- dim(X)[1]
  for(g in 1:ngen){
    Xs <- X[,,1] + X[,,2]
    Y <- Xs%*%beta + rnorm(n=N)
    oY <- order(Y, decreasing=T)
    sF <-sample(oY,
                N,
                replace=T,
                prob=p)
    XF <- X[sF,,]
    sM <-sample(oY,
                N,
                replace=T,
                prob=p)
    XM <- X[sM,,]
    CrossF <- rpois(N, lambda=  lambda_cross)
    CrossM <- rpois(N, lambda=  lambda_cross)
    XF1 <- XF[,,1]
    XF2 <- XF[,,2]
    XM1 <- XM[,,1]
    XM2 <- XM[,,2]

    for(i in 1:N){
      splits <- diff(c(0,sort(sample(2:(n-1),CrossF[i],replace = F)),n))
      index <- ((rep(1:(CrossF[i]+1), splits)+(runif(1)<0.5) ) %%2)==T
      XF1[i, index] <- XF2[i, index]
      splits <- diff(c(0,sort(sample(2:(n-1),CrossM[i],replace = F)),n))
      index <- ((rep(1:(CrossM[i]+1), splits)+(runif(1)<0.5) ) %%2)==T
      XM1[i, index] <- XM2[i, index]
    }
    if(runif(1)>0.5){
      X[,,1] <- XF1
      X[,,2] <- XM1
    }else{
      X[,,1] <- XM1
      X[,,2] <- XF1

    }
  }
  return(X)
}


###
#' simulation of Xs through selection
#' input
#' @param beta         - (n x 1) coeffients
#' @param X            - (N x n x 2) original coeffients
#' @param ngen         - (int)  number of generations to simulate
#' @param p            - (N x 1) probability of generating offspring
#' @param lambda_cross - (double) intensty of number of crossing
#'
#' @export
###
simulation_GeneDataSelection <-function(beta, X, ngen, p, lambda_cross){

  n <- dim(X)[2]
  N <- dim(X)[1]
  for(g in 1:ngen){
    Xs <- X[,,1] + X[,,2]
    Y <- Xs%*%beta + rnorm(n=N)
    oY <- order(Y, decreasing=T)
    index <- runif(N)< 0.5
    sF <-sample(oY[index],
                N,
                replace=T,
                prob=p[1:sum(index)])
    XF <- X[sF,,]
    sM <-sample(oY[index==0],
                N,
                replace=T,
                prob=p[1:sum(1-index)])
    XM <- X[sM,,]
    CrossF <- rpois(N, lambda=  lambda_cross)
    CrossM <- rpois(N, lambda=  lambda_cross)
    XF1 <- XF[,,1]
    XM1 <- XM[,,2]

    crossSelect(XF1, XM1, CrossF);


    if(runif(1)>0.5){
      X[,,1] <- XF1
      X[,,2] <- XM1
    }else{
      X[,,1] <- XM1
      X[,,2] <- XF1

    }
  }
  return(X)
}



#' mixing two populations
#'
#' @param X1            - (N1 x n x 2) original X pop1
#' @param X2            - (N2 x n x 2) original X pop2
#' @param ngen         - (int)  number of generations to simulate
#' @param lambda_cross - (double) intensty of number of crossing
#' @export
simulation_mixingR <- function(X1,
                              X2,
                              lambda_cross,
                              generations = 1
){
  N1 <- dim(X1)[1]
  N2 <- dim(X2)[1]


  N <- N1 + N2
  X <- array(runif(n*N*2), dim=c(N, n, 2))
  X[1:N1,,]             <- X1
  X[(N1+1):(N1+N2), , ] <- X2
  Z <- X
  Z[1:N1,,]             <- 0
  Z[(N1+1):(N1+N2), , ] <- 1
  for(g in 1:generations){

    sF <-sample(N,
                replace=T)
    XF <- X[sF,,]
    ZF <- Z[sF,,]

    sM <-sample(N,
                replace=T)
    XM <- X[sM,,]
    ZM <- Z[sM,,]
    CrossF <- rpois(N, lambda=  lambda_cross)
    SF <- sample(1:2,N,replace=T)
    SM <- sample(1:2,N,replace=T)


    for(i in 1:N){
      X[i,,1] <- XF[i,,SF[i]]
      X[i,,2] <- XM[i,,SM[i]]
      Z[i,,1] <- ZF[i,,SF[i]]
      Z[i,,2] <- ZM[i,,SM[i]]
      splits <- diff(c(0,sort(sample(2:(n-1),CrossF[i],replace = F)),n))
      index <- ((rep(1:(CrossF[i]+1), splits)+(runif(1)<0.5) ) %%2) ==T
      XM_i <-  X[i,index,2]
      XF_i <-  X[i,index,1]
      ZM_i <-  Z[i,index,2]
      ZF_i <-  Z[i,index,1]
      X[i, index, 1] <- XM_i
      X[i, index, 2] <- XF_i
      Z[i, index, 1] <- ZM_i
      Z[i, index, 2] <- ZF_i
    }
  }
  return(list(X=X,Z=Z))
}

#' mixing two populations
#'
#' @param X1            - (N1 x n x 2) original X pop1
#' @param X2            - (N2 x n x 2) original X pop2
#' @param ngen         - (int)  number of generations to simulate
#' @param lambda_cross - (double) intensty of number of crossing
#' @export
simulation_mixing <- function(X1,
                              X2,
                              lambda_cross,
                              generations = 1
){
  N1 <- dim(X1)[1]
  N2 <- dim(X2)[1]


  N <- N1 + N2
  X <- array(runif(n*N*2), dim=c(N, n, 2))
  X[1:N1,,]             <- X1
  X[(N1+1):(N1+N2), , ] <- X2
  Z <- X
  Z[1:N1,,]             <- 0
  Z[(N1+1):(N1+N2), , ] <- 1
  X1g <- X[,,1]
  X2g <- X[,,2]
  Z1g <- matrix(as.integer(0), nrow = N, ncol = n)
  Z2g <- matrix(as.integer(0), nrow = N, ncol = n)
  Z1g <- as.matrix(Z[,,1])
  Z2g <- as.matrix(Z[,,2])
  mode(Z1g) <-  "integer"
  mode(Z2g) <-  "integer"
  mixing_population(X1g,
                    X2g,
                    Z1g,
                    Z2g,
                    generations,
                    lambda_cross);
  X[,,1] <- X1g
  X[,,2] <- X2g
  Z[,,1] <- Z1g
  Z[,,2] <- Z2g
  return(list(X=X,Z=Z))
}
