###
# Figures for  t-values for trait 975,1190
#
# Date: 2018-12-30
###
save.fig= F
library(RcppEigen)
library(ggplot2)
library(bigstep)
seed = 4
set.seed(seed)
alpha_type = 2
library(PolyMixed)
# parameters------------------------------------------------------------
chr <- 10
# chromosomes
m <- 150                # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
q <- chr * m             # traits
n <- 200                 # observations
d <- 1                   # distance between markers (cM)
sigma <- 1               # sd for error term
tau <- 0.01              # sd for polygenic effect
mu_sd  <- 0.007          # ask gosia
cis <- c(0.5, 0.1)       # cis (mean, sd)
qtl <- c(1, 525, 1190)   # qtl (beta; position: marker, trait or traits)
qtl2 <- c(-1, 375, 979)
qtl3 <- c(1, 601, 62)
crit.liberal <- qnorm(1 - 0.05 / (M*2) * 5) # Biometrics correction
# ----------------------------------------------------------------------
if(alpha_type == 1){
  alpha = 0.05
}else{
  alpha = 0.05 / q
}

# design matrix and traits
X <- simulateDesignMatrix(chr, markers, n, d, seed = seed)

ys <- matrix(0, nrow=n, ncol=q)
traits <- 1:q

mus             <- rnorm(q, mean=0      , sd = mu_sd)
betas           <- rnorm(q, mean = cis[1], sd = cis[2])
for(trait in traits){
  set.seed(seed + trait)
  beta <- rep(0,dim(X)[2])
  beta[trait] <- betas[trait]
  ys[,trait]           <- X%*%(beta + mus[trait] + rnorm(q, 0, sd=tau)) + rnorm(n)
}
dupl <- findDuplicate(X)
SVDX = svd(X)



count <- 0
t_s <- matrix(0, nrow=q,ncol=q)
Design_fixed <- vector(length=q, mode="numeric")
Et_s <- matrix(0, nrow=q,ncol=q)
unitv   <- as.matrix(rep(1,n))
selected <- matrix(0, nrow=q,ncol=q)
XXt <- X%*%t(X)
trXXt <- sum(diag(XXt))
tXX_one <- t(X)%*%cbind(X,unitv)
tBB_full  <- t(cbind(X,unitv))%*%cbind(X,unitv)
tXX_one2 <- t(tXX_one)%*%tXX_one
for(trait in traits){
  if(trait%%10==0)
    cat('t = ',trait,'\n')


  X_t <- as.matrix(X[,trait])


  Ebeta <- rep(0,q)
  ES2   <- rep(0,q)
  beta <- rep(0,dim(X)[2])
  beta[trait] <- betas[trait]
  Xbeta <- X%*%beta
  X1    <- X%*%rep(1,q)
  EY <- Xbeta + X1*mus[trait]
  Z1 <- X1
  t.val <- rep(0,q)
  for(test in 1:q){
    if(test == trait){
      B = cbind(X_t,unitv)
      #Temp2 <- tXX_one[,c(test,q+1)]
      tBB <- tBB_full[c(test,q+1),c(test,q+1)]
      TT2   <- tXX_one2[c(test,q+1),c(test,q+1)]
    }else{
      B = cbind(X[,test],X_t,unitv)
      tBB <- tBB_full[c(test,trait,q+1),c(test,trait,q+1)]
      TT2   <- tXX_one2[c(test,trait,q+1),c(test,trait,q+1)]
    }
    #tBB <- t(B)%*%B
    #res <- qr(tBB)
    if(dim(B)[2]==3){
      iBB <- matrix(0,3,3)
      iBB[1,1] <- tBB[3,3] * tBB[2,2] - tBB[2,3]^2
      iBB[1,2] <- tBB[2,3]*tBB[1,3] - tBB[3,3] * tBB[1,2]
      iBB[2,1] <- iBB[1,2]
      iBB[1,3] <- (tBB[2,3] * tBB[1,2]  - tBB[2,2]*tBB[1,3])
      iBB[3,1] <- iBB[1,3]
      iBB[2,2] <- tBB[3,3] * tBB[1,1] - tBB[1,3]^2
      iBB[2,3] <- tBB[1,2] * tBB[1,3] - tBB[1,1] * tBB[2,3]
      iBB[3,2] <- iBB[2,3]
      iBB[3,3] <- tBB[1,1] * tBB[2,2] - tBB[1,2]^2
      detBB <- tBB[1,1]  * iBB[1,1] +
               tBB[2,1]  * iBB[1,2] +
                tBB[3,1] * iBB[1,3]
    }else{
      iBB <- tBB
      iBB[1,1] <- tBB[2,2]
      iBB[2,2] <- tBB[1,1]
      iBB[1,2] = -iBB[1,2]
      iBB[2,1] = iBB[1,2]
      detBB <- tBB[1,1]*tBB[2,2] - 2*tBB[1,2]

    }
    if(abs(detBB) <10^-10){
      Ebeta[test]  <- 0
      ES2[test]    <- Inf
    }else{
    iBB  <-  iBB/detBB
    #iBB  <- solve(tBB)
    #PB   <- B%*%iBB%*%t(B)
    BY  <- t(B)%*%EY
    mu_B <- iBB%*%BY
    Design_fixed[test] <- (iBB%*%(t(B)%*%Z1))[1]
    #Temp <- t(X)%*%(B%*%iBB)
    #Temp <- Temp2%*%iBB

    #Sigma_B <- tau^2 * (t(Temp)%*%Temp) + iBB[1,1] * 1
    Sigma_B <- tau^2 * iBB%*%TT2%*%iBB + iBB[1,1] * 1

    Ebeta[test] <- sqrt(Sigma_B[1,1] * 2/pi) * exp(-mu_B[1]^2/(2*Sigma_B[1,1])) -mu_B[1]*(1-2*pnorm(mu_B[1]/sqrt(Sigma_B[1,1])))
    p <- dim(B)[2]
    #ERSS <- (n-p)*1 + sum(EY^2)  - sum(BY*mu_B) + tau^2 * (trXXt - c(PB)%*%c(XXt))
    ERSS <- (n-p)*1 + sum(EY^2)  - sum(BY*mu_B) + tau^2 * (trXXt - c(iBB)%*%c(TT2))

    ES2[test] <- iBB[1,1]*ERSS/(n-p)
    beta_hat <- iBB%*%t(B)%*%ys[,trait]
    sd <- sum((ys[,trait]-B%*%beta_hat)^2)/(n-p)
    t.val[test] <-beta_hat[1]/sqrt(sd*iBB[1,1])
    }
  }
  ET <- Ebeta/sqrt(ES2)

  ##
  # calculating all t-values expect trait
  ##
  #X_mt <- X[,-trait]
  #y    <- ys[,trait]
  #t_mt <- apply(X_mt, 2, function(x){fit <- fastLmPure(cbind(x,X_t,unitv), y); fit$coefficients[1] / fit$stderr[1]})

  ##
  # calculating t-value for trait
  ##
  #fit_t <- fastLmPure(as.matrix(X_t,unitv), y)
  #if(trait==1){
  #  t.val <- c( fit_t$coefficients[1] / fit_t$stderr[1], t_mt)
  #}else if(trait==q){
  #  t.val <- c(t_mt,fit_t$coefficients[1] / fit_t$stderr[1])
  #}else{
  #  t.val <- c(t_mt[1:(trait-1)], fit_t$coefficients[1] / fit_t$stderr[1], t_mt[trait:length(t_mt)])
  #}


  t_s[, trait]       <- t.val
  t_s[ET==0, trait]  <- 0
  Et_s[, trait]      <- ET
  crit.simple <- calculateCrit(abs(t_s[,trait]), markers, alpha = alpha )
  selected[, trait]  <- abs(t_s[,trait]) > crit.simple

}
t_s[ET==0]=0
if(save.fig==FALSE){
  x11()
}else{
  pdf("Figure6_meantrait.pdf", w = 5, h = 3)
}
crit.ber <- calculateCrit(rowMeans(t_s), markers, alpha = 0.05 / q)
fig <- plotTrait(rowMeans(abs(t_s)), markers, crit.ber, addline = rowMeans(abs(Et_s)))
addline2 <- apply(abs(t_s),1, function(x) quantile(x,c(0.75)))
fig <- fig + geom_line(aes(Marker, addline2), col = "blue", size = 0.3)
#fig <- fig + geom_line(aes(Marker, Design_fixed/sd(Design_fixed)), col = "green", size = 0.3)
crit.ber <- calculateCrit(rowMeans(t_s), markers, alpha = 0.05 )
fig <- fig + geom_hline(yintercept = crit.ber)
fig <- fig + theme(axis.title.y=element_text(size=20), axis.title.x = element_text(size=20))
print(fig)
if(save.fig)
  dev.off()


if(save.fig){
  pdf(paste("Figure6_simple_2008_alpha",5,".pdf",sep=''), w = 5.21, h = 3)
}else{
  x11(w=5,h=3)
}


selected_simple <- matrix(0, nrow  = q, ncol  = q)
for(trait in traits){
  t_s[Et_s[,trait]==0, trait]  <- 0
  crit.simple <- calculateCrit(abs(t_s[,trait]), markers, alpha = alpha )
  selected_simple[abs(t_s[,trait]) > crit.simple, trait] = 1
}
print(plotHot(t(selected_simple), crit = 0.5, markers)+
        theme(panel.border = element_rect(fill = NA),axis.title.y=element_text(size=20), axis.title.x = element_text(size=20)))

if(save.fig)
  dev.off()
