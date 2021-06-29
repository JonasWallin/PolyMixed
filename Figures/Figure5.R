###
# test how positive correlates with size of polygenic effect
###
save.fig=T
library(bigstep)
library(PolyMixed)
library(RcppEigen)
library(ggplot2)
library(bigstep)
set.seed(6)
table <- 2        # which table
q     <- 100     # number of simulations
#thin =  how often to observe X (i.e. cm distnace)
#type =  1 - pure measurement error (\sigma^2I )
#        2 -  mixed effect     (\sigma^2I + XX^T \tau )
#signal = 1 - strong signal , 0 weak signal, -1 no signal

if(table==1){
  thin    = 1
  signal  = 1
}else{
  thin = 1
  signal = 1
  active_propotion = 0.2
}
nPcs     = c(0     )
use_mu0s = c(FALSE )
use_taus = c(FALSE)
# parameters----------------------------------------------------------------
chr <- 10                 # chromosomes
m <- 150                 # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
ns <- c(200)              # observations
d <- 1                   # distance between markers (cM)
sigma <- 1               # sd for error term
tau <- 0.01              # sd for polygenic effect
mu <- rep(0.004, q)      # mean for polygenic effect


t.mix.forward.known <- matrix(0, q, M)
t.mix.forward <- matrix(0, q, M)
sigma.est.forward <- numeric(q)
tau.est.forward <- numeric(q)
crit.mix.forward <- numeric(q)
ignoreTau <- numeric(q) # LRT test

nc = NULL
qtl.pos <- c(151, 824, 1274)
mag = 0.35
if(signal==1){
  mag = 0.5
}else if(signal==-1){
  mag = 0
}
qtl  <- c(mag, qtl.pos[1], 1) # qtl (beta; position: marker, trait or traits)
qtl2 <- c(mag, qtl.pos[2], 1)
qtl3 <- c(-mag, qtl.pos[3], 1)
beta.true <- c(qtl[1], qtl2[1], qtl3[1])

# --------------------------------------------------------------------------

for(n in ns){
  string_out <- paste('\n\n',n)
  Xs <- list()
  ys <- list()
  betas <- list()
  dupls <- list()
  SVDXs <- list()
  for(i in 1:q) {
    # Should we allow for duplitacte?
    X <- simulateDesignMatrix(chr, m, n, d)

    Z <- rnorm(dim(X)[2])
    if(table==1){
      beta_          <- rep(0, dim(X)[2])
      beta_[qtl.pos] <- beta.true
      V <- rgamma(dim(X)[2],1)
      y <- X%*%(beta_ + mu[1] + tau * sqrt(V)*Z) + sigma * rnorm(dim(X)[1])
    }else{
      index    = sample(1:dim(X)[2],ceiling(dim(X)[2]*active_propotion))
      V        = rep(0, dim(X)[2])
      V[index] = 1
      beta_    = V*mu[1]/active_propotion + tau/active_propotion * sqrt(V)*Z
      beta_[qtl.pos] <- beta.true
      y <- X%*%beta_ + sigma * rnorm(dim(X)[1])
    }

    #y <- simulateTraits(X, q = 1, mu, tau, qtl = 0, qtl2 = 0, qtl3 = 0)

    thin_index <- seq(1,dim(X)[2],by = thin)
    X          <- X[,thin_index]
    Xs[[i]]    <- X
    ys[[i]]    <- y
    betas[[i]] <- beta_
    dupls[[i]] <- findDuplicate(X)
    SVDXs[[i]] <- svd(X, nu = dim(X)[1])


  }
  for(count in 1:length(nPcs)){
    nPC <- nPcs[count]
    use_mu0 <- use_mu0s[count]
    use_tau <- use_taus[count]
    model = 'simple'
    if(use_tau)
      model = 'mixed'


    find.mix.backward <- list()
    beta.mix.backward <- list()
    find.mix.forward  <- list()
    beta.mix.forward  <- list()
    find.mix.tvalues  <- list()
    tau.forward       <- rep(0, q)
    sigma.forward     <- rep(0, q)
    mu.forward        <- rep(0, q)
    tau.backward      <- rep(0, q)
    sigma.backward    <- rep(0, q)
    mu.backward       <- rep(0, q)
    pb <- txtProgressBar(min = 0, max = q, style = 3)
    crit.vec <- rep(0, q)
    for(i in 1:q) {
      X <- Xs[[i]]
      y <- ys[[i]]
      dupl <- dupls[[i]]
      SVDX <- SVDXs[[i]]
      # mixed model with forward and known tau and sigma
      MixGeneObj <- SetupGeneMix('Y ~ 1',
                                 data = data.frame(Y=y),
                                 X=X,
                                 SVDX = SVDX,
                                 meanMuOff = !use_mu0,
                                 tauOff = !use_tau,
                                 nPC = nPC)
      MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                              markers = markers/thin,
                                              dCM = thin,
                                              dupl = dupl)
      #res <- mixedModelForward(X, y, markers = markers, dupl = dupl, liberal = TRUE, SVDX = SVDX)
      t.mix.forward.known[i, ] <- MixGeneObj$t
      t.mix.forward[i, ]       <- MixGeneObj$t
      if(length(MixGeneObj$find)==0){
        find.mix.forward[[i]]    <- list()
        beta.mix.forward[[i]]    <- list()

      }else{
        find.mix.forward[[i]]    <- thin_index[MixGeneObj$find]
        beta.mix.forward[[i]]    <- MixGeneObj$beta[grep('X',names(MixGeneObj$beta))]

      }
      if(0){
        if(MixGeneObj$tau>0){
          MixGeneObj <- SetupGeneMix('Y ~ 1',
                                     data = data.frame(Y=y),
                                     X=X,
                                     SVDX = SVDX,
                                     meanMuOff = !use_mu0,
                                     tauOff = !use_tau,
                                     nPC = nPC)
          MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                                  markers = markers/thin,
                                                  dCM = thin,
                                                  dupl = dupl)
        }
      }
      tau.est.forward[i]       <- MixGeneObj$tau
      sigma.est.forward[i]     <- MixGeneObj$sigma
      ignoreTau[i]             <- MixGeneObj$tauOff
      crit.mix.forward[i]      <- MixGeneObj$crit
      #    print(c(length(find.mix.forward[[i]])))

      setTxtProgressBar(pb, i)
    }


    #
  }
}
select_true <- c()
pos_range=5
for(pos in qtl.pos){
  select_true <- c(select_true, pos + -pos_range:pos_range)
}
k <- 58
false_result <- setdiff(find.mix.forward[[k]],select_true)
beta_false <- false_result
gamma_false <- c()
i <- 1
for( false_i in false_result){
  beta_false[i] <- beta.mix.forward[[k]][find.mix.forward[[k]]==false_result[i]]
  index <- intersect(false_i + -pos_range:pos_range,1:length(betas[[k]]))
  max_i <- min(which(abs(betas[[k]][index])== max(abs(betas[[k]][index])))) #largest of the true values
  gamma_false <- c(gamma_false, betas[[k]][index[max_i]])
  i <- i + 1
}
gammas <- betas[[k]][-qtl.pos]
gammas <- gammas[abs(gammas)>0]
gammas <- c(gammas,0)
vec <- cbind(rank(gammas), gammas)
index <- which(gammas%in% gamma_false)
if(save.fig)
  pdf('Figure5.pdf')
plot(vec[-index,1],vec[-index,2],pch=20,
     ylab=expression(gamma),
     xlab='Index in the ordered sequence of polygenic effects',
     ylim = c(1.05*min(c(vec[,2],beta_false)),1.05*max(c(vec[,2],beta_false))),
     xlim = c(0, length(gammas)),
     cex=0.3,
     col='blue')
abline(a=0,b=0)
points(vec[index,1], beta_false, pch=5,cex=1, col='red')
#abline(0)
#qtl.pos
if(save.fig)
  dev.off()


