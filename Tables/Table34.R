##
# power simulation for mixed effect model
# to create Table 3: set.seed(2) q=1000, table =2 ,gives sparse mixed effect
# to create Table 4: set.seed(2) q=1000, table =1 ,gives laplace distribution
# D: 2020-02-10
##

rm(list=ls())
save.file=F
graphics.off()
library(bigstep)
library(PolyMixed)
library(RcppEigen)
library(ggplot2)
library(bigstep)
set.seed(2)
table <- 2        # which table
q     <- 1000     # number of simulations
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
nPcs     = c(0     ,3     ,10    ,0    ,0     ,3     ,10)
use_mu0s = c(FALSE ,FALSE ,FALSE ,TRUE ,FALSE ,FALSE ,FALSE)
use_taus = c(FALSE ,FALSE ,FALSE ,TRUE ,TRUE  ,TRUE  ,TRUE)
# parameters----------------------------------------------------------------
chr <- 10                # chromosomes
m <- 150                 # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
ns <- c(200,400)         # observations
d <- 1                   # distance between markers (cM)
sigma <- 1               # sd for error term
tau <- 0.01              # sd for polygenic effect
mu <- rep(0.004, q)      # mean for polygenic effect




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

t.mix.forward.known <- matrix(0, q, M)
t.mix.forward <- matrix(0, q, M)
sigma.est.forward <- numeric(q)
tau.est.forward <- numeric(q)
crit.mix.forward <- numeric(q)
ignoreTau <- numeric(q) # LRT test
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

    summ.forward <- summarizeFind.new(find.mix.forward, beta.mix.forward, beta.true, qtl.pos, maxCM = 15)
    #summ.back    <- summarizeFind.new(find.mix.backward, beta.mix.backward, beta.true, qtl.pos, maxCM = 15)

    #cat('for (n = ',n,'):\n')
    #print(round(summ.forward, 2))
    output <- list(summ.forward    = summ.forward,
                   tau.forward     = tau.est.forward,
                   sigma.est       = sigma.est.forward)
    string_out <- paste(string_out,'&',model,paste('& ',round(summ.forward[1,1:5],3),sep='',collapse = ' '),sep=' ')
    string_out <- paste(string_out,'& ',round(mean(tau.est.forward)/sqrt(thin),3),sep='')
    string_out <- paste(string_out,' \\\\ \n& ')
    if(nPC > 0)
      string_out <- paste(string_out,'$(PC=',nPC,')$',sep='')
    if(use_mu0 == FALSE & use_tau == TRUE & nPC == 0)
      string_out <- paste(string_out,'$(\\mu=0)$',sep='')
    for(s in 1:3){
      string_out <- paste(string_out,' & (',round(summ.forward[2,s],2),',',round(summ.forward[3,s],2),')',sep='')
    }
    string_out <- paste(string_out,' & & \\\\ \n',sep='')
    cat('\n')
    cat(string_out)
    if(save.file){
      if(table==1){
       write(string_out, file = paste("data/table_laplace_",table,".txt",sep=''), append=TRUE)
      }else{
        write(string_out, file = paste("data/table_subset_",table,".txt",sep=''), append=TRUE)
      }
    }
    string_out <- ''
  }

}
