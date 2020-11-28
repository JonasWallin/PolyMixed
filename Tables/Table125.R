##
# power simulation for mixed effect model
# to create Table 1: set.seed(2) q=1000, ns = c(200,400), type = 2, strong =1, use_mu0   = TRUE
#            for mixed-PC: use_mu0   = FALSE, nPC = number of Prinicpal components
#            for reg-PC  : use_mu0   = FALSE, use_tau=FALSE nPC = number of Prinicpal components
# to create Table 2: set.seed(2) q=1000, ns = c(200,400), type = 1, strong =1, use_mu0   = TRUE
#            for mixed-PC: use_mu0   = FALSE, nPC = number of Prinicpal components
#            for reg-PC  : use_mu0   = FALSE, use_tau=FALSE nPC = number of Prinicpal components
# to create Table 5: set.seed(2) q=1000, ns = c(200,400), type = 1, strong =0, use_mu0   = TRUE
#            for mixed-PC: use_mu0   = FALSE, nPC = number of Prinicpal components
#            for reg-PC  : use_mu0   = FALSE, use_tau=FALSE nPC = number of Prinicpal components
# D: 2019-04-01
##

save.file=F
graphics.off()
library(bigstep)
library(PolyMixed)
library(RcppEigen)
library(ggplot2)
library(bigstep)
set.seed(2)

# parameters----------------------------------------------------------------
type = 2                 # 1 - pure measurement error (\sigma^2I )
                         # 2 -  mixed effect     (\sigma^2I + XX^T \tau )
signal = 1               # 1 - strong signal , 0 weak signal, -1 no signal
model = 'simple'         # (no effect currently) 'simple', 'constant' (not implimented), 'PC'
nPC   = -1                # how many PC to use
use_mu0   = TRUE         # estimate mu_0
use_tau   = TRUE         #
q <- 1000                # number of simulations
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
if(model=='constant')
  nc <- chr
for(n in ns){
  if(type > 0){
    qtl.pos <- c(151, 825, 1275)
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
  }else{
    qtl.pos   <- c()
    beta.true <- c()
  }
  # --------------------------------------------------------------------------


  find.mix.backward <- list()
  beta.mix.backward <- list()
  find.mix.forward  <- list()
  beta.mix.forward  <- list()
  find.mix.tvalues  <- list()
  t.mix.forward.known <- matrix(0, q, M)
  t.mix.forward <- matrix(0, q, M)
  sigma.est.forward <- numeric(q)
  tau.est.forward <- numeric(q)
  crit.mix.forward <- numeric(q)
  ignoreTau <- numeric(q) # LRT test
  tau.forward       <- rep(0, q)
  sigma.forward     <- rep(0, q)
  mu.forward        <- rep(0, q)
  tau.backward      <- rep(0, q)
  sigma.backward    <- rep(0, q)
  mu.backward       <- rep(0, q)
  pb <- txtProgressBar(min = 0, max = q, style = 3)
  crit.vec <- rep(0, q)
  for(i in 1:q) {
    # Should we allow for duplitacte?
    X <- simulateDesignMatrix(chr, m, n, d)
    SVDX = svd(X)
    dupl <- findDuplicate(X)
    if(type == 2){
      y <- simulateTraits(X, q = 1, mu, tau, qtl = qtl, qtl2 = qtl2, qtl3 = qtl3)
    }else if(type==1){
      beta_          <- rep(0, dim(X)[2])
      beta_[qtl.pos] <- beta.true
      y <- X%*%beta_+ sigma * rnorm(dim(X)[1])
    }else if(type == 0){
        betas =  tau * rnorm(dim(X)[2])
        y <- X%*%(mu[1] + betas) + sigma * rnorm(dim(X)[1])
        #y <- simulateTraits(X, q = 1, mu, tau, qtl = 0, qtl2 = 0, qtl3 = 0)

    }
    # mixed model with forward and known tau and sigma
    MixGeneObj <- SetupGeneMix('Y ~ 1',
                               data = data.frame(Y=y),
                               X=X,
                               SVDX = SVDX,
                               meanMuOff = !use_mu0,
                               tauOff = !use_tau,
                               nPC = nPC)
    MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                    markers = markers,
                                    dupl = dupl)
    #res <- mixedModelForward(X, y, markers = markers, dupl = dupl, liberal = TRUE, SVDX = SVDX)
    t.mix.forward.known[i, ] <- MixGeneObj$t
    t.mix.forward[i, ]       <- MixGeneObj$t
    if(length(MixGeneObj$find)==0){
      find.mix.forward[[i]]    <- list()
      beta.mix.forward[[i]]    <- list()

    }else{
      find.mix.forward[[i]]    <- MixGeneObj$find
      beta.mix.forward[[i]]    <- MixGeneObj$beta[grep('X',names(MixGeneObj$beta))]

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

  cat('for (n = ',n,'):\n')
  print(round(summ.forward, 2))
  output <- list(summ.forward    = summ.forward,
                 tau.forward     = tau.est.forward,
                 sigma.est       = sigma.est.forward)
  if(save.file){
    if(nPC==0){
      save(output, file=paste("data/result_v2_",type,'_',n,'_',model,".RData",sep=""))
    }else{
      save(output, file=paste("data/result_v2_",type,'_',n,'_PC_',nPC,".RData",sep=""))
    }
  }
}
