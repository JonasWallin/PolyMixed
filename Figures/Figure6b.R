  ###
# Figures for  t-values for trait 975,1190
#
# Date: 2018-12-30
###
rm(list=ls())
graphics.off()
save.fig= F
library(RcppEigen)
library(ggplot2)
library(bigstep)
seed = 4
set.seed(seed)
alpha_type = 2 # 1-0.05 , 2- 0.05/traits
library(PolyMixed)
# parameters------------------------------------------------------------
chr <- 10                # chromosomes
m <- 150                # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
q <- chr * m             # traits
n <- 200                 # observations
d <- 1                   # distance between markers (cM)
sigma <- 1               # sd for error term
tau <- 0.01              # sd for polygenic effect
mu <- rnorm(q, 0, 0.007) # mean for polygenic effect
cis <- c(0.5, 0.1)       # cis (mean, sd)

crit.liberal <- qnorm(1 - 0.05 / (M*2) * 5) # Biometrics correction
# ----------------------------------------------------------------------
if(alpha_type == 1){
  alpha = 0.05
}else{
  alpha = 0.05 / q
}


# heritability
(mu.mean <- mean(abs(mu))) # average |mu| should be comparable with tau
cov_m <- function(i, j) exp(-2 * abs(i - j) * d/100)
sum.cov <- 0
for (i in 1:m) for (j in 1:m) sum.cov <- sum.cov + cov_m(i, j)
1 - sigma^2 / (chr*sum.cov*mu.mean^2 + M*tau^2 + sigma^2) # H^2

# design matrix and traits
X <- simulateDesignMatrix(chr, markers, n, d, seed = seed)
y <- matrix(0, nrow=n, ncol=q)
mus             <- rnorm(q, mean=0      , sd = 0.007)
betas           <- rnorm(q, mean = cis[1], sd = cis[2])
for(trait in 1:q){
  set.seed(seed + trait)
  beta <- rep(0,dim(X)[2])
  beta[trait] <- betas[trait]
  y[,trait]           <- X%*%(beta + mus[trait] + rnorm(q, 0, sd=tau)) + rnorm(n)
}
#
#for(trait in traits){
#  set.seed(seed + trait)
#  beta <- rep(0,dim(X)[2])
#  beta[trait] <- betas[trait]
#  y[,trait]           <- X%*%(beta + mus[trait] + rnorm(q, 0, sd=tau)) + rnorm(n)
#}
dupl <- findDuplicate(X)

traits <- 1:q
count <- 0
selected        <- matrix(0, nrow  = q, ncol  = q)
selected_simple <- matrix(0, nrow  = q, ncol  = q)
SVDX = svd(X)
for(trait in traits){
  cat('*')
  count = count +1
  MixGeneObj <- SetupGeneMix('Y ~ 1',
                             data = data.frame(Y=y[, trait]),
                             X=X,
                             SVDX = SVDX,
                             meanMuOff = F,
                             tauOff = F)
  MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                          markers,
                                          dupl = dupl,
                                          alpha = alpha)
  selected[trait,MixGeneObj$find] <- 1
  # simple model
  MixGeneObj <- SetupGeneMix('Y ~ 1',
                             data = data.frame(Y=y[, trait]),
                             X=X,
                             SVDX = SVDX,
                             meanMuOff = T,
                             tauOff = T)
  MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                          markers,
                                          dupl = dupl,
                                          alpha = alpha)
    selected_simple[trait, MixGeneObj$find] <- 1
}




if(save.fig){
  pdf(paste("Figure6_hotspot_mixture_alpha",alpha_type,".pdf",sep=''), w = 5.21, h = 3)
}else{
  x11(w=5,h=3)
}
print(plotHot(selected, crit = 0.5, markers)+
      theme(axis.title.y=element_text(size=20), axis.title.x = element_text(size=20))
      )
if(save.fig)
  dev.off()

