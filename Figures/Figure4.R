###
# Figures for  t-values for trait 975,1190
#
# Date: 2018-12-30
###
graphics.off()
save.fig= T
alpha_res = 0.33
#source("../codes/misc/3-exp_matrix.R")
#source("../codes/misc/4-help_functions.R")

library(RcppEigen)
library(ggplot2)
library(bigstep)
library(reshape2)
seed  = 7
set.seed(seed)
col_v = 1.5
library(PolyMixed)
# parameters------------------------------------------------------------
chr <- 10                # chromosomes
m <- 150                 # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
q <- chr * m             # traits
n <- 200                 # observations
d <- 1                   # distance between markers (cM)
sigma <- 1               # sd for error term
tau <- 0.01              # sd for polygenic effect
mu <- 0.004              # mean for polygenic effect
beta.true <- c(0.5,0.5,-0.5)
qtl.pos   <- c(151, 825, 1275)
crit.liberal <- qnorm(1 - 0.05 / (M*2) * 5) # Biometrics correction
# ----------------------------------------------------------------------

# heritability
(mu.mean <- mean(abs(mu))) # average |mu| should be comparable with tau
cov_m <- function(i, j) exp(-2 * abs(i - j) * d/100)
sum.cov <- 0
for (i in 1:m) for (j in 1:m) sum.cov <- sum.cov + cov_m(i, j)
1 - sigma^2 / (chr*sum.cov*mu.mean^2 + M*tau^2 + sigma^2) # H^2

# design matrix and traits
X <- simulateDesignMatrix(chr, markers, n, d, seed = seed)
SVDX = svd(X)
dupl <- findDuplicate(X)
mareker_size = 20
tick_size   = 12
traits <- c(1190)
count <- 0
text_annonate1 <- c(151)
text_annonate2 <- c(825)
text_annonate3 <- c(1275)
for(trait in traits){
  count = count +1
  beta_          <- mu + tau * rnorm(dim(X)[2])
  beta_[qtl.pos] <- beta.true
  y <- X%*%beta_
  y <- y + sigma * rnorm(dim(X)[1])
  for(nPC in c(-1, 0, 3, 10)){

    MixGeneObj <- SetupGeneMix('Y ~ 1',
                               data = data.frame(Y=y),
                               X=X,
                               SVDX = SVDX,
                               meanMuOff = nPC>=0,
                               nPC = nPC)
    MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                            markers = markers,
                                            dupl = dupl)
    tvalues <- rep(0, dim(X)[2])
    beta_mixed <- rep(0, dim(X)[2])
    beta_mixed[MixGeneObj$find] <- MixGeneObj$beta[grep('X',names(MixGeneObj$beta))]
    MixGeneObj <- mixedModel(MixGeneObj, find = MixGeneObj$find, dupl = dupl, estPar = F)
    while(sum(abs(MixGeneObj$t[dupl[,1]]-MixGeneObj$t[dupl[,2]]), na.rm=T)>0){
      MixGeneObj$t[dupl[,1]] <- MixGeneObj$t[dupl[,2]]
    }
    if(save.fig){
      if(nPC==-1){
        pdf(paste("Figure4_mixed_non_cond_",trait,".pdf",sep=''), w = 5, h = 3)
      }else{
        pdf(paste("Figure4_PC_non_cond_",nPC,"_",trait,".pdf",sep=''), w = 5, h = 3)
      }
    }else{
      x11(w=5,h=3)
    }

    fig <- plotTrait(abs(MixGeneObj$t), markers, MixGeneObj$crit)
    for(f in MixGeneObj$find){
      fig <- fig + geom_vline(xintercept = f,
                              color = "red",
                              size=col_v,
                              alpha = alpha_res)
    }
    print(fig +
            annotate("text", x =  text_annonate3[count], y = -0.2 , label = "QTL3") +
            annotate("text", x =  text_annonate1[count], y = -0.2 , label = "QTL1") +
            annotate("text", x =  text_annonate2[count], y = -0.2 , label = "QTL2") +
            theme(axis.title.y=element_text(size=mareker_size),
                  axis.title.x = element_text(size=mareker_size),
                  plot.title = element_text(size=mareker_size),
                  axis.text.x = element_text(size=tick_size),
                  axis.text.y = element_text(size=tick_size)))
    if(save.fig)
      dev.off()



    ###
    # pure t-values
    ###

    MixGeneObj_0 <- mixedModel(MixGeneObj, find = NULL, dupl = dupl, estPar = F)
    if(save.fig){
      if(nPC==-1){
        pdf(paste("Figure4_mixed_tval_",trait,".pdf",sep=''), w = 5, h = 3)
      }else{
        pdf(paste("Figure4_PC_tval_",nPC,"_",trait,".pdf",sep=''), w = 5, h = 3)
      }
    }else{
      x11(w=5,h=3)
    }
    while(sum(abs(MixGeneObj_0$t[dupl[,1]]-MixGeneObj_0$t[dupl[,2]]))>0){
      MixGeneObj_0$t[dupl[,1]] <- MixGeneObj_0$t[dupl[,2]]
    }

    print(plotTrait(abs(MixGeneObj_0$t), markers, MixGeneObj$crit)+
            annotate("text", x =  text_annonate3[count], y = -0.2 , label = "QTL3") +
            annotate("text", x =  text_annonate1[count], y = -0.2 , label = "QTL1") +
            annotate("text", x =  text_annonate2[count], y = -0.2 , label = "QTL2")  +
            theme(axis.title.y=element_text(size=mareker_size),
                  axis.title.x = element_text(size=mareker_size),
                  plot.title = element_text(size=mareker_size),
                  axis.text.x = element_text(size=tick_size),
                  axis.text.y = element_text(size=tick_size)))
    if(save.fig)
      dev.off()


    ###
    # conditioning on other chromosone
    ###
    for(i in 1:chr){
      chr_index <- (m*(i-1)+1):(m*i)
      find_i <- setdiff(MixGeneObj$find, chr_index)
      MixGeneObj_i <- mixedModel(MixGeneObj, find = find_i, dupl = dupl, estPar = F)
      while(sum(abs(MixGeneObj_i$t[dupl[,1]]-MixGeneObj_i$t[dupl[,2]]))>0){
        MixGeneObj_i$t[dupl[,1]] <- MixGeneObj_i$t[dupl[,2]]
      }
      tvalues[chr_index]        <- MixGeneObj_i$t[chr_index]
    }
    if(save.fig){
      if(nPC==-1){
        pdf(paste("Figure4_mixed_",trait,".pdf",sep=''), w = 5, h = 3)
      }else{
        pdf(paste("Figure4_PC_",nPC,"_",trait,".pdf",sep=''), w = 5, h = 3)
      }
    }else{
      x11(w=5,h=3)
    }
    fig <- plotTrait(abs(tvalues), markers, MixGeneObj$crit)
    for(f in MixGeneObj$find){
      fig <- fig + geom_vline(xintercept = f,
                              color = "red",
                              size=col_v,
                              alpha = alpha_res)
    }
    print(fig +
            annotate("text", x =  text_annonate3[count], y = -0.2 , label = "QTL3") +
            annotate("text", x =  text_annonate1[count], y = -0.2 , label = "QTL1") +
            annotate("text", x =  text_annonate2[count], y = -0.2 , label = "QTL2")  +
            theme(axis.title.y=element_text(size=mareker_size),
                  axis.title.x = element_text(size=mareker_size),
                  plot.title = element_text(size=mareker_size),
                  axis.text.x = element_text(size=tick_size),
                  axis.text.y = element_text(size=tick_size)))
    if(save.fig)
      dev.off()
  }

  MixGeneObj <- SetupGeneMix('Y ~ 1',
                             data = data.frame(Y=y),
                             X=X,
                             SVDX = SVDX,
                             meanMuOff = T,
                             tauOff = T)
  MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                          markers,
                                          dupl = dupl)
  tvalues.simple <- rep(0, dim(X)[2])
  beta_simple <- rep(0, dim(X)[2])
  beta_simple[MixGeneObj$find] <- MixGeneObj$beta[2:length(MixGeneObj$beta)]

  t <-  MixGeneObj$t
  t[is.na(t)] <- 0
  while(sum(abs(t[dupl[,1]]-t[dupl[,2]]), na.rm=T)>0)
    t[dupl[,1]] <- t[dupl[,2]]
  if(save.fig){
    pdf(paste("Figure4_simple_noncond_",trait,".pdf",sep=''), w = 5, h = 3)
  }else{
    x11(w=5,h=3)
  }
  crit <- calculateCrit(t, markers)

  fig <- plotTrait(abs(t), markers, crit)
  for(f in MixGeneObj$find){
    fig <- fig + geom_vline(xintercept = f,
                            color = "red",
                            size=col_v,
                            alpha = alpha_res)
  }
  print(fig+
          annotate("text", x =  text_annonate3[count], y = -0.2 , label = "QTL3") +
          annotate("text", x =  text_annonate1[count], y = -0.2 , label = "QTL1") +
          annotate("text", x =  text_annonate2[count], y = -0.2 , label = "QTL2") +
          theme(axis.title.y=element_text(size=mareker_size),
                axis.title.x = element_text(size=mareker_size),
                plot.title = element_text(size=mareker_size),
                axis.text.x = element_text(size=tick_size),
                axis.text.y = element_text(size=tick_size)))
  if(save.fig)
    dev.off()



  ###
  # pure t-values
  ###

  MixGeneObj_0 <- mixedModel(MixGeneObj, find = NULL, dupl = dupl, estPar = F)
  if(save.fig){
    pdf(paste("Figure4_simple_tvalues_",trait,".pdf",sep=''), w = 5, h = 3)
  }else{
    x11(w=5,h=3)
  }
  while(sum(abs(MixGeneObj_0$t[dupl[,1]]-MixGeneObj_0$t[dupl[,2]]))>0){
    MixGeneObj_0$t[dupl[,1]] <- MixGeneObj_0$t[dupl[,2]]
  }
  print(plotTrait(abs(MixGeneObj_0$t), markers, MixGeneObj$crit)+
          annotate("text", x =  text_annonate3[count], y = -0.2 , label = "QTL3") +
          annotate("text", x =  text_annonate1[count], y = -0.2 , label = "QTL1") +
          annotate("text", x =  text_annonate2[count], y = -0.2 , label = "QTL2") +
          theme(axis.title.y=element_text(size=mareker_size),
                axis.title.x = element_text(size=mareker_size),
                plot.title = element_text(size=mareker_size),
                axis.text.x = element_text(size=tick_size),
                axis.text.y = element_text(size=tick_size)))
  if(save.fig)
    dev.off()

  ###
  # conditioning on other chromosone
  ###
  for(i in 1:chr){
    chr_index <- (m*(i-1)+1):(m*i)
    find_i <- setdiff(MixGeneObj$find, chr_index)
    MixGeneObj_i <- mixedModel(MixGeneObj, find = find_i, dupl = dupl, estPar = F)
    while(sum(abs(MixGeneObj_i$t[dupl[,1]]-MixGeneObj_i$t[dupl[,2]]))>0){
      MixGeneObj_i$t[dupl[,1]] <- MixGeneObj_i$t[dupl[,2]]
    }
    tvalues[chr_index]        <- MixGeneObj_i$t[chr_index]
  }
  if(save.fig){
    pdf(paste("Figure4_simple_",trait,".pdf",sep=''), w = 5, h = 3)
  }else{
    x11(w=5,h=3)
  }
  fig <- plotTrait(abs(tvalues), markers, MixGeneObj$crit)
  for(f in MixGeneObj$find){
    fig <- fig + geom_vline(xintercept = f,
                            color = "red",
                            size= col_v,
                            alpha = alpha_res)
  }
  print(fig +
          annotate("text", x =  text_annonate3[count], y = -0.2 , label = "QTL3") +
          annotate("text", x =  text_annonate1[count], y = -0.2 , label = "QTL1") +
          annotate("text", x =  text_annonate2[count], y = -0.2 , label = "QTL2") +
          theme(axis.title.y=element_text(size=mareker_size),
                axis.title.x = element_text(size=mareker_size),
                plot.title = element_text(size=mareker_size),
                axis.text.x = element_text(size=tick_size),
                axis.text.y = element_text(size=tick_size)))
  if(save.fig)
    dev.off()


}

