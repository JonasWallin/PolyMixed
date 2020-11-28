###
# Figure 3 in Ghost QTL and hotspots in experimental crosses - novel approach for modeling polygenic effects
# to produce Figure 3 a) set params_known=F
# to produce Figure 3 b) set params_known=T
# Testing the FWER for different tau's
# D: 2019-01-03
###
#rm(list=ls())
load.data = T
save.fig = F
save.file = F
plot_data = T
params_known = F
set.seed(1)
graphics.off()
library(PolyMixed)
library(bigstep)
library(RcppEigen)
library(ggplot2)
library(bigstep)
library(foreach)
library(doParallel)
#set.seed(10)
ns <- c(100,200,300,400,500)         # observations
taus <- c(0,0.01,0.05)     # sd for polygenic effect

# parameters----------------------------------------------------------------
q <- 200                # number of simulations
chr <- 10                # chromosomes
m <- 150                 # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
d <- 1                   # distance between markers (cM)
sigma <- 1               # sd for error term

mu <- rep(0.004, q)      # mean for polygenic effect

qtl.pos   <- c()
beta.true <- c()
result <- matrix(NA, nrow = length(ns), ncol = length(taus))
colnames(result) <- paste('tau_',taus,sep='')
rownames(result) <- paste('n_',ns,sep='')

if(load.data==F){
  find.mix.forward <- list()
  beta.mix.forward <- list()
  for(i in 1:length(taus)){
    find.mix.forward[[i]] <- list()
    beta.mix.forward[[i]] <- list()
  }
  n_count <- 0
  tau_in   = NULL
  sigma_in = NULL
  dupl_vec <- matrix(0, nrow=length(ns), ncol = q)

  for(n in ns){
    n_count <- n_count + 1
    pb <- txtProgressBar(min = 0, max = q, style = 3)
    for(i in 1:q) {
      # Should we allow for duplitacte?
      X <- simulateDesignMatrix(chr, m, n, d)
      SVDX = svd(X)
      dupl <- findDuplicate(X)
      dupl_vec[n_count,i] <- length(dupl[,1])
      count <- 0
      Z_beta<- rnorm(dim(X)[2])
      E     <- rnorm(dim(X)[1])
      for(tau in taus){
        count <- count + 1
        y <- X%*%(mu[1] + tau *Z_beta) + sigma * E
        if(params_known==F | tau>0){
        MixGeneObj <- PolyMixed::SetupGeneMix('Y ~ -1',
                                              data = data.frame(Y=y),
                                              X=X,
                                              SVDX = SVDX,
                                              meanMuOff = F,
                                              tauOff = F)
        }else{

          MixGeneObj <- PolyMixed::SetupGeneMix('Y ~ -1',
                                                data = data.frame(Y=y),
                                                X=X,
                                                SVDX = SVDX,
                                                meanMuOff = F,
                                                tauOff = T)
        }
        if(params_known==T){
          MixGeneObj$sigma <- sigma
          if(tau >0){
            MixGeneObj$tau   <- tau
            MixGeneObj$theta0 <- log(c(sigma, tau))
          }else{
            MixGeneObj$theta0 <- log(c(sigma))
          }

        }
        res <- PolyMixed::mixedModelForwardBackward(MixGeneObj,
                                         markers = markers,
                                         liberal = TRUE,
                                         dupl   = dupl,
                                         estParam = !params_known )

        if(is.null(res$find)){
          find.mix.forward[[count]][[i]]    <- list()
          beta.mix.forward[[count]][[i]]    <- list()

        }else{
          find.mix.forward[[count]][[i]]    <- res$find
          beta.mix.forward[[count]][[i]]    <- res$beta

        }
      }

      setTxtProgressBar(pb, i)
    }
    cat('\n')
    count <- 0
    for(tau in taus){

      count <- count + 1
      summ.forward <- summarizeFind.new(find.mix.forward[[count]], beta.mix.forward[[count]], beta.true, qtl.pos, maxCM = 15)
      #summ.back    <- summarizeFind.new(find.mix.backward, beta.mix.backward, beta.true, qtl.pos, maxCM = 15)

      result[n_count, count] <- summ.forward[3]

      cat('for (n = ',n,' tau  =',taus[count],'):\n')
      print(round(summ.forward, 3))
    }
  }
  if(save.file){
    if(params_known==F){
      save(result, file=paste("FWER_FIGURE_sim_",q,".RData",sep=""))
    }else{
      save(result, file=paste("FWER_FIGURE_noest_sim_",q,".RData",sep=""))
    }

  }
}else{
  if(params_known==F){
    load(file=paste("FWER_FIGURE_sim_",q,".RData",sep=""))
  }else{
    load(file=paste("FWER_FIGURE_noest_sim_",q,".RData",sep=""))
  }


}
if(plot_data){
  ns <- ns[-1]
  result<-result[-1,]
  n_n   <- length(ns)
  n_tau <- length(taus)
  plot_data <- c()
  p <- 0.05
  for(i in 1:n_tau){
    plot_data <- rbind(plot_data,
                       cbind(result[,i],
                             rep(taus[i],
                                 n_n),
                                 ns,
                             qbinom(p = p/2  , q, prob=0.05)/q,
                             qbinom(p = 1-p/2, q, prob=0.05)/q
                             ))
  }
  colnames(plot_data) <- c('FWER','tau','n','l','u')
  plot_data <- data.frame(plot_data)
  plot_data$tau <- factor(plot_data$tau)
  g<-ggplot(plot_data, aes(x=n, y=FWER, colour=tau)) + geom_line() +
    geom_hline(yintercept = qbinom(p = p/2  , q, prob=0.05)/q,
               width=25,
               color= 'black',
               linetype="dashed") +
    geom_hline(yintercept = qbinom(p = 1-p/2  , q, prob=0.05)/q,
               width=25,
               color= 'black',
               linetype="dashed") +
    geom_point() + ylim(0, 0.09)
  g<- g +   theme(axis.title.y=element_text(size=20),
                  axis.title.x = element_text(size=20),
                  plot.title = element_text(size=20),
                  legend.title.align=0.5) + labs(color=expression(tau))
  print(g)
  if(save.fig){
    if(params_known){
      ggsave('Figure3a.pdf', plot = g, width = 16, height = 16,  units = "cm")
    }else{
      ggsave('Figure3b.pdf', plot = g, width = 16, height = 16,  units = "cm")
    }
  }
}
