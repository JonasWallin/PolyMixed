##
# File for creating the second figure in
# Ghost QTL and hotspots in experimental crosses - novel approach for modeling polygenic effects
# Creating the convergence plots for the expontial correlation
#
# D: 2019-05-04
##
save.fig = F

set.seed(3)
# parameters----------------------------------------------------------------
q <- 1                   # number of simulations
chr <- 10                # chromosomes
m <- 150                 # markers on chromosome
M <- chr * m             # all markers
markers <- rep(m, chr)   # markers on each chromosome
n <- 400                 # observations
d <- 1                   # distance between markers (cM)
sigma <- 1               # sd for error term
tau <- 0.0# sd for polygenic effect
mu <- rep(0.004, q)      # mean for polygenic effect
qtl.pos <- c(151, 181, 1275)
signal <- 0.5
qtl  <- c(signal, qtl.pos[1], 1) # qtl (beta; position: marker, trait or traits)
qtl2 <- c(signal, qtl.pos[2], 1)
qtl3 <- c(-0, qtl.pos[3], 1)
# --------------------------------------------------------------------------

#source("../codes/misc/3-exp_matrix.R")
library(RcppEigen)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(PolyMixed)
#source("../codes/misc/4-help_functions.R")



#i = 1
Xm <- PolyMixed::simulateDesignMatrix(chr, m, n, d)
count <- 0
plist <- list()
taus <- c(0.01,0.2)
ns   <- c(100, 400)
for(tau in taus){
  for(n in ns){
    count <- count + 1
    X <- Xm[1:n,]
    y <- X%*%rnorm(dim(X)[2], sd = tau) + rnorm(dim(X)[1], sd = sigma)
    dupl <- PolyMixed::findDuplicate(X)


    SVDX <- svd(X)
    MixGeneObj <- PolyMixed::SetupGeneMix('Y ~ -1',
                                         data = data.frame(Y=y),
                                         X=X,
                                         SVDX = SVDX,
                                         meanMuOff = T,
                                         tauOff = F)
    k <- 5
    MixGeneObj$theta0 <- log(c(sigma, tau))
    MixGeneObj <- PolyMixed::mixedModel(MixGeneObj, find = NULL, dupl = dupl, estPar = F)
    M <- sum(markers) - length(markers)
    cor.t <- numeric(k)
    cumm <- cumsum(markers)
    t <- MixGeneObj$t
    for (i in 1:k) {
      t1 <- t[-cumm][1:(M-i)]
      t2 <- t[-(cumm + i)][(i+1):M]
      cor.t[i] <- cor(t1, t2)
    }
    x    <-  d*(1:max(sum(cor.t > 0.3), 2))
    cor.tx <- cor.t[x]
    beta <- -coef(lm(log(cor.tx) ~ x - 1))
    p <- ggplot(data = data.frame(x = 1:length(cor.t), y = log(cor.t)), aes(x, y))
    p <- p + geom_abline(intercept = 0, slope = -beta,colour='blue') + geom_point()
    p <- p + labs(x='Distance in cM', y='log(Cor)')
    p <- p + theme(axis.title.y=element_text(size=20), axis.title.x = element_text(size=20))
    p <- p + ggtitle(bquote(tau~'='~.(tau)~","~n~'='~.(n)))
    p <- p + theme(plot.title = element_text(hjust = 0.5,size=20))
    plist[[count]] <- p
  }
}
if(save.fig)
  ggpubr::ggexport(marrangeGrob(plist, nrow=2, ncol=2, top=NULL), filename ="OU.pdf")
#ggsave("OU.pdf",marrangeGrob(plist, nrow=2, ncol=2))
