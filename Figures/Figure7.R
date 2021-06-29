###
#  real data analysis,
#  Creating figures for the real data analysis
#
# D: 2019-01-10
###
rm(list=ls())
graphics.off()
save.fig = T
library(RcppEigen)
library(ggplot2)
library(bigstep)
library(grid)
library(PolyMixed)
mareker_size = 30
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##
# setting up data
##
# data
data(zeng)
# X
X <- as.matrix(rbind(zeng$d1[, 1:45], zeng$d2[, 1:45]))
#X <- as.matrix(rbind(d1[, 1:45], d2[, 1:45], d3[, 1:45], d4[, 1:45]))
X <- X[, -(1:6)] # odrzucam 1. chromosom X
X[X == 9] <- NA
X <- replaceNA(X) - 1
n <- nrow(X)
M <- ncol(X)
# Y
y <- as.matrix(c(zeng$d1[, 46], zeng$d2[, 46]), ncol = 1)
#y <- as.matrix(c(d1[, 46], d2[, 46], d3[, 46], d4[, 46]), ncol = 1)
y <- scale(y)
# map
markers <- c(6, 16, 23)[-1] # bez 1. chromosomu
chr <- length(markers)
len <- c(0, 3.60, 10.60, 9.20, 17.20, 18.70, 0, 6.98, 10.10, 4.94, 6.51, 6.19,
         20.46, 12.78, 3.90, 4.55, 7.49, 30.02, 16.85, 4.34, 3.71, 7.03, 0, 4.99,
         9.34, 6.97, 7.44, 14.46, 6.79, 3.55, 6.32, 11.86, 4.58, 6.85, 6.35,
         11.79, 12.88, 9.15, 3.30, 7.98, 13.09, 10.04, 3.70, 9.79, 3.43)[-(1:6)]
len_cum <- unlist(tapply(len, rep(1:chr, markers), cumsum), use.names = FALSE)
by <- 2
res <- pseudomarkers(X, markers, len_cum, by)
X <- res$X
M <- ncol(X)
markers <- res$markers
q <- 1
findDuplicate2 <- function(X) {
  X[X < -0.7] <- -1
  X[X > 0.7] <- 1
  X[X > -0.3 & X < 0.3] <- 0
  dupl.ind <- which(duplicated(X, MARGIN = 2)) # ktore kolumny sie powtarzaja
  dupl.repl <- which(duplicated(X, MARGIN = 2, fromLast = TRUE)) # z tymi sie powtarzaja
  dupl <- cbind(dupl.ind, dupl.repl)
  return(dupl)
}
dupl <- findDuplicate2(X)


##
# estimating simple mode
##



unitv   <- as.matrix(rep(1,n))
t.simple <- apply(X, 2, function(x){fit <- fastLmPure(cbind(x,unitv), y); fit$coefficients[1] / fit$se[1]})


crit.simple <- calculateCrit(t.simple, markers)
find.simple <- which(abs(t.simple) > crit.simple)

if(save.fig==F)
  x11()
print(plotTrait(abs(t.simple), markers, crit.simple) +
        ylim(0, 12.5) +
        theme(axis.title.y=element_text(size=mareker_size),
              axis.title.x = element_text(size=mareker_size),
              plot.title = element_text(size=mareker_size)))

if(save.fig)
  ggsave('Figure7_simple.pdf')
###
# mixed model
###
MixGeneObj <- SetupGeneMix('Y ~ 1',
                           data = data.frame(Y=y),
                           X=X,
                           meanMuOff = F,
                           tauOff = F)
MixGeneObj <- mixedModelForwardBackward(MixGeneObj,
                                        markers,
                                        dupl = dupl)

t

if(save.fig==F)
  x11()

print(plotTrait(abs(MixGeneObj$t), markers, MixGeneObj$crit)+
        ylim(0, 12.5) +
        theme(axis.title.y=element_text(size=mareker_size),
              axis.title.x = element_text(size=mareker_size),
              plot.title = element_text(size=mareker_size))
      )

if(save.fig)
  ggsave('Figure7_mixed.pdf')
###
# CIMMmodel
###

t.CIM       <- CIMModel(X, y, window = 10)
crit.CIM    <- calculateCrit(t.CIM, markers)
find.CIM    <- which(abs(t.CIM) > crit.CIM)
if(save.fig==F)
  x11()


print(plotTrait(abs(t.CIM), markers, crit.CIM) +
        ylim(0, 12.5) +
        theme(axis.title.y=element_text(size=mareker_size),
              axis.title.x = element_text(size=mareker_size),
              plot.title = element_text(size=mareker_size)))
if(save.fig)
  ggsave('Figure7_CIM.pdf')
