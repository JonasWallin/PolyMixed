##
# creating the fist figure in
#  Ghost QTL and hotspots in experimental crosses - novel approach for modeling polygenic effects
##
library(PolyMixed)
save.fig = F
data('Figure1')

t2    <- 0:1500
X     <- Figure1$X
Y     <- Figure1$Y
gamma <- Figure1$gamma



par(mfrow=c(2,2))
plot(t2, type='l', xlab='Marker');
lim<-150*seq(1:9);
abline(v=lim, col='red', lty=2)
crit<-qnorm(1-0.05/(2*1500))
abline(h=c(crit,-crit), col='green',lwd=2)
her<-var(X%*%gamma)/var(Y);

per<-(apply(X,1,'mean')+1)/2;
c2<-cor(X,per);
plot(c2, type='l', xlab='Marker', ylab='Cor(X,U)')
abline(v=lim, col='red', lty=2)

save(t2,c2, file='intro.Rdata')


plot(t2~gamma, xlab=expression(gamma), ylab='t')
plot(t2~c2, xlab='Cor(X,U)', ylab='t')

if(save.fig){
  dev.print(file="Figure1",dev=pdf, width=15, height=15)
  dev.off()
}



