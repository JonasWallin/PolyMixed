##
# creating the fist figure in
#  Ghost QTL and hotspots in experimental crosses - novel approach for modeling polygenic effects
##
library(PolyMixed)
library(latex2exp)
save.fig = T
data('Figure1')

t2    <- 0:1500
X     <- Figure1$X
Y     <- Figure1$Y
gamma <- Figure1$gamma

ttes<-function(x,y){
  obj1<-lm(y~x);
  u<-summary(obj1);
  return(u$coefficients[2,3])}
t2<-apply(X,2,ttes, y=Y);

#x11()
par(mar=c(5.1, 5.1, 4.1, 2.1))
par(mfrow=c(2,2))
plot(t2, type='l', xlab='Marker', ylab='t', cex.lab=2);
lim<-150*seq(1:9);
abline(v=lim, col='red', lty=2)
crit<-qnorm(1-0.05/(2*1500))
abline(h=c(crit,-crit), col='green',lwd=2)
her<-var(X%*%gamma)/var(Y);

per<-(apply(X,1,'mean')+1)/2;
c2<-cor(X,per);
y_lab <- TeX("Cor($X$,U)")
plot(c2, type='l', xlab='Marker', ylab=y_lab, cex.lab=2)
abline(v=lim, col='red', lty=2)

plot(t2~gamma, xlab=expression(paste(gamma)), ylab='t',cex.lab=2)
plot(t2~c2, xlab=y_lab, ylab='t',cex.lab=2)

if(save.fig){
  dev.print(file="Figure1.pdf",dev=pdf, width=15, height=15)
  dev.off()
}




