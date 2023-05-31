source("Methods.R")

###############
#Example Figure
###############

set.seed(1)
Exp=rep(10,1000)
s2a=0.5
true_s=inv_tri_gamma(s2a)
true_r=exp(digamma(inv_tri_gamma(s2a)))
Lambda=Exp*rgamma(1000,true_s,true_r)
Lambda[1:100]=Lambda[1:100]*exp(1)
Obs=rpois(1000,Lambda)
Z_FE=(Obs-Exp)/sqrt(Exp)

par(mfrow=c(1,2))

#Asymptotic EN Figure
empirical_null=asymptotic_en(Obs,Exp)
var_alpha=empirical_null$var_alpha
hist(Z_FE[101:1000],col="lightgray",freq=F,xlab="Null Z-Scores",main="(a) Conventional Asymptotic Method",xlim=c(-15,15),breaks=25,cex.lab=1.5,cex.main=1.5,ylim=c(0,0.23))
x=seq(-15,15,0.01)
lines(dnorm(x,median(Z_FE),sqrt(1+var_alpha*10))~x,lwd=3,col=2)
legend("topright",col=2,lwd=3,lty=1,legend="Estimated Null Distribution")

#Proposed EN Figure
empirical_null=exact_en(Obs,Exp)
s=empirical_null$s
r=empirical_null$r
xx=barplot(prop.table(table(Obs[101:1000])),axis.lty=1,col="lightgray",xlab="Null Observed Counts",ylab="Density",main="(b) Proposed Exact Method",cex.lab=1.5,cex.main=1.5)
lines(dnbinom(sort(unique(Obs[101:1000])),s,(r/10)/(1+r/10))~xx,lwd=3,col=2)
legend("topright",col=2,lwd=3,lty=1,legend="Estimated Null Distribution")

