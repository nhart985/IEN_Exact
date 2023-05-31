source("Methods.R")
library(ggplot2)
library(ggpubr)
library(data.table)
library(readxl)

####################
#Dialysis 3-year SMR
####################

CMS=fread("DFR_Data_FY2022.csv")
CMS=as.data.frame(CMS)

E=as.numeric(CMS$fyexdz_f)
O=as.numeric(CMS$fydeaz_f)
O=O[!is.na(O)]
E=E[!is.na(E)]

#Flagging
flag_exact=-flag_naive_exact(O,E)
flag_asymp_en=-flag_asymptotic_en(O,E)
flag_inv_en=-flag_invert_en(O,E)
flag_et_en=-flag_exact_en(O,E)
table(flag_exact)
100*prop.table(table(flag_exact))
table(flag_asymp_en)
100*prop.table(table(flag_asymp_en))
table(flag_inv_en)
100*prop.table(table(flag_inv_en))
table(flag_et_en)
100*prop.table(table(flag_et_en))
table(flag_et_en,flag_exact)
table(flag_et_en,flag_asymp_en)
table(flag_et_en,flag_inv_en)

#Histogram for SMR
Z=(O-E)/sqrt(E)
low=1
up=2
approx=(low+up)/2
en=asymptotic_en(O,E)
A=-1.96*sqrt(1+en$var_alpha*E)
B=1.96*sqrt(1+en$var_alpha*E)
par(mfrow=c(1,2))
hist(Z[(E > low & E < up) & (A < Z & Z < B)],col="lightgray",freq=F,breaks=10,main="(a) Conventional Asymptotic Method",ylim=c(0,0.46),xlab="Naive Z-Scores from Null Intervals",ylab="Density",cex.axis=1.25,cex.main=1.5,cex.lab=1.25)
x=seq(-1.5,2,0.01)
lines(x,dnorm(x,mean=0,sd=sqrt(1+approx*en$var_alpha)),lwd=3,col=2)
legend("topright",col=2,lwd=3,lty=1,legend="Estimated Null Distribution")
en=exact_en(O,E)
xx=barplot(prop.table(table(O[(E > low & E < up) & (A < Z & Z < B)])),col="lightgray",axis.lty=1,ylim=c(0,0.46),main="(b) Proposed Exact Method",xlab="Observed Counts from Null Intervals",ylab="Density",cex.axis=1.25,cex.main=1.5,cex.lab=1.25)
x=seq(0,4,1)
lines(dnbinom(x,en$s,(en$r/approx)/(1+en$r/approx))~xx,lwd=3,col=2,type="o")
legend("topright",col=2,lwd=3,lty=1,legend="Estimated Null Distribution")

#Funnel Plot
funnel_pois=function(O,E,N=NULL) {
  Exp=seq(min(E),max(E),length.out=1000)
  #Naive Exact 
  temp_upper=(1-ppois(qpois(0.975,Exp)-1,Exp)-0.025)/dpois(qpois(0.975,Exp),Exp)
  temp_lower=(ppois(qpois(0.025,Exp),Exp)-0.025)/dpois(qpois(0.025,Exp),Exp)
  upper_naive=(qpois(0.975,Exp)+temp_upper-0.5)/Exp
  lower_naive=(qpois(0.025,Exp)-temp_lower+0.5)/Exp
  #Asymptotic
  en=asymptotic_en(O,E,N)
  upper_asymptotic=1+qnorm(0.975)*sqrt((1+en$var_alpha*Exp)/Exp)
  lower_asymptotic=1-qnorm(0.975)*sqrt((1+en$var_alpha*Exp)/Exp)
  #Inversion
  en=invert_en(O,E,N)
  o_u=qpois(pmin(pnorm(qnorm(0.975)*sqrt(1+en$var_alpha*Exp)),0.999999999),Exp)
  o_l=qpois(pnorm(-qnorm(0.975)*sqrt(1+en$var_alpha*Exp)),Exp)
  temp_upper=(1-ppois(o_u-1,Exp)-pnorm(-1.96*sqrt(1+en$var_alpha*Exp)))/dpois(o_u,Exp)
  temp_lower=(ppois(o_l,Exp)-pnorm(-1.96*sqrt(1+en$var_alpha*Exp)))/dpois(o_l,Exp)
  upper_invert=(o_u+temp_upper-0.5)/Exp
  lower_invert=(o_l-temp_lower+0.5)/Exp
  #Exact IEN NB
  en=exact_en(O,E,N)
  shape=en$s
  rate=(en$r/Exp)/(1+(en$r/Exp))
  temp_upper=(1-pnbinom(qnbinom(0.975,shape,rate)-1,shape,rate)-0.025)/dnbinom(qnbinom(0.975,shape,rate),shape,rate)
  temp_lower=(pnbinom(qnbinom(0.025,shape,rate),shape,rate)-0.025)/dnbinom(qnbinom(0.025,shape,rate),shape,rate)
  upper_exact=(qnbinom(0.975,shape,rate)+temp_upper-0.5)/Exp
  lower_exact=(qnbinom(0.025,shape,rate)-temp_lower+0.5)/Exp
  df=data.frame(method=rep(c("Naive","Asymptotic","Inversion","Exact"),each=length(Exp)),
                upper=c(upper_naive,upper_asymptotic,upper_invert,upper_exact),
                lower=c(lower_naive,lower_asymptotic,lower_invert,lower_exact),
                exp=Exp)
  return(df)
}

get_data=function(cutoff){
  plot_data=funnel_pois(O,E)
  plot_data=plot_data[plot_data$exp < cutoff,]
  plot_data$method=factor(plot_data$method,levels=c("Naive","Asymptotic","Inversion","Exact"))
  plot_data$method=factor(plot_data$method,labels=c("Naive Exact","Asymptotic Normal IEN",
                                                    "P-value Inversion IEN", "Exact Poisson IEN"))
  plot_data2=data.frame(E=E, SMR=O/E)
  plot_data2=plot_data2[plot_data2$E < cutoff,]
  plot_data2=plot_data2[plot_data2$SMR < 7,]
  return(list(plot_data,plot_data2))
}

pl=get_data(Inf)
plot_data=pl[[1]]; plot_data2=pl[[2]]

g=ggplot()
g=g+geom_point(data=plot_data2,mapping=aes(x=E,y=SMR),colour="gray")
g=g+geom_line(data=plot_data,mapping=aes(x=exp,y=upper,colour=method,linetype=method),size=1)+geom_line(data=plot_data,mapping=aes(x=exp,y=lower,colour=method,linetype=method),size=1)
g=g+ylab("SMR")+xlab("Effective Sample Size")
g=g+theme_classic()
g=g+theme(legend.position="top",legend.title=element_blank())
g=g+scale_colour_manual(values=c("Naive Exact"="black","Asymptotic Normal IEN"="blue",
                                 "P-value Inversion IEN"="dark green","Exact Poisson IEN"="dark orange"))
g=g+scale_linetype_manual(values=c("Naive Exact"="solid","Asymptotic Normal IEN"="dotted",
                                   "P-value Inversion IEN"="dashed","Exact Poisson IEN"="dotdash"))
g=g+scale_y_continuous(breaks=c(0,3,6),labels=c(0,3,6))
g=g+theme(text=element_text(size=20))
g=g+ggtitle("(a)")
g=g+geom_hline(yintercept=1,linetype="dashed",size=1)
g=g+theme(legend.key.width=unit(1.25,"cm"))
g=g+guides(colour=guide_legend(ncol=2))

pl=get_data(15)
plot_data=pl[[1]]; plot_data2=pl[[2]]

library(ggplot2)
g2=ggplot()
g2=g2+geom_point(data=plot_data2,mapping=aes(x=E,y=SMR),colour="gray")
g2=g2+geom_line(data=plot_data,mapping=aes(x=exp,y=upper,colour=method,linetype=method),size=1)+geom_line(data=plot_data,mapping=aes(x=exp,y=lower,colour=method,linetype=method),size=1)
g2=g2+ylab("SMR")+xlab("Effective Sample Size")
g2=g2+theme_classic()
g2=g2+theme(legend.position="top",legend.title=element_blank())
g2=g2+scale_colour_manual(values=c("Naive Exact"="black","Asymptotic Normal IEN"="blue",
                                   "P-value Inversion IEN"="dark green","Exact Poisson IEN"="dark orange"))
g2=g2+scale_linetype_manual(values=c("Naive Exact"="solid","Asymptotic Normal IEN"="dotted",
                                     "P-value Inversion IEN"="dashed","Exact Poisson IEN"="dotdash"))
g2=g2+theme(text=element_text(size=20))
g2=g2+ggtitle("(b)")
g2=g2+scale_y_continuous(breaks=c(0,3,6,9),labels=c(0,3,6,9))
g2=g2+geom_hline(yintercept=1,linetype="dashed",size=1)
g2=g2+theme(legend.key.width=unit(1.25,"cm"))
g2=g2+guides(colour=guide_legend(ncol=2))

ggarrange(g,g2,common.legend=T)


#######################
#Transplant Acceptance
#######################

CTR=as.data.frame(read_excel("csrs_final_tables_2111_KI.xls",
                             sheet="Table B11 & Figures B10-B14"))
N=as.numeric(CTR$OA_OVERALL_OFFERS_CENTER)
O=as.numeric(CTR$OA_OVERALL_ACCEPTS_CENTER)
E=as.numeric(CTR$OA_OVERALL_EXP_ACCEPTS_CENTER)

#Flagging
N=N[!is.na(N)]
O=O[!is.na(O)]
E=E[!is.na(E)]
flag_et_pb=flag_naive_exact(O,E,N)
flag_asymp_en=flag_asymptotic_en(O,E)
flag_inv_en=flag_invert_en(O,E)
flag_et_en=flag_exact_en(O,E,N)
table(flag_et_pb)
100*prop.table(table(flag_et_pb))
table(flag_asymp_en)
100*prop.table(table(flag_asymp_en))
table(flag_inv_en)
100*prop.table(table(flag_inv_en))
table(flag_et_en)
100*prop.table(table(flag_et_en))
table(flag_et_en,flag_et_pb)
table(flag_et_en,flag_asymp_en)
table(flag_et_en,flag_inv_en)

pl=get_data(Inf)
plot_data=pl[[1]]; plot_data2=pl[[2]]

g=ggplot()
g=g+geom_point(data=plot_data2,mapping=aes(x=E,y=SMR),colour="gray")
g=g+geom_line(data=plot_data,mapping=aes(x=exp,y=upper,colour=method,linetype=method),size=1)+geom_line(data=plot_data,mapping=aes(x=exp,y=lower,colour=method,linetype=method),size=1)
g=g+geom_hline(yintercept=1,linetype="longdash",size=1,colour="gray48")+ylim(-0.5,3.5)
g=g+ylab("SAR")+xlab("Effective Sample Size")
g=g+theme_classic()
g=g+theme(legend.position="top",legend.title=element_blank())
g=g+scale_colour_manual(values=c("Naive Exact"="black","Asymptotic Normal IEN"="blue",
                                 "P-value Inversion IEN"="dark green","Exact Poisson IEN"="dark orange"))
g=g+scale_linetype_manual(values=c("Naive Exact"="solid","Asymptotic Normal IEN"="dotted",
                                   "P-value Inversion IEN"="dashed","Exact Poisson IEN"="dotdash"))
g=g+theme(text=element_text(size=20),legend.key.width=unit(1.25,"cm"))
g=g+guides(colour=guide_legend(ncol=2))
g







