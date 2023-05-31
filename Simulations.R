source("Methods.R")
m=function(x) {return(mean(x,na.rm=T))}

#Estimation and Flagging for Count Data
simulation=function(effect,ntilde,nit=5000,I=100,s2a=0.5,dist="G") {
  true_s=inv_tri_gamma(s2a)
  true_r=exp(digamma(true_s))
  s2a_est1=vector("numeric")
  s2a_est2=vector("numeric")
  s=vector("numeric")
  r=vector("numeric")
  flag=matrix(0,nrow=nit,ncol=4)
  for(i in 1:nit) {
    E=runif(I,1,100)
    E[1]=ntilde
    gamma=rep(0,I)
    gamma[(I-0.1*I+1):(I-0.05*I)]=-1
    gamma[(I-0.05*I+1):I]=1
    gamma[1]=effect
    if(dist=="G") {
      lambda=rgamma(1000,true_s,true_r)
    } else if(dist=="LN") {
      lambda=rlnorm(1000,0,sqrt(s2a))
    } else if(dist=="W"){
      a=pi/sqrt(6)/sqrt(s2a)
      lambda=rweibull(1000,a,exp(-digamma(1)/a))
    }
    O=rpois(I,E*lambda*exp(gamma))
    flag_et=flag_naive_exact(O,E)[1]
    flag_asymp_en=flag_asymptotic_en(O,E)[1]
    flag_inv_en=flag_invert_en(O,E)[1]
    flag_et_en=flag_exact_en(O,E)[1]
    flag[i,]=c(flag_et,flag_asymp_en,flag_inv_en,flag_et_en)
    s2a_est1[i]=asymptotic_en(O,E)$var_alpha
    s2a_est2[i]=invert_en(O,E)$var_alpha
    ests=exact_en(O,E)
    s[i]=ests$s
    r[i]=ests$r
  }
  return(list(flag=flag,s=s,r=r,true_s=true_s,true_r=true_r,s2a_est1=s2a_est1,s2a_est2=s2a_est2,s2a=s2a))
}

#Sensitivity Analysis for Binary Data
simulation2=function(effect,nit=5000,I=100,s2a=0.5,norm=-1) {
  library(poibin)
  true_s=inv_tri_gamma(s2a)
  true_r=exp(digamma(inv_tri_gamma(s2a)))
  s2a_est1=vector("numeric")
  s2a_est2=vector("numeric")
  s=vector("numeric")
  r=vector("numeric")
  flag=matrix(0,nrow=nit,ncol=4)
  natl_props=vector("numeric")
  for(i in 1:nit) {
    N=100
    ct=rep(1:I,each=N)
    alpha=rep(log(rgamma(I,true_s,true_r)),each=N)
    Xmeans=rep(rnorm(I),each=N)
    X=rnorm(I*N,Xmeans)
    gamma=rep(0,I)
    gamma[(I-0.1*I+1):(I-0.05*I)]=-1
    gamma[(I-0.05*I+1):I]=1
    gamma[1]=effect
    gamma=rep(gamma,each=N)
    E=tapply(plogis(norm+X),ct,sum)
    m=plogis(norm+alpha+gamma+X)
    O=vector("numeric")
    for(c in 1:I) {
      O[c]=rpoibin(1,m[ct==c])
    }
    pr=data.frame(ct=ct,prs=plogis(norm+X))
    N=rep(N,I)
    flag_et=flag_naive_exact(O,E,N,pr)[1]
    flag_asymp_en=flag_asymptotic_en(O,E,N,pr)[1]
    flag_inv_en=flag_invert_en(O,E,N,pr)[1]
    flag_et_en=flag_exact_en(O,E,N,pr)[1]
    flag[i,]=c(flag_et,flag_asymp_en,flag_inv_en,flag_et_en)
    s2a_est1[i]=asymptotic_en(O,E,N,pr)$var_alpha
    s2a_est2[i]=invert_en(O,E,N,pr)$var_alpha
    ests=exact_en(O,E,N,pr)
    s[i]=ests$s
    r[i]=ests$r
    natl_props[i]=sum(O)/sum(N)
  }
  return(list(flag=flag,natl_props=mean(natl_props),res1=mean(s2a_est1),res2=mean(s2a_est2)))
}

combine=function(l1,l2) {
  n=length(l1)
  for(i in 1:n) {
    l1[[i]]=rbind(l1[[i]],l2[[i]])
  }
  return(l1)
}

#######################################################
#Estimation and False Flagging by Effective Sample Size
#######################################################
ntildes=seq(5,100,5)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=combine,.errorhandling="remove") %dopar% {
  simul=simulation(0,ntilde=ntildes[i])
  list(abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a)),
       100*abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a))/c(simul$true_s,simul$true_r,simul$s2a,simul$s2a),
       apply(simul$flag!=0,2,mean))
}
results=as.data.frame(result[[3]])
names(results)=c("Naive Exact","Asymptotic IEN","Inversion IEN","Exact IEN")
results$ntildes=ntildes
write.csv(results,"False_Flagging_Ntilde.csv",row.names=F)

results=as.data.frame(result[[2]])
names(results)=c("s % Bias","r % Bias","s2a % Bias Asymptotic IEN","s2a % Bias P-value IEN")
results$ntildes=ntildes
write.csv(results,"Estimation_Ntilde.csv",row.names=F)

#################################################################
#Estimation and False Flagging by Unobserved Confounding Variance
#################################################################
s2as=seq(0.05,1,0.05)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=combine,.errorhandling="remove") %dopar% {
  simul=simulation(0,10,s2a=s2as[i])
  list(abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a)),
       100*abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a))/c(simul$true_s,simul$true_r,simul$s2a,simul$s2a),
       apply(simul$flag!=0,2,mean))
}
results=as.data.frame(result[[3]])
names(results)=c("Naive Exact","Asymptotic IEN","Inversion IEN","Exact IEN")
results$s2as=s2as
write.csv(results,"False_Flagging_s2as.csv",row.names=F)

results=as.data.frame(result[[2]])
names(results)=c("s % Bias","r % Bias","s2a % Bias Asymptotic IEN","s2a % Bias P-value IEN")
results$s2as=s2as
write.csv(results,"Estimation_s2as.csv",row.names=F)

###############################################
#True Flagging (Lower Tail) by Effect Magnitude
###############################################
gammas=-seq(0.5,10,0.5)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=combine,.errorhandling="remove") %dopar% {
  simul=simulation(effect=gammas[i],10)
  list(abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a)),
       100*abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a))/c(simul$true_s,simul$true_r,simul$s2a,simul$s2a),
       apply(simul$flag < 0,2,mean))
}
results=as.data.frame(result[[3]])
names(results)=c("Naive Exact","Asymptotic IEN","Inversion IEN","Exact IEN")
results$gammas=gammas
write.csv(results,"True_Flagging_Low_Gammas.csv",row.names=F)


###############################################
#True Flagging (Upper Tail) by Effect Magnitude
###############################################
gammas=seq(0.05,2.5,length.out=20)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=combine,.errorhandling="remove") %dopar% {
  simul=simulation(effect=gammas[i],10)
  list(abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a)),
       100*abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a))/c(simul$true_s,simul$true_r,simul$s2a,simul$s2a),
       apply(simul$flag > 0,2,mean))
}
results=as.data.frame(result[[3]])
names(results)=c("Naive Exact","Asymptotic IEN","Inversion IEN","Exact IEN")
results$gammas=gammas
write.csv(results,"True_Flagging_Up_Gammas.csv",row.names=F)

####################################################
#True Flagging (Lower Tail) by Effective Sample Size
####################################################
ntildes=seq(5,100,5)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=combine,.errorhandling="remove") %dopar% {
  simul=simulation(-2,ntilde=ntildes[i])
  list(abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a)),
       100*abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a))/c(simul$true_s,simul$true_r,simul$s2a,simul$s2a),
       apply(simul$flag < 0,2,mean))
}
results=as.data.frame(result[[3]])
names(results)=c("Naive Exact","Asymptotic IEN","Inversion IEN","Exact IEN")
results$ntildes=ntildes
write.csv(results,"True_Flagging_Low_Ntilde.csv",row.names=F)

####################################################
#True Flagging (Upper Tail) by Effective Sample Size
####################################################
ntildes=seq(5,100,5)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=combine,.errorhandling="remove") %dopar% {
  simul=simulation(2,ntilde=ntildes[i])
  list(abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a)),
       100*abs(c(mean(simul$s)-simul$true_s,mean(simul$r)-simul$true_r,mean(simul$s2a_est1)-simul$s2a,mean(simul$s2a_est2)-simul$s2a))/c(simul$true_s,simul$true_r,simul$s2a,simul$s2a),
       apply(simul$flag > 0,2,mean))
}
results=as.data.frame(result[[3]])
names(results)=c("Naive Exact","Asymptotic IEN","Inversion IEN","Exact IEN")
results$ntildes=ntildes
write.csv(results,"True_Flagging_Up_Ntilde.csv",row.names=F)

#########################################
#Robustness (Gamma,Lognormal,Exponential)
#########################################
ntildes=seq(5,100,5)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=rbind,.errorhandling="remove") %dopar% {
  simul1=simulation(0,ntilde=ntildes[i],dist="G")
  simul2=simulation(0,ntilde=ntildes[i],dist="LN")
  simul3=simulation(0,ntilde=ntildes[i],dist="W")
  c(apply(simul1$flag!=0,2,mean)[4],apply(simul2$flag!=0,2,mean)[4],apply(simul3$flag!=0,2,mean)[4])
}
results=as.data.frame(result)
names(results)=c("Gamma","Lognormal","Weibull")
results$ntildes=ntildes
write.csv(results,"False_Flagging_Ntilde_Distribution.csv",row.names=F)

s2as=seq(0.05,1,0.05)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=rbind,.errorhandling="remove") %dopar% {
  simul1=simulation(0,10,s2a=s2as[i],dist="G")
  simul2=simulation(0,10,s2a=s2as[i],dist="LN")
  simul3=simulation(0,10,s2a=s2as[i],dist="W")
  c(apply(simul1$flag!=0,2,m)[4],apply(simul2$flag!=0,2,m)[4],apply(simul3$flag!=0,2,m)[4])
}
results=as.data.frame(result)
names(results)=c("Gamma","Lognormal","Weibull")
results$s2as=s2as
write.csv(results,"False_Flagging_s2as_Distribution.csv",row.names=F)

######################
#Robustness Binary FFP
######################
norms=seq(-6,6,0.5)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=rbind) %dopar% {
  simul=simulation2(0,norm=norms[i])
  c(apply(simul$flag!=0,2,mean),simul$natl_props,simul$res1,simul$res2)
}
results=as.data.frame(result)
names(results)=c("Naive Exact","Asymptotic IEN","Inversion IEN","Exact IEN","National Proportion")
write.csv(results,"False_Flagging_Binary.csv",row.names=F)


######################
#Robustness Binary TFP
######################
norms=seq(-6,6,0.5)
B=20 
cl=parallel::makeCluster(tail(numbers::divisors(B)[numbers::divisors(B)<=parallel::detectCores()],1))
doParallel::registerDoParallel(cl)
`%dopar%`=foreach::`%dopar%`
result=foreach::foreach(i=iterators::icount(B),.combine=rbind) %dopar% {
  simul=simulation2(-2,norm=norms[i])
  c(apply(simul$flag < 0,2,mean),simul$natl_props)
}
results=as.data.frame(result)
names(results)=c("Naive Exact","Asymptotic IEN","Inversion IEN","Exact IEN","National Proportion")
write.csv(results,"True_Flagging_Binary.csv",row.names=F)



