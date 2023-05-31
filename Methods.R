
#Inverse Tri-Gamma Function
inv_tri_gamma=function(var_alpha) {
  fn=function(x,var_alpha) {
    return(var_alpha-trigamma(x))
  }
  return(uniroot(fn,c(0,1e10),var_alpha=var_alpha)$root)
}

#Mid P-Value Poisson 
mid_pval_P=function(O,E) {
  temp=0.5*(ppois(O,E)+ppois(O-1,E))
  pval=2*pmin(temp,1-temp)
  return(pval)
}

#Mid P-Value Poisson-Binomial 
mid_pval_PB=function(O,E,N,pr=NULL) {
  library(poibin)
  if(is.null(pr)) {
    pval=mid_pval_P(O,E)
  } else {
    temp=vector("numeric")
    for(i in 1:length(O)) {
      temp[i]=0.5*(ppoibin(O[i],pr$prs[pr$ct==i])+ppoibin(O[i]-1,pr$prs[pr$ct==i]))
    }
    pval=2*pmin(temp,1-temp)
  }
  return(pval)
}

#Mid P-Value Negative-Binomial
mid_pval_NB=function(O,E,s,r,N=NULL,pr=NULL) {
  if(is.null(N)) {
    ntilde=E
  } else {
    if(is.null(pr)) {
      ntilde=E*(1-E/N)
    } else {
      f=function(x) {return(sum(x*(1-x)))}
      ntilde=sapply(split(pr$prs,pr$ct),f)
    }
  }
  temp=0.5*(pnbinom(O,s,(r/ntilde)/(1+r/ntilde))+pnbinom(O-1,s,(r/ntilde)/(1+(r/ntilde))))
  pval=2*pmin(temp,1-temp)
  return(pval)
}

#Asymptotic Normal EN Estimation
asymptotic_en=function(O,E,N=NULL,pr=NULL) {
  library(MASS)
  if(is.null(N)) {
    ntilde=E
  } else {
    if(is.null(pr)) {
      ntilde=E*(1-E/N)
    } else {
      f=function(x) {return(sum(x*(1-x)))}
      ntilde=sapply(split(pr$prs,pr$ct),f)
    }
  }
  Z_FE=(O-E)/sqrt(ntilde)
  var_alpha0=pmax((rlm(Z_FE~1)$s^2-1)/mean(ntilde),0.0001)
  A=-1.96*sqrt(1+var_alpha0*ntilde)
  B=1.96*sqrt(1+var_alpha0*ntilde)
  negloglik=function(par,Z_FE,O,ntilde,A,B,p_null) {
    indices_null=which(Z_FE > A & Z_FE < B)
    I_null=length(indices_null)
    indices_non_null=which(Z_FE < A | Z_FE > B)
    Z_FE_null=Z_FE[indices_null]
    ntilde_null=ntilde[indices_null]
    ntilde_non_null=ntilde[indices_non_null]
    A_non_null=A[indices_non_null]
    B_non_null=B[indices_non_null]
    Q=pnorm(B_non_null,0,sqrt(1+par*ntilde_non_null))-pnorm(A_non_null,0,sqrt(1+par*ntilde_non_null))
    loglik=I_null*log(p_null)+sum(log(dnorm(Z_FE_null,0,sqrt(1+par*ntilde_null))))+sum(log(1-p_null*Q))
    return(-loglik)
  }
  p_null_grid=seq(0.50,0.99,0.01)
  var_alpha_hat=vector("numeric")
  negloglik_vals=vector("numeric")
  for(i in 1:length(p_null_grid)) {
    est=optim(par=var_alpha0,fn=negloglik,Z_FE=Z_FE,O=O,ntilde=ntilde,A=A,B=B,p_null=p_null_grid[i],
              method="Brent",lower=0,upper=100)
    var_alpha_hat[i]=est$par
    negloglik_vals[i]=est$value
  }
  i_best=which.min(negloglik_vals)
  p_null_best=p_null_grid[i_best]
  var_alpha_best=var_alpha_hat[i_best]
  return(list(var_alpha=var_alpha_best,p_null=p_null_best))
}

#Invert P-Value EN
invert_en=function(O,E,N=NULL,pr=NULL) {
  library(MASS)
  if(is.null(N)) {
    ntilde=E
    pvals=mid_pval_P(O,E)/2
  } else {
    if(is.null(pr)) {
      ntilde=E*(1-E/N)
    } else {
      f=function(x) {return(sum(x*(1-x)))}
      ntilde=sapply(split(pr$prs,pr$ct),f)
    }
    pvals=mid_pval_PB(O,E,N,pr)/2
  }
  Z_FE=rep(0,length(pvals))
  Z_FE[!pvals %in% c(0,1)]=ifelse(O[!pvals %in% c(0,1)] < E[!pvals %in% c(0,1)], qnorm(pvals[!pvals %in% c(0,1)]), qnorm(1-pvals[!pvals %in% c(0,1)]))
  Z_FE[pvals==0 & O < E]=-38
  Z_FE[pvals==0 & O > E]=38
  var_alpha0=pmax((rlm(Z_FE~1)$s^2-1)/mean(ntilde),0.0001)
  A=-1.96*sqrt(1+var_alpha0*ntilde)
  B=1.96*sqrt(1+var_alpha0*ntilde)
  negloglik=function(par,Z_FE,O,ntilde=ntilde,A,B,p_null) {
    indices_null=which(Z_FE > A & Z_FE < B)
    I_null=length(indices_null)
    indices_non_null=which(Z_FE < A | Z_FE > B)
    Z_FE_null=Z_FE[indices_null]
    ntilde_null=ntilde[indices_null]
    ntilde_non_null=ntilde[indices_non_null]
    A_non_null=A[indices_non_null]
    B_non_null=B[indices_non_null]
    Q=pnorm(B_non_null,0,sqrt(1+par*ntilde_non_null))-pnorm(A_non_null,0,sqrt(1+par*ntilde_non_null))
    loglik=I_null*log(p_null)+sum(log(dnorm(Z_FE_null,0,sqrt(1+par*ntilde_null))))+sum(log(1-p_null*Q))
    return(-loglik)
  }
  p_null_grid=seq(0.50,0.99,0.01)
  var_alpha_hat=vector("numeric")
  negloglik_vals=vector("numeric")
  for(i in 1:length(p_null_grid)) {
    est=optim(par=var_alpha0,fn=negloglik,Z_FE=Z_FE,O=O,ntilde=ntilde,A=A,B=B,p_null=p_null_grid[i],
              method="Brent",lower=0,upper=100)
    var_alpha_hat[i]=est$par
    negloglik_vals[i]=est$value
  }
  i_best=which.min(negloglik_vals)
  p_null_best=p_null_grid[i_best]
  var_alpha_best=var_alpha_hat[i_best]
  return(list(var_alpha=var_alpha_best,p_null=p_null_best))
}

#Starting Values for the Exact EN
get_starts=function(O,E,N=NULL,pr=NULL) {
  if(is.null(N)) {
    var_alpha0=asymptotic_en(O,E)$var_alpha
    s0=inv_tri_gamma(var_alpha0)
    r0=exp(digamma(s0))
  } else {
    var_alpha0=asymptotic_en(O,E,N,pr)$var_alpha
    s0=1/var_alpha0
    r0=1/var_alpha0
  }
  return(list(s0=s0,r0=r0))
}

#Exact EN Estimation
exact_en=function(O,E,N=NULL,pr=NULL) {
  if(is.null(N)) {
    ntilde=E
  } else {
    if(is.null(pr)) {
      ntilde=E*(1-E/N)
    } else {
      f=function(x) {return(sum(x*(1-x)))}
      ntilde=sapply(split(pr$prs,pr$ct),f)
    }
  }
  starts=get_starts(O,E,N,pr)
  s0=starts$s0
  r0=starts$r0
  A=qnbinom(0.025,s0,(r0/ntilde)/((r0/ntilde)+1))
  B=qnbinom(0.975,s0,(r0/ntilde)/((r0/ntilde)+1))
  negloglik=function(par,O,ntilde,A,B,p_null) {
    s=par[1]
    r=par[2]
    indices_null=which(O >= A & O <= B)
    I_null=length(indices_null)
    indices_non_null=which(O < A | O > B)
    O_null=O[indices_null]
    ntilde_null=ntilde[indices_null]
    ntilde_non_null=ntilde[indices_non_null]
    A_non_null=A[indices_non_null]
    B_non_null=B[indices_non_null]
    Q=pnbinom(B_non_null,s,(r/ntilde_non_null)/((r/ntilde_non_null)+1))-pnbinom(A_non_null-1,s,(r/ntilde_non_null)/((r/ntilde_non_null)+1))
    loglik=I_null*log(p_null)+sum(log(dnbinom(O_null,s,(r/ntilde_null)/((r/ntilde_null)+1))))+sum(log(1-p_null*Q))
    return(-loglik)
  }
  p_null_grid=seq(0.50,0.99,0.01)
  s_hat=vector("numeric")
  r_hat=vector("numeric")
  negloglik_vals=vector("numeric")
  for(i in 1:length(p_null_grid)) {
    est=optim(par=c(s0,r0),fn=negloglik,O=O,ntilde=ntilde,A=A,B=B,p_null=p_null_grid[i])
    param=est$par
    s_hat[i]=param[1]
    r_hat[i]=param[2]
    negloglik_vals[i]=est$value
  }
  i_best=which.min(negloglik_vals)
  p_null_best=p_null_grid[i_best]
  s_best=s_hat[i_best]
  r_best=r_hat[i_best]
  return(list(s=s_best,r=r_best,p_null=p_null_best))
}

#Flag Naive Exact Test
flag_naive_exact=function(O,E,N=NULL,pr=NULL) {
  if(is.null(N)) {
    pvals=mid_pval_P(O,E)
  } else {
    pvals=mid_pval_PB(O,E,N,pr)
  }
  flag=rep(0,length(O))
  flag[pvals < 0.05 & O > E]=1
  flag[pvals < 0.05 & O < E]=-1
  return(flag)
}

#Flag Asymptotic Normal EN Test
flag_asymptotic_en=function(O,E,N=NULL,pr=NULL) {
  if(is.null(N)) {
    ntilde=E
  } else {
    if(is.null(pr)) {
      ntilde=E*(1-E/N)
    } else {
      f=function(x) {return(sum(x*(1-x)))}
      ntilde=sapply(split(pr$prs,pr$ct),f)
    }
  }
  en=asymptotic_en(O,E,N,pr)
  Z_FE=(O-E)/sqrt(ntilde)
  Z_EN=Z_FE/sqrt(1+en$var_alpha*ntilde)
  flag=rep(0,length(Z_EN))
  flag[Z_EN > 1.96]=1
  flag[Z_EN < -1.96]=-1
  return(flag)
}

#Flag Poisson Invert P-Value Test
flag_invert_en=function(O,E,N=NULL,pr=NULL) {
  if(is.null(N)) {
    ntilde=E
    pvals=mid_pval_P(O,E)/2
  } else {
    if(is.null(pr)) {
      ntilde=E*(1-E/N)
    } else {
      f=function(x) {return(sum(x*(1-x)))}
      ntilde=sapply(split(pr$prs,pr$ct),f)
    }
    pvals=mid_pval_PB(O,E,N,pr)/2
  }
  en=invert_en(O,E,N,pr)
  Z_FE=rep(0,length(pvals))
  Z_FE[!pvals %in% c(0,1)]=ifelse(O[!pvals %in% c(0,1)] < E[!pvals %in% c(0,1)], qnorm(pvals[!pvals %in% c(0,1)]), qnorm(1-pvals[!pvals %in% c(0,1)]))
  Z_FE[pvals==0 & O < E]=-38
  Z_FE[pvals==0 & O > E]=38
  Z_EN=Z_FE/sqrt(1+en$var_alpha*ntilde)
  flag=rep(0,length(Z_EN))
  flag[Z_EN > 1.96]=1
  flag[Z_EN < -1.96]=-1
  return(flag)
}

#Flag Exact EN Test
flag_exact_en=function(O,E,N=NULL,pr=NULL) {
  en=exact_en(O,E,N,pr)
  pvals=mid_pval_NB(O,E,en$s,en$r,N,pr)
  flag=rep(0,length(pvals))
  flag[pvals < 0.05 & O > E]=1
  flag[pvals < 0.05 & O < E]=-1
  return(flag)
}
