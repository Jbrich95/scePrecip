#Negative log-Likelihood for free fits - single conditioning site 

#Delta-Laplace margins
free_nll_DL=function(X,Y,par,c,v=NULL){
  
  alpha=par[1]
  beta=par[2]
  sig=par[3]
  mu=par[4]
  delta=par[5]
  
  X1=X[Y>c]
  X2=X[Y <= c]
  Y=Y[Y>c]

  #alpha and beta Positive
   if(alpha < 0 || alpha > 1  || beta > 1 || beta < 0 || delta < 0){return(1e10)}
  
  if(sig<=0){return(1e10)}
  
  
  negloglik <- -sum(ddlaplace(Y,alpha*X1+mu*X1^beta,sig*X1^beta,delta,log=T))
  negloglik <- negloglik-sum(
    log(
      apply(cbind(alpha*X2+mu*X2^beta,sig*X2^beta),1,function(x){
        pdlaplace(c,x[1],x[2],delta)
      })
    )
  )
  
  
  
  
  if(is.finite(negloglik)){
    return(negloglik)
  }else{return(1e10)}
}

#Gaussian margins
free_nll_Gauss=function(X,Y,par,c,v=NULL){
  
  alpha=par[1]
  beta=par[2]
  sig=par[3]
  mu=par[4]

  
  X1=X[Y>c]
  X2=X[Y <= c]
  Y=Y[Y>c]
  
  #alpha and beta Positive
  if(alpha < 0 || alpha > 1  || beta > 1 || beta < 0){return(1e10)}
  
  if(sig<=0){return(1e10)}
  
  negloglik <- -sum(dnorm(Y,alpha*X1+mu*((X1)^beta),sig*((X1)^beta),log=T))-sum(pnorm(c,alpha*X2+mu*((X2)^beta),sig*((X2)^beta),log=T))

  
  if(is.finite(negloglik)){
    return(negloglik)
  }else{return(1e10)}
}

#Gaussian margins, fixed beta
free_nll_Gauss_fixedbeta=function(X,Y,par,beta,c,v=NULL){
  
  alpha=par[1]
  sig=par[2]
  mu=par[3]
  
  
  X1=X[Y>c]
  X2=X[Y <= c]
  Y=Y[Y>c]
  
  #alpha and beta Positive
  if(alpha < 0 || alpha > 1  || beta > 1 || beta < 0){return(1e10)}
  
  if(sig<=0){return(1e10)}
  
  negloglik <- -sum(dnorm(Y,alpha*X1+mu*((X1)^beta),sig*((X1)^beta),log=T))-sum(pnorm(c,alpha*X2+mu*((X2)^beta),sig*((X2)^beta),log=T))
  
  
  if(is.finite(negloglik)){
    return(negloglik)
  }else{return(1e10)}
}
