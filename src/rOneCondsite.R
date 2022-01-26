source("src/DL_funcs.R")

#Simulate one field from fitted model for {Y(s)}|(Y(s_O)>v)

##Inputs:
# 
# For n observed fields with d sampling locations
#
#  par: parameter estimates for AI spatial model, see spatial_fit.R
#  c.vec: d-vector of censoring thresholds on Laplace scale
#  coords: d x 2 matrix of lon/lat coordinates
#  v: exceedance threshold, must satissty v >= u for u defined in spatial_fit.R
#  Cond.ind: index of row in coords that gives conditioning site
#  n.sims: number of replications to produce


rOneCondsite=function(par,v,c.vec,coords,Cond.ind,n.sims){
  
  p=dim(coords)[1]
  
  KA1=par[1]
  KA2=par[2]
  KB1=par[3]
  KB2=par[4]
  KB3=1
  
  KM1=par[5]
  KM2=par[6]
  KM3=par[7]
  KS1=par[8]
  KS2=par[9]
  KR1=par[10]
  KR2=par[11]
  KD1=par[12]
  KD2=par[13]
  KD3=par[14]
  KD4=1
  
  re.coord=coords[c(Cond.ind,(1:p)[-Cond.ind]),]
  re.c=c.vec[c(Cond.ind,(1:p)[-Cond.ind])]
  
  h.pairs=rdist.earth(re.coord,miles=F)
  diag(h.pairs)=rep(0,p)
  Cor=Matern(h.pairs,range=KR1,smoothness = KR2)
  h.Cond=h.pairs[1,]
  rm(h.pairs)
  
  alpha=exp(-(abs(h.Cond)/KA1)^KA2)
  alpha=alpha[-1]
  
  beta=KB3*exp(-(abs(h.Cond)/KB1)^KB2)
  beta=beta[-1]
  
  sig=sqrt(2)*(1-exp(-(h.Cond/KS1)^KS2))
  sig=sig[-1]
  mu=KM1*h.Cond^{KM2}*exp(-h.Cond/KM3)
  mu=mu[-1]
  delta=1+(KD1*h.Cond^{KD2}-KD4)*exp(-h.Cond/KD3)
  delta=delta[-1]
  
  Cor=(Cor[-1,-1]-Cor[1,-1]%*%t(Cor[1,-1]))/((1-Cor[1,-1]^2)^(1/2)%*%(1-t(Cor[-1,1])^2)^(1/2))
  if(p > 2){
    Sigma=diag(sig)%*%Cor%*%diag(sig)
  }else{
    Sigma=sig^2*Cor
  }
  X0=v+rexp(n.sims)
  Z=rmvdlaplace(n=n.sims,dim=length(sig),mu=mu,sigmad=sqrt(diag(Sigma)),Sigma=Sigma,delta=delta)
  
  X=Z
  for(i in 1:n.sims){
    X[i,]=X0[i]*alpha+(X0[i]^beta)*Z[i,]
    for(j in 1:(p-1)){
      
      if(X[i,j]<re.c[j+1]){X[i,j]= -Inf}
      
    }
  }
  
  if(n.sims==1){
    if(Cond.ind!=1 & Cond.ind !=p){
      return(c(X[1:(Cond.ind-1)],X0,X[(Cond.ind):(p-1)]))
    }else if(Cond.ind==1){
      return(c(X0,X[(Cond.ind):(p-1)]))
      
    }else if(Cond.ind==p){
      return(c(X[1:(Cond.ind-1)],X0))
      
    }
  }else{
    if(Cond.ind!=1 & Cond.ind !=p){
      return(cbind(X[,1:(Cond.ind-1)],X0,X[,(Cond.ind):(p-1)]))
    }else if(Cond.ind==1){
      return(cbind(X0,X[,(Cond.ind):(p-1)]))
      
    }else if(Cond.ind==p){
      return(cbind(X[,1:(Cond.ind-1)],X0))
      
    }
  }  
  
}
