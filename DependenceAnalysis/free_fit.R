#Import all required packages
source("RequiredPackages.R")


source("src/DL_funcs.r")
source("src/free_nll.R")

##Load required Rdata
load("Data/Data.Rdata")
rm(Data)
load("MarginalAnalysis/Laplace_Data.Rdata")

##Loaded Objects:
#  
# For n observed fields with d sampling locations
#
#  Dat_Lap: n x d matrix of data on Laplace scale See marginal_transform.R
#  c.vec: d-vector of censoring thresholds on Laplace scale
#  coords: d x 2 matrix of lon/lat coordinates

#Get indices of conditioning sites. We pick five randomly.
cond_inds=sample(1:nrow(coords),5)

#We estimate the ``free" parameter estimates of alpha, beta, sigma, mu and delta. See Figure 2 of the paper.

#Set exceedance threshold u. We take u as the 98% Laplace quantile
u=qlaplace(0.98)
#Include Delta-Laplace
alpha.free<-beta.free<-sig.free<-mu.free<-delta.free<-matrix(nrow=length(cond_inds),ncol=nrow(coords))

plot(1,type="n",ylab=expression(beta(h)),xlab=expression(h),main="",xlim=c(0,max(rdist.earth(coords,miles=F))),ylim=c(0,1))
for(j in 1:length(cond_inds)){ #for all conditioning sites
  
  exceeds=Dat_Lap[Dat_Lap[,cond_inds[j]]>u,]
  
for(i in 1:length(alpha.free)){ # for all other sites
    if(i == cond_inds[j]){
      alpha.free[j,i]=1
      beta.free[j,i]=1
      sig.free[j,i]=0
      mu.free[j,i]=0
      delta.free[j,i]=1
    }else{
    #We perform an iterative optimisation procedure by first fixing (beta = 0.5, delta = 2) and then allowing beta to vary with delta = 2 fixed, and finally allowing all parameters to vary.
    #This tends to improve optimisation speed and efficacy
      opt=optim(free_nll_Gauss_fixedbeta,X=exceeds[,cond_inds[j]],Y=exceeds[,i],par=c(0.1,0.5,0.1),c=c.vec[i],beta=0.5)
      opt=optim(free_nll_Gauss_fixedbeta,X=exceeds[,cond_inds[j]],Y=exceeds[,i],par=opt$par,method="BFGS",c=c.vec[i],beta=0.5)
      
    opt=optim(free_nll_Gauss,X=exceeds[,cond_inds[j]],Y=exceeds[,i],par=c(opt$par[1],0.5,opt$par[2:3]),c=c.vec[i])
    opt=optim(free_nll_Gauss,X=exceeds[,cond_inds[j]],Y=exceeds[,i],par=opt$par,method="BFGS",c=c.vec[i])
      
    opt=optim(free_nll_DL,X=exceeds[,cond_inds[j]],Y=exceeds[,i],par=c(opt$par,2),c=c.vec[i])
    opt=optim(free_nll_DL,X=exceeds[,cond_inds[j]],Y=exceeds[,i],par=opt$par,method="BFGS",c=c.vec[i])
    
    
    
    #print(opt)
    alpha.free[j,i]=opt$par[1]
    beta.free[j,i]=opt$par[2]
    sig.free[j,i]=opt$par[3]
    mu.free[j,i]=opt$par[4]
    delta.free[j,i]=opt$par[5]
    }
  points(rdist.earth(coords,miles=F)[cond_inds[j],i],beta.free[j,i],col="green")
  
  print(i)
}
}

#Overlay all estimates
plot(1,type="n",ylab=expression(alpha(h)),xlab=expression(h),main="",xlim=c(0,max(rdist.earth(coords,miles=F))),ylim=c(0,1))
for(j in 1:length(cond_inds)) points(rdist.earth(coords,miles=F)[cond_inds[j],],alpha.free[j,])

plot(1,type="n",ylab=expression(beta(h)),xlab=expression(h),main="",xlim=c(0,max(rdist.earth(coords,miles=F))),ylim=c(0,1))
for(j in 1:length(cond_inds)) points(rdist.earth(coords,miles=F)[cond_inds[j],],beta.free[j,])

plot(1,type="n",ylab=expression(sigma(h)),xlab=expression(h),main="",xlim=c(0,max(rdist.earth(coords,miles=F))),ylim=range(sig.free,na.rm=T))
for(j in 1:length(cond_inds)) points(rdist.earth(coords,miles=F)[cond_inds[j],],sig.free[j,])

plot(1,type="n",ylab=expression(mu(h)),xlab=expression(h),main="",xlim=c(0,max(rdist.earth(coords,miles=F))),ylim=range(mu.free,na.rm=T))
for(j in 1:length(cond_inds)) points(rdist.earth(coords,miles=F)[cond_inds[j],],mu.free[j,])

plot(1,type="n",ylab=expression(delta(h)),xlab=expression(h),main="",xlim=c(0,max(rdist.earth(coords,miles=F))),ylim=c(0,max(delta.free,na.rm=T)))
for(j in 1:length(cond_inds)) points(rdist.earth(coords,miles=F)[cond_inds[j],],delta.free[j,])

