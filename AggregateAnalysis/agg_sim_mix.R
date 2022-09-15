
## Note that this scripts follows the methodology described in Richards et al. (2022b)

#Also note, a lot of RAM is required if n.reals is large! Approximately 100GB for application in paper!

## This requires full spatial and marginal fits for both convective and nonconvective rainfall!
print(file.exists("Data/conv.Rdata","MarginalAnalysis/convProbNoRain.Rdata","MarginalAnalysis/convGPDfits.Rdata","DependenceAnalysis/convfullspatfit.Rdata",
                  "Data/nonconv.Rdata","MarginalAnalysis/nonconvProbNoRain.Rdata","MarginalAnalysis/nonconvGPDfits.Rdata","DependenceAnalysis/nonconvfullspatfit.Rdata"))

# Import all required packages
source("RequiredPackages.R")

source("src/rOneCondsite.R")
source("src/spatial_fit_funcs.R")

cbPalette <- c("#999999","#E69F00","#56B4E9","#009E73",
               "#F0E442","#0072B2","#D55E00","#CC79A7") # Color palette for plot

# Set simulation hyperparameters
n_s=500; tau=27.5  #Defining S_s
n_c=1250 #Defining A_c for nonconvective rainfall
n.reals<-5.5e5 #number of realisations
bprime_conv_scale=8; bprime_front_scale=12 #Defining no. of simulates required for importance sampling algorithm

##  First simulate nonconvective events

load("Data/nonconv.Rdata")

size_N=nrow(Data)

# Define aggregate region


plot(coords,xlab="",ylab="")
#Aggregate regions are defined as circular and centered within the domain. 
agg.radius = 42.5

centre.ind=c(mean(coords[,1]),mean(coords[,2]))
centre.ind=which(rdist(rbind(centre.ind,coords))[1,-1]==min(rdist(rbind(centre.ind,coords))[1,-1]))

centre.coords=coords[centre.ind,]
S_s.inds=which(rdist.earth(coords,miles=F)[centre.ind,]<agg.radius+tau)

set.seed(2)
S_s.inds=c(S_s.inds,sample((1:dim(coords)[1])[-S_s.inds],n_s))
S_s=coords[S_s.inds,]

points(S_s,col=cbPalette[1])
A.inds=which(rdist.earth(S_s,miles=F)[which(rdist(rbind(centre.coords,S_s))[1,-1]==min(rdist(rbind(centre.coords,S_s))[1,-1])),]<agg.radius)

A=S_s[A.inds,]
points(A,col=cbPalette[8])

sub.Data=Data[,S_s.inds]
d=length(S_s.inds)


#Make A_c coords
lon_range=range(coords[,1])
lat_range=range(coords[,2])
new.coords=matrix(NA,nrow=n_c,ncol=2)
for(i in 1:n_c){
  boo=0
  while(boo ==0){
    temp=c(rnorm(1,centre.coords[1],1.7),rnorm(1,centre.coords[2],1.7))
    if( (temp[1]> lon_range[1] & temp[1] < lon_range[2]) & (temp[2]> lat_range[1] & temp[2] < lat_range[2]) ){
      boo=0
    }else{boo = 1}
  }
  new.coords[i,]=temp
}
plot(coords,ylim=range(new.coords[,2]),xlim=range(new.coords[,1]),asp=1,xlab="",ylab="")
points(S_s,col=cbPalette[1])
points(new.coords,col=cbPalette[8])

A_c=rbind(S_s,new.coords)

load("MarginalAnalysis/nonconvProbNoRain.Rdata")
load("MarginalAnalysis/nonconvGPDfits.Rdata")
load("DependenceAnalysis/nonconvfullspatfit.Rdata")

c.vec = qlaplace(prob_zero) #Censoring thresholds

#Set exceedance threshold u. We take u as the 98% Laplace quantile
u=qlaplace(0.99)

#For simulations, require v >= u. We set v = u
v=u

lambda=0.005 #Exceedance (above q) probability

v.star=rep(0,length(S_s.inds)) #Transform v to the sitewise original margins
for(i in 1:length(v.star)){
  
  q=as.numeric(gpd_pred[S_s.inds[i],3])
  d=prob_zero[i]
  if(plaplace(v)>=1-lambda){
    v.star[i]=q+qgpd(d=(plaplace(v)-(1-lambda))/lambda,loc=0, scale=gpd_pred[S_s.inds[i],1],
                     shape=gpd_pred[S_s.inds[i],2])
  }else if(plaplace(v)<1-lambda ){
    prob_scaled=(plaplace(v)-d)/(1-d-lambda)*mean(sub.Data[,i][sub.Data[,i]>0]<q)
    
    
    v.star[i]=quantile(sub.Data[,i][sub.Data[,i]>0],prob=prob_scaled)
  }
}
#Calculate the proportion of observations  X(s)|max(X(s))>v or similarly Y(s)|max(Y(s) - v.star) >0
comp.inds=apply(sub.Data,1,function(x){
  if(sum(x>v.star) >= 1){
    return(1)}else{return(0)}
})
n.comp=sum(comp.inds)
print(n.comp)

boo=rbinom(n=n.reals,size=1,prob=1-n.comp/dim(sub.Data)[1]  ) # bernoulli replications to determine if the i-th replication of Y(s) comes from  Y(s)|max(X(s))<=v or Y(s)|max(X(s))>v

#We require b replications of Y(s)|max(X(s))>v. We simulate bprime=bprime_front_scale * b  and then thin the replications using importance sampling

b=bprime_front_scale*sum(boo==0)
size.A_C=dim(A_c)[1]
Sim_Mat=matrix(NA,nrow=b,ncol=size.A_C)
import.prob=rep(0,b)
Cond.sites=rep(0,b)
n.cond.ind=1:size.A_C
for(i in 1:b){
  Cond.sites[i]=sample(1:size.A_C,1) #randomly sample conditioning sites
}
for(i in 1:size.A_C){
  n.cond.ind[i]=sum(Cond.sites==i)
}


#Much faster if parallelised
ncores<-5
#ncores <- detectCores()
cl=makeCluster(ncores)
setDefaultCluster(cl=cl)
invisible(clusterEvalQ (cl , library("evd")))
invisible(clusterEvalQ (cl , library("fields")))
invisible(clusterEvalQ (cl , library("mvnfast")))

clusterExport(cl , "c.vec")
clusterExport(cl , "rOneCondsite_nonconv")
clusterExport(cl , "rmvdlaplace")
clusterExport(cl , "qdlaplace")
clusterExport(cl , "v")
clusterExport(cl , "S_s.inds")
clusterExport(cl , "deltafunc")
clusterExport(cl , "beta_nc_func")
clusterExport(cl , "beta_nc_funcderv")

simpar=spat_pars[-c(length(spat_pars)-1,length(spat_pars))] #Remove anisotropy parameters and perform transformation
clusterExport(cl , "simpar")

tcoord=apply(A_c,1,function(x){
  anisotransform(x,theta=spat_pars[length(spat_pars)],L=spat_pars[length(spat_pars)-1])
})
tcoord=t(tcoord)

clusterExport(cl , "tcoord")

A_C.c.vec=c(c.vec[S_s.inds],rep(-Inf,n_c))
clusterExport(cl , "A_C.c.vec")

tmp=cbind(1:size.A_C,n.cond.ind)

#Draw n.cond.ind replications of X(s)|X(s_O)>v for each conditioning site in A_c
Sim_Mat_list=parApply(cl=cl,tmp,1,function(x){
  if(as.numeric(x[2]) > 0){
    Sim=rOneCondsite_nonconv(par=simpar,v=v,c.vec=A_C.c.vec,coord=tcoord,n.sims=as.numeric(x[2]),Cond.ind=as.numeric(x[1]))
  }else{
    Sim=NULL
  }
  return(Sim)
  
})
stopCluster(cl)

Sim_Mat=matrix(NA,nrow=b,ncol=size.A_C)
int=1
for(i in 1:size.A_C){
  if(n.cond.ind[i]!=0){
    Sim_Mat[int:(int+n.cond.ind[i]-1),]=Sim_Mat_list[[i]]
  }
  int=int+n.cond.ind[i]
}
rm(Sim_Mat_list)


#Calculate importance weights
import.prob=apply(Sim_Mat,1,function(x){
  if(sum(x > v) > 0){
    return(1/sum(x>v))
  }else{
    return(0)
  }
})

import.prob[is.infinite(import.prob)]=0

import.prob=import.prob/sum(import.prob)


#subset Sims importance sampling
samp.inds = sample(1:length(import.prob),size=sum(boo==0),prob=import.prob)
sub.Sim=Sim_Mat[samp.inds,]
rm(Sim_Mat)

sub.Sim=sub.Sim[,A.inds]

Sim_Mat_U=plaplace(sub.Sim) # On uniform margins
rm(sub.Sim)

#Sim on original margins with orginal coordinate system
Sim_Orig=apply(rbind(A.inds,Sim_Mat_U),2,function(x){
  coord_ind=x[1]
  q=as.numeric(gpd_pred[coord_ind,3])
  
  probs=x[-1]
  
  ind.above=which(probs >= 1-lambda)
  
  if(length(ind.above)!=0){
    prob.above=probs[ind.above]
    orig.above=q+qgpd(p=(prob.above-(1-lambda))/lambda,loc=0, scale=gpd_pred[coord_ind,1],
                      shape=gpd_pred[coord_ind,2])
  }
  p0=prob_zero[coord_ind]
  
  bulk.inds=which(probs > p0 & probs < (1-lambda))
  prob_scaled=(probs[bulk.inds]-p0)/(1-p0-lambda)*mean(Data[,coord_ind][Data[,coord_ind]>0]<q)
  
  orig.bulk=quantile(Data[,coord_ind][Data[,coord_ind]>0],prob=prob_scaled)
  
  all.orig=rep(0,length(probs))
  all.orig[bulk.inds]=orig.bulk
  
  if(length(ind.above)!=0){
    all.orig[ind.above]=orig.above
  }
  
  return(all.orig)
  
})
rm(Sim_Mat_U)


# Aggregate

max.exceed.inds=which(comp.inds==1)

#Full Simulate
Y_N=matrix(nrow=n.reals,ncol=length(A.inds))

sub.Orig.belowmax=sub.Data[-max.exceed.inds,A.inds] #We sample from this data if boo==1, i.e, max(X(s))<v

orig.inds=sample(size=sum(boo==1),x=1:dim(sub.Orig.belowmax)[1],replace=T) # Original data sample indices

Y_N[1:sum(boo==1),]=sub.Orig.belowmax[orig.inds,]
Y_N[-c(1:sum(boo==1)),]=Sim_Orig #Use simulates from SCE simulation

rm(Sim_Orig) #Free-up space

R_NA=rowMeans(Y_N) #Derive sample of R_{\mathcal{N},\mathcal{A}}
hist(R_NA)

save(R_NA,A.inds,file="AggregateAnalysis/R_NA.Rdata")



### Now simulate convective events


load("Data/conv.Rdata")
load("MarginalAnalysis/convProbNoRain.Rdata")
load("MarginalAnalysis/convGPDfits.Rdata")
load("DependenceAnalysis/convfullspatfit.Rdata")

A_c= S_s # A_c and S_s are equivalent for convective rainfall

c.vec = qlaplace(prob_zero) #Censoring thresholds

#Set exceedance threshold u. We take u as the 98% Laplace quantile
u=qlaplace(0.96)

#For simulations, require v >= u. We set v = u
v=u

lambda=0.005 #Exceedance (above q) probability

v.star=rep(0,length(S_s.inds)) #Transform v to the sitewise original margins
for(i in 1:length(v.star)){
  
  q=as.numeric(gpd_pred[S_s.inds[i],3])
  d=prob_zero[i]
  if(plaplace(v)>=1-lambda){
    v.star[i]=q+qgpd(d=(plaplace(v)-(1-lambda))/lambda,loc=0, scale=gpd_pred[S_s.inds[i],1],
                     shape=gpd_pred[S_s.inds[i],2])
  }else if(plaplace(v)<1-lambda ){
    prob_scaled=(plaplace(v)-d)/(1-d-lambda)*mean(sub.Data[,i][sub.Data[,i]>0]<q)
    
    
    v.star[i]=quantile(sub.Data[,i][sub.Data[,i]>0],prob=prob_scaled)
  }
}
#Calculate the proportion of observations  X(s)|max(X(s))>v or similarly Y(s)|max(Y(s) - v.star) >0
comp.inds=apply(sub.Data,1,function(x){
  if(sum(x>v.star) >= 1){
    return(1)}else{return(0)}
})
n.comp=sum(comp.inds)
print(n.comp)

boo=rbinom(n=n.reals,size=1,prob=1-n.comp/dim(sub.Data)[1]  ) # bernoulli replications to determine if the i-th replication of Y(s) comes from  Y(s)|max(X(s))<=v or Y(s)|max(X(s))>v

#We require b replications of Y(s)|max(X(s))>v. We simulate bprime=bprime_front_scale * b  and then thin the replications using importance sampling

b=bprime_conv_scale*sum(boo==0)
size.A_C=dim(A_c)[1]
Sim_Mat=matrix(NA,nrow=b,ncol=size.A_C)
import.prob=rep(0,b)
Cond.sites=rep(0,b)
n.cond.ind=1:size.A_C
for(i in 1:b){
  Cond.sites[i]=sample(1:size.A_C,1) #randomly sample conditioning sites
}
for(i in 1:size.A_C){
  n.cond.ind[i]=sum(Cond.sites==i)
}


#Much faster if parallelised
ncores<-5
#ncores <- detectCores()
cl=makeCluster(ncores)
setDefaultCluster(cl=cl)
invisible(clusterEvalQ (cl , library("evd")))
invisible(clusterEvalQ (cl , library("fields")))
invisible(clusterEvalQ (cl , library("mvnfast")))

clusterExport(cl , "c.vec")
clusterExport(cl , "rOneCondsite_conv")
clusterExport(cl , "rmvdlaplace")
clusterExport(cl , "qdlaplace")
clusterExport(cl , "v")
clusterExport(cl , "S_s.inds")
clusterExport(cl , "deltafunc")

simpar=spat_pars[-c(length(spat_pars)-1,length(spat_pars))] #Remove anisotropy parameters and perform transformation
clusterExport(cl , "simpar")

tcoord=apply(A_c,1,function(x){
  anisotransform(x,theta=spat_pars[length(spat_pars)],L=spat_pars[length(spat_pars)-1])
})
tcoord=t(tcoord)

clusterExport(cl , "tcoord")


tmp=cbind(1:size.A_C,n.cond.ind)

#Draw n.cond.ind replications of X(s)|X(s_O)>v for each conditioning site in A_c
Sim_Mat_list=parApply(cl=cl,tmp,1,function(x){
  if(as.numeric(x[2]) > 0){
    Sim=rOneCondsite_conv(par=simpar,v=v,c.vec=c.vec,coord=tcoord,n.sims=as.numeric(x[2]),Cond.ind=as.numeric(x[1]))
  }else{
    Sim=NULL
  }
  return(Sim)
  
})
stopCluster(cl)

Sim_Mat=matrix(NA,nrow=b,ncol=size.A_C)
int=1
for(i in 1:size.A_C){
  if(n.cond.ind[i]!=0){
    Sim_Mat[int:(int+n.cond.ind[i]-1),]=Sim_Mat_list[[i]]
  }
  int=int+n.cond.ind[i]
}
rm(Sim_Mat_list)


#Calculate importance weights
import.prob=apply(Sim_Mat,1,function(x){
  if(sum(x > v) > 0){
    return(1/sum(x>v))
  }else{
    return(0)
  }
})

import.prob[is.infinite(import.prob)]=0

import.prob=import.prob/sum(import.prob)


#subset Sims importance sampling
samp.inds = sample(1:length(import.prob),size=sum(boo==0),prob=import.prob)
sub.Sim=Sim_Mat[samp.inds,]
rm(Sim_Mat)

sub.Sim=sub.Sim[,A.inds]

Sim_Mat_U=plaplace(sub.Sim) # On uniform margins
rm(sub.Sim)

#Transform to original margins with original coordinate system
Sim_Orig=apply(rbind(A.inds,Sim_Mat_U),2,function(x){
  coord_ind=x[1]
  q=as.numeric(gpd_pred[coord_ind,3])
  
  probs=x[-1]
  
  ind.above=which(probs >= 1-lambda)
  
  if(length(ind.above)!=0){
    prob.above=probs[ind.above]
    orig.above=q+qgpd(p=(prob.above-(1-lambda))/lambda,loc=0, scale=gpd_pred[coord_ind,1],
                      shape=gpd_pred[coord_ind,2])
  }
  p0=prob_zero[coord_ind]
  
  bulk.inds=which(probs > p0 & probs < (1-lambda))
  prob_scaled=(probs[bulk.inds]-p0)/(1-p0-lambda)*mean(Data[,coord_ind][Data[,coord_ind]>0]<q)
  
  orig.bulk=quantile(Data[,coord_ind][Data[,coord_ind]>0],prob=prob_scaled)
  
  all.orig=rep(0,length(probs))
  all.orig[bulk.inds]=orig.bulk
  
  if(length(ind.above)!=0){
    all.orig[ind.above]=orig.above
  }
  
  return(all.orig)
  
})
rm(Sim_Mat_U)


# Aggregate

max.exceed.inds=which(comp.inds==1)

#Full Simulate
Y_C=matrix(nrow=n.reals,ncol=length(A.inds))

sub.Orig.belowmax=sub.Data[-max.exceed.inds,A.inds] #We sample from this data if boo==1, i.e, max(X(s))<v

orig.inds=sample(size=sum(boo==1),x=1:dim(sub.Orig.belowmax)[1],replace=T) #Indices for empirical obs

Y_C[1:sum(boo==1),]=sub.Orig.belowmax[orig.inds,]
Y_C[-c(1:sum(boo==1)),]=Sim_Orig #Use extreme fields

R_CA=rowMeans(Y_C) #Derive sample of R_{\mathcal{C},\mathcal{A}}
hist(R_CA)

save(R_CA,A.inds,file="AggregateAnalysis/R_CA.Rdata")

# Combine to get a sample from R_A

p_C= dim(Data)[1]/(dim(Data)[1]+size_N) #Prob of convective rainfall occuring

boo_conv=rbinom(n=n.reals,size=1,prob=p_C ) #Should the i-th sample of R_A come from R_NA or R_CA

nconv_inds=sample(1:n.reals,size=sum(boo_conv==0),replace=F) #Which realistions of Y_N should be used?
conv_inds=sample(1:n.reals,size=sum(boo_conv),replace=F) #Which realistions of Y_C should be used?

Y=rbind(Y_N[nconv_inds,],Y_C[conv_inds,]) #Samples from Y_M

R_A=rowMeans(Y)
save(R_A,A.inds,file="AggregateAnalysis/R_A.Rdata")
