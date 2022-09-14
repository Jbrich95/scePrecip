
## Note that this scripts follows the methodology described in Richards et al. (2022a)

#Import all required packages
source("RequiredPackages.R")

source("src/rOneCondsite.R")
source("src/spatial_fit_funcs.R")


##Load required Rdata
load("Data/Data.Rdata")
load("MarginalAnalysis/ProbNoRain.Rdata")
load("MarginalAnalysis/GPDfits.Rdata")
load("DependenceAnalysis/fullspatfit.Rdata")

##Loaded Objects:
#  #  
# For n observed fields with d sampling locations
#
#  Data: n x d matrix of hourly precipitation rate (mm/hour). We set all values <= 1e-5 to 0.
#  coords: d x 2 matrix of lon/lat coordinates
#  gpd_pred: a d x 3 vector of marginal parameter estimates. Each row corresponds to estimates at one site.
#  gpd_pred[,1]: d-vector of estimated GPD scale parameters
#  gpd_pred[,2]: d-vector of estimated GPD shape parameters
#  gpd_pred[,3]: d-vector of estimated GPD exceedance threshold
#  prob_zero: d-vector of estimated probability of no rain at each site
#  spat_pars: Vector of 16 parameter estimates for the full spatial AI model, see spatial_fit.R

c.vec = qlaplace(prob_zero) #Censoring thresholds

#Set exceedance threshold u. We take u as the 98% Laplace quantile
u=qlaplace(0.98)

#For simulations, require v >= u. We set v = u
v=u

simpar=spat_pars[-c(15,16)] #Remove anisotropy parameters and perform transformation

tcoord=t(apply(coords,1,function(x){
  anisotransform(x,theta=spat_pars[15],L=spat_pars[16])
}))

lambda=0.005 #Exceedance (above q) probability

v.star=rep(0,nrow(coords)) #Transform v to the sitewise original margins
for(i in 1:length(v.star)){
  
  q=as.numeric(gpd_pred[i,3])
  p=prob_zero[i]
  if(plaplace(v)>=1-lambda){
    v.star[i]=q+qgpd(p=(plaplace(v)-(1-lambda))/lambda,loc=0, scale=gpd_pred[i,1],
                    shape=gpd_pred[i,2])
  }else if(plaplace(v)<1-lambda ){
    prob_scaled=(plaplace(v)-p)/(1-p-lambda)*mean(Data[,i][Data[,i]>0]<q)
    

    v.star[i]=quantile(Data[,i][Data[,i]>0],prob=prob_scaled)
  }
}
#Calculate the proportion of observations  X(s)|max(X(s))>v or similarly Y(s)|max(Y(s) - v.star) >0
comp.inds=apply(Data,1,function(x){
  if(sum(x>v.star) >= 1){
    return(1)}else{return(0)}
})

n.comp=sum(comp.inds)
print(n.comp)
n.full=5e5 #the desired number of replications of R_\mathcal{A}
boo=rbinom(n=n.full,size=1,prob=1-n.comp/dim(Data)[1]  ) # bernoulli replications to determine if the i-th replication of Y(s) comes from  Y(s)|max(X(s))<=v or Y(s)|max(X(s))>v

#We require N=sum(boo==0) replications of Y(s)|max(X(s))>v. We simulate Nprime=5N and then thin the replications using importance sampling

N=sum(boo==0)
Nprime=5*N
p=dim(coords)[1]

Sim_Mat=matrix(NA,nrow=Nprime,ncol=p)
import.prob<-Cond.sites<-rep(0,Nprime)
n.cond.ind=1:p
for(i in 1:Nprime){
  Cond.sites[i]=sample(1:p,1) #randomly sample conditioning sites
}
for(i in 1:p){
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

clusterExport(cl , "rOneCondsite")
clusterExport(cl , "rmvdlaplace")
clusterExport(cl , "qdlaplace")
clusterExport(cl , "tcoord")
clusterExport(cl , "simpar")
clusterExport(cl , "c.vec")
clusterExport(cl , "v")

#Draw n.cond.ind replications of X(s)|X(s_O)>v for each conditioning site
Sim_Mat_list=parApply(cl=cl,(cbind(1:p,n.cond.ind)),1,function(x){
  if(as.numeric(x[2]) > 0){
  Sim=rOneCondsite(par=simpar,v=v,c.vec=c.vec,coords=tcoord,n.sims=as.numeric(x[2]),Cond.ind=as.numeric(x[1]))
  }else{
    Sim=NULL
  }
  return(Sim)
})
stopCluster(cl)

Sim_Mat=matrix(NA,nrow=Nprime,ncol=p)
int=1
for(i in 1:p){
  if(n.cond.ind[i]!=0){
    Sim_Mat[int:(int+n.cond.ind[i]-1),]=Sim_Mat_list[[i]]
  }
  int=int+n.cond.ind[i]
}
rm(Sim_Mat_list)


#Calculate importance weights
import.prob=apply(Sim_Mat,1,function(x){
  1/sum(x>v)
})

import.prob=import.prob/sum(import.prob)

#subset Sims importance sampling
samp.inds = sample(1:length(import.prob),size=sum(boo==0),prob=import.prob)

sub.Sim=Sim_Mat[samp.inds,]

Sim_U=plaplace(sub.Sim) #Transform to Uniform margins

#Back-transform to original margins


Sim_Orig=apply(rbind(1:nrow(coords),Sim_U),2,function(x){
  
  coord_ind=x[1]
  q=as.numeric(gpd_pred[coord_ind,3])
  
  probs=x[-1]
  
  ind.above=which(probs >= 1-lambda)
  
  if(length(ind.above)!=0){
    prob.above=probs[ind.above]
    orig.above=q+qgpd(p=(prob.above-(1-lambda))/lambda,loc=0, scale=gpd_pred[coord_ind,1],
                      shape=gpd_pred[coord_ind,2])
  }
  p=prob_zero[coord_ind]
  
  bulk.inds=which(probs > p & probs < (1-lambda))
  prob_scaled=(probs[bulk.inds]-p)/(1-p-lambda)*mean(Data[,coord_ind][Data[,coord_ind]>0]<q)
  
  orig.bulk=quantile(Data[,coord_ind][Data[,coord_ind]>0],prob=prob_scaled)
  
  all.orig=rep(0,length(probs))
  all.orig[bulk.inds]=orig.bulk
  
  if(length(ind.above)!=0){
    all.orig[ind.above]=orig.above
  }
  
  return(all.orig)
})

# Define aggregate region \mathcal{A} - Here we use a circular region with radius 16km centered in the domain

centre.ind=c(mean(coords[,1]),mean(coords[,2]))
centre.ind=which(rdist(rbind(centre.ind,coords))[1,-1]==min(rdist(rbind(centre.ind,coords))[1,-1]))
agg.inds=which(rdist.earth(rbind(coords),miles=F)[centre.ind,]<16)
             
plot(coords,xlab="",ylab='')
points(coords[agg.inds,],col="red")             
                       

              
max.exceed.inds=which(comp.inds==1)
#Full Simulate - Y(s):s \in \mathcal{A}
Y=matrix(nrow=n.full,ncol=length(agg.inds))

sub.Orig.belowmax=Data[-max.exceed.inds,agg.inds] #We sample from this data if boo==1, i.e, max(X(s))<v

orig.inds=sample(size=sum(boo==1),x=1:dim(sub.Orig.belowmax)[1],replace=T)

Y[1:sum(boo==1),]=sub.Orig.belowmax[orig.inds,]
Y[-c(1:sum(boo==1)),]=Sim_Orig[,agg.inds]

R_A=rowMeans(Y) #Derive sample of R_\mathcal{A}
hist(R_A)

save(R_A,agg.inds,file="AggregateAnalysis/R_A.Rdata")
