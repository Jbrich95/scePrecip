
## Note that this scripts follows the methodology described in Richards et al. (2022b)

#Also note, a lot of RAM is required if n_real is large! Approximately 100GB for application in paper!

## This requires full spatial and marginal fits for both convective and nonconvective rainfall!
print(file.exists("Data/conv.Rdata","MarginalAnalysis/convProbNoRain.Rdata","MarginalAnalysis/convGPDfits.Rdata","DependenceAnalysis/convfullspatfit.Rdata", "MarginalAnalysis/convLaplace_Data.Rdata",
                  "Data/nonconv.Rdata","MarginalAnalysis/nonconvProbNoRain.Rdata","MarginalAnalysis/nonconvGPDfits.Rdata","DependenceAnalysis/nonconvfullspatfit.Rdata","MarginalAnalysis/nonconvLaplace_Data.Rdata"))

# Import all required packages
source("RequiredPackages.R")

source("src/rOneCondsite.R")
source("src/spatial_fit_funcs.R")

# Set simulation hyperparameters
n_s=500; tau=27.5  #Defining S_s
n_c=1250 #Defining A_c for nonconvective rainfall
n_real<-5.5e5 #number of realisations
bprime_conv_scale=8; bprime_front_scale=12 #Defining no. of simulates required for importance sampling algorithm

#First simulate nonconvective events

load("Data/nonconv.Rdata")

size_N=nrow(Data)

# Define aggregate region


cbPalette <- c("#999999","#E69F00","#56B4E9","#009E73",
               "#F0E442","#0072B2","#D55E00","#CC79A7") # Color palette for plot

plot(coords,xlab="",ylab="")
#Aggregate regions are defined as circular and centred within the domain. 
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
p=length(S_s.inds)


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
  p=prob_zero[i]
  if(plaplace(v)>=1-lambda){
    v.star[i]=q+qgpd(p=(plaplace(v)-(1-lambda))/lambda,loc=0, scale=gpd_pred[S_s.inds[i],1],
                     shape=gpd_pred[S_s.inds[i],2])
  }else if(plaplace(v)<1-lambda ){
    prob_scaled=(plaplace(v)-p)/(1-p-lambda)*mean(sub.Data[,i][sub.Data[,i]>0]<q)
    
    
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

boo=rbinom(n=n_real,size=1,prob=1-n.comp/dim(sub.Data)[1]  ) # bernoulli replications to determine if the i-th replication of Y(s) comes from  Y(s)|max(X(s))<=v or Y(s)|max(X(s))>v

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

Sim_Mat_U=plaplace(sub.Sim)
rm(sub.Sim)

#Sim on Original margins with orginal coordinate system
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
rm(Sim_Mat_U)


# Aggregate

max.exceed.inds=which(comp.inds==1)

#Full Simulate
Y_N=matrix(nrow=n_real,ncol=length(A.inds))

sub.Orig.belowmax=sub.Data[-max.exceed.inds,A.inds] #We sample from this data if boo==1, i.e, max(X(s))<v

orig.inds=sample(size=sum(boo==1),x=1:dim(sub.Orig.belowmax)[1],replace=T)

Y_N[1:sum(boo==1),]=sub.Orig.belowmax[orig.inds,]
Y_N[-c(1:sum(boo==1)),]=Sim_Orig

rm(Sim_Orig) #Free-up space

R_NA=rowMeans(Y_N) #Derive sample of R_{\mathcal{N},\mathcal{A}}
hist(R_NA)

save(R_NA,A.inds,file="AggregateAnalysis/R_NA.Rdata")



### Now simulate convective events


q=0.96

load("Data/conv.Rdata")
load("MarginalAnalysis/convProbNoRain.Rdata")
load("MarginalAnalysis/convGPDfits.Rdata")
load("DependenceAnalysis/convfullspatfit.Rdata")


Data=Data.Mat_all
rm(Data.Mat_all)
coords=coords
Data[Data<1e-5]=0
#Mod2 - scaleDL

rOneCondSite=function(par,v,c.vec,coord,Cond.ind,n.sims){
  
  p=dim(coord)[1]
  
  KA1=par[1]
  KA2=par[2]
  KB1=par[3]
  KB2=par[4]
  KB3=par[5]
  
  KM1=par[6]
  KM2=par[7]
  KM3=par[8]
  KS1=par[9]
  KS2=par[10]
  lambda=par[11]
  kappa=par[12]
  KD1=par[13]
  KD2=par[14]
  KD3=par[15]
  KD4=par[16]
  
  re.coord=coord[c(Cond.ind,(1:p)[-Cond.ind]),]
  re.c=c.vec[c(Cond.ind,(1:p)[-Cond.ind])]
  
  h.pairs=rdist.earth(re.coord,miles=F)
  diag(h.pairs)=rep(0,p)
  Cor=Matern(h.pairs,range=lambda,smoothness = kappa)
  h.Cond=h.pairs[1,]
  rm(h.pairs)
  
  alpha=exp(-(abs(h.Cond)/KA1)^KA2)
  alpha=alpha[-1]
  
  beta=KB3*exp(-(abs(h.Cond)/KB1)^KB2)
  beta=beta[-1]
  
  sig=sqrt(2)*(1-exp(-(h.Cond/KS1)^KS2))
  sig=sig[-1]
  mu=KM1*h.Cond^{KM2-1}*exp(-h.Cond/KM3)
  mu=mu[-1]
  delta=1+(KD1*h.Cond^{KD2-1}-KD4)*exp(-h.Cond/KD3)
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

plot(coords)

centre.ind=which(coords[,1]==centre.coords[1])
S_s.inds=which(rdist.earth(coords,miles=F)[centre.ind,]<max_agg_range+tau)

points(coords[S_s.inds,],col="red")

S_s=coords[S_s.inds,]

set.seed(2)
S_s.inds=c(S_s.inds,sample((1:dim(coords)[1])[-S_s.inds],n_s))

S_s=coords[S_s.inds,]

A.inds=which(rdist.earth(S_s,miles=F)[which(rdist(rbind(centre.coords,S_s))[1,-1]==min(rdist(rbind(centre.coords,S_s))[1,-1])),]<max_agg_range)

A=S_s[A.inds,]
points(coords[S_s.inds,],col="blue")
points(A,col="green")

sub.Data=Data[,S_s.inds]
rm(Data)
p=length(S_s.inds)
#Sim level on Laplace - Check on original scale
v=qlaplace(q)

GPDscale=gpd_pred[,1]
GPDshape=gpd_pred[,2]
GPDthresh=gpd_pred[,3]


load("Conv_Laplace_Data_Large4.Rdata")
unif=function(x) rank(x)/(length(x)+1)
rm(Dat_Lap)

#Sim level on Laplace - Check on original scale

#Subset by max > v

v.star=rep(0,p)

for(i in 1:p){
  
  threshq=1-sum(sub.Data[,i]>GPDthresh[S_s.inds[i]])/length(sub.Data[,i])
  thresh=GPDthresh[S_s.inds[i]]
  if(plaplace(v)>=threshq){
    
    v.star[i]=thresh+qgpd(p=(plaplace(v)-threshq)/(1-threshq),loc=0, scale=GPDscale[S_s.inds[i]],
                          shape=GPDshape[S_s.inds[i]])
  }else if(plaplace(v)<threshq ){
    scale.p = (plaplace(v)-plaplace(c.vec[S_s.inds[i]]))/(1-plaplace(c.vec[S_s.inds[i]]))
    v.star[i]=quantile(sub.Data[,i][sub.Data[, i]>0],prob=scale.p)
    
  }  
  
}
print(v.star[1])

comp.inds=apply(sub.Data,1,function(x){
  if(sum(x>v.star) >= 1){
    return(1)}else{return(0)}
})

n.comp=sum(comp.inds)
print(n.comp)

plot(coords)

points(S_s,col="red")

boo=rbinom(n=n_real,size=1,prob=1-n.comp/dim(sub.Data)[1]  )


b=bprime_conv_scale*sum(boo==0)

Sim_Mat=matrix(NA,nrow=b,ncol=p)
import.prob=rep(0,b)
Cond.sites=rep(0,b)
n.cond.ind=1:p
for(i in 1:b){
  Cond.sites[i]=sample(1:p,1)
}
for(i in 1:p){
  n.cond.ind[i]=sum(Cond.sites==i)
}
#One Cond site


if(jobid!=0){
  load(paste0("Boot_opts/Conv_bootopt1_",jobid,".Rdata"))
  
  
  pars = bootopt$par
  
  pars=c(pars[1:4],1,pars[-c(1:4,15:16)],1,pars[15:16])
  simpar=pars[-c(17,18)]
}else{
  load(paste("Conv_FullOpt_newd_1q960_4.Rdata",sep=""))
  
  
  pars = optconv_full_aniso$par
  simpar=pars[-c(17,18)]
}
anisotransform=function(s,theta,L){
  return(matrix(c(1,0,0,1/L),nrow=2,ncol=2,byrow=T)%*%matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow=2,ncol=2,byrow=T)%*%as.matrix(s))
  
}

##Sim over sub domain
cl=makeCluster(20)
setDefaultCluster(cl=cl)
invisible(clusterEvalQ (cl , library("evd")))
invisible(clusterEvalQ (cl , library("fields")))
invisible(clusterEvalQ (cl , library("mvnfast")))

clusterExport(cl , "c.vec")
clusterExport(cl , "simpar")
clusterExport(cl , "rOneCondSite")
clusterExport(cl , "rmvdlaplace")
clusterExport(cl , "qdlaplace")
clusterExport(cl , "v")
clusterExport(cl , "S_s.inds")
clusterExport(cl , "deltafunc")

tcoord=apply(S_s,1,function(x){
  anisotransform(x,theta=pars[17],L=pars[18])
})
tcoord=t(tcoord)
clusterExport(cl , "tcoord")

temp=cbind(1:p,n.cond.ind)
Sim_Mat_list=parApply(cl=cl,temp,1,function(x){
  if(as.numeric(x[2]) > 0){
    Sim=rOneCondSite(par=simpar,v=v,c.vec=c.vec[S_s.inds],coord=tcoord,n.sims=as.numeric(x[2]),Cond.ind=as.numeric(x[1]))
  }else{
    Sim=NULL
  }
  return(Sim)
  
})
stopCluster(cl)

Sim_Mat=matrix(NA,nrow=b,ncol=p)
int=1
for(i in 1:p){
  if(n.cond.ind[i]!=0){
    Sim_Mat[int:(int+n.cond.ind[i]-1),]=Sim_Mat_list[[i]]
  }
  int=int+n.cond.ind[i]
}
rm(Sim_Mat_list)


import.prob=apply(Sim_Mat,1,function(x){
  1/sum(x>v)
})
import.prob[is.infinite(import.prob)]=0
import.prob=import.prob/sum(import.prob)


#subset Sims importance sampling
samp.inds = sample(1:length(import.prob),size=sum(boo==0),prob=import.prob)
sub.Sim=Sim_Mat[samp.inds,]
rm(Sim_Mat)

sub.Sim=sub.Sim[,A.inds]
Sim_Mat_U=plaplace(sub.Sim)
rm(sub.Sim)

#Sim on Original margins with orginal coordiante system
Sim_Orig=apply(rbind(A.inds,Sim_Mat_U),2,function(x){
  
  sub.coord.ind=x[1]
  
  
  thresh=GPDthresh[S_s.inds[sub.coord.ind]]
  threshq=1-sum(sub.Data[,sub.coord.ind]>thresh)/length(sub.Data[,sub.coord.ind])
  
  probs=x[-1]
  
  ind.above=which(probs > threshq)
  
  if(length(ind.above)!=0){
    prob.above=probs[ind.above]
    orig.above=as.numeric(thresh)+qgpd(p=(prob.above-threshq)/(1-threshq),loc=0, scale=GPDscale[S_s.inds[sub.coord.ind]],
                                       shape=GPDshape[S_s.inds[sub.coord.ind]])
  }
  ps=plaplace(c.vec[S_s.inds[sub.coord.ind]])
  
  non.zero.inds=which(probs > ps)
  probs.below=(probs[non.zero.inds]-ps)/(1-ps)
  # 
  
  orig.below=quantile(sub.Data[,sub.coord.ind][sub.Data[,sub.coord.ind]>0],prob=probs.below)
  
  all.orig=rep(0,length(probs))
  all.orig[non.zero.inds]=orig.below
  # if(length(ind.above)!=0){
  #  all.orig[ind.above]=orig.above
  #}
  return(all.orig)
})
rm(Sim_Mat_U)


################################################3
##



# plot(qrdata,qrsim,asp=1,main="Above max")
# abline(a=0,b=1,col="red")

comp.inds=which(comp.inds==1)

#Full Simulate
Full_Sim2=matrix(nrow=n_real,ncol=length(A.inds))

sub.Orig.belowmax=sub.Data[-comp.inds,1:length(A.inds)]

orig.inds=sample(size=sum(boo==1),x=1:dim(sub.Orig.belowmax)[1],replace=T)

Full_Sim2[1:sum(boo==1),]=sub.Orig.belowmax[orig.inds,]

Full_Sim2[-c(1:sum(boo==1)),]=Sim_Orig

rm(Sim_Orig)


sub.Data=sub.Data[,A.inds]

centre.ind=which(A[,1]==centre.coords[1])

p.max=1-1/dim(sub.Data)[1]
probs=seq(0.8,p.max,length=5000)
agg_boot=matrix(nrow=length(probs),ncol  =4)


#pdf(file=paste("CombinedAgg_Conv_4_v",version,".#pdf",sep=""),width=10,height=10)
par(mfrow=c(2,2))
for(it in 1:4){
  
  sub.agg.inds=which(rdist.earth(A,miles=F)[centre.ind,]<radii[it])
  
  # plot(coords)
  # points(S_s,col="red")
  # points(S_s[agg.inds,],col="green")
  R_Orig=rowSums(sub.Data[,sub.agg.inds])/length(sub.agg.inds)
  R_Sim=rowSums(Full_Sim2[,sub.agg.inds])/length(sub.agg.inds)
  
  
  qrdata=quantile(R_Orig,prob=probs)
  qrsim=quantile(R_Sim,prob=probs)
  title=paste("Region",it)
  agg_boot[,it]=qrsim
  
  plot(qrdata,qrsim,ylim=range(qrdata,qrsim),xlim=range(qrdata,qrsim),main="",xlab="",ylab="")
  mtext(side=3, title,cex = 1.5,line=0.5)
  mtext(side=2, "Model",cex = 1.5, line=2.6)
  mtext(side=1, "Empirical",cex = 1.5, line=2.6)
  abline(a=0,b=1,col="red")
}
#dev.off()
save(agg_boot, file = paste("Agg_boot/agg_size_rl",jobid,"_conv.Rdata",sep=""))


par(mfrow=c(1,1))
lat_range=range(coords[,2])+c(-0.6,0.6)
lon_range=range(coords[,1])+c(-0.6,0.6)
plotmap     <- map('worldHires', xlim=lon_range,ylim=lat_range,interior=F,plot=FALSE)


plot(coords,ylab="",xlab="",ylim=lat_range,xlim=lon_range)
points(plotmap,type="l")
points(S_s,col="orange")
for(it in 4:1){
  
  sub.agg.inds=which(rdist.earth(A,miles=F)[centre.ind,]<radii[it])
  
  points(A[sub.agg.inds,],col=it+1,pch=20)
}

rm(sub.Data)
rm(sub.Orig.belowmax)
load("All_Original_Data_Large_Masked4.Rdata")


prob_is_frontal= size_N/43200

n.all=min(n_real,n_real)

num_front=rbinom(n=n.all,size=1,prob=prob_is_frontal )

front_inds=sample(1:n_real,size=sum(num_front),replace=F)
conv_inds=sample(1:n_real,size=sum(num_front==0),replace=F)

Combined_Sim=rbind(Y_N[front_inds,],Full_Sim2[conv_inds,])

centre.ind=which(coords[,1]==centre.coords[1])
agg.inds=which(rdist.earth(coords,miles=F)[centre.ind,]<max_agg_range)
sub.Data=Data.Mat_all[,agg.inds]
agg.coords=coords[agg.inds,]
centre.ind=which(agg.coords[,1]==centre.coords[1])

#save(Combined_Sim,file="Com_Sim.Rdata")

#pdf(file=paste("CombinedAgg_Both_4_v",version,".#pdf",sep=""),width=10,height=10)
p.max=1-1/dim(sub.Data)[1]
probs=seq(0.8,p.max,length=5000)



agg_boot=matrix(nrow=length(probs),ncol  =4)
front_exceed_prob=matrix(nrow=2,ncol=4)
agg_Lambda0.95=matrix(nrow=2,ncol=4)
agg_Lambda0.99=matrix(nrow=2,ncol=4)

par(mfrow=c(2,2))
for(it in 1:4){
  
  sub.agg.inds=which(rdist.earth(agg.coords,miles=F)[centre.ind,]<radii[it])
  
  # plot(coords)
  # points(S_s,col="red")
  # points(S_s[agg.inds,],col="green")
  R_Orig=rowSums(sub.Data[,sub.agg.inds])/length(sub.agg.inds)
  R_Sim=rowSums(Combined_Sim[,sub.agg.inds])/length(sub.agg.inds)
  
  
  
  R0.99=quantile(R_Sim,0.99)
  event_inds=which(R_Sim > R0.99)
  print(paste("Region",it,"Prob_Frontal for R > R_0.99"))
  
  print(mean(event_inds <= length(front_inds)))
  front_exceed_prob[1,it]=mean(event_inds <= length(front_inds))
  
  R0.995=quantile(R_Sim,0.995)
  event_inds=which(R_Sim > R0.995)
  print(paste("Region",it,"Prob_Frontal for R > R_0.995"))
  front_exceed_prob[2,it]=mean(event_inds <= length(front_inds))
  
  print(mean(event_inds <= length(front_inds)))
  
  qrdata=quantile(R_Orig,prob=probs)
  qrsim=quantile(R_Sim,prob=probs)
  
  agg_boot[,it]=qrsim
  
  
  print("Agg")
  
  pmin=0.95
  n_p=43200*(1-pmin)
  ps=pmin+(1:n_p)/(n_p+1)*(1-pmin)
  qrdata=quantile(R_Orig,prob=ps)
  qrsim=quantile(R_Sim,prob=ps)
  print(paste0("pmin = ",pmin))
  print(mean(abs(qrdata-qrsim)))
  print(mean((qrdata-qrsim)^2))
  
  agg_Lambda0.95[1,it]=mean(abs(qrdata-qrsim))
  agg_Lambda0.95[2,it]=mean((qrdata-qrsim)^2)
  
  pmin=0.99
  n_p=43200*(1-pmin)
  ps=pmin+(1:n_p)/(n_p+1)*(1-pmin)
  qrdata=quantile(R_Orig,prob=ps)
  qrsim=quantile(R_Sim,prob=ps)
  print(paste0("pmin = ",pmin))
  print(mean(abs(qrdata-qrsim)))
  print(mean((qrdata-qrsim)^2))
  
  agg_Lambda0.99[1,it]=mean(abs(qrdata-qrsim))
  agg_Lambda0.99[2,it]=mean((qrdata-qrsim)^2)
}
#dev.off()
save(agg_boot, file = paste("Agg_boot/agg_size_rl",jobid,".Rdata",sep=""))
save(front_exceed_prob, file = paste("Agg_boot/agg_front_exceed_prob",jobid,".Rdata",sep=""))
save(agg_Lambda0.99, file = paste("Agg_boot/agg_Lambdaq99_",jobid,".Rdata",sep=""))
save(agg_Lambda0.95, file = paste("Agg_boot/agg_Lambdaq95_",jobid,".Rdata",sep=""))
