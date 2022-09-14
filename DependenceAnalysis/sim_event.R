# Import all required packages
source("RequiredPackages.R")

source("src/rOneCondsite.R")
source("src/spatial_fit_funcs.R")

# Are you using the mixture model as defined in Richards et al. (2022b)? If so, set mix.boo==T. 
# Do you want the convective or non-convective fit? For the former, set conv.boo==T and for the latter, set conv.boo==F.
mix.boo <- T; conv.boo <-T #If using the model of Richards et al. (2022a), set mix.boo==F.

## Load required Rdata
if(mix.boo==F){load("Data/Data.Rdata");load("MarginalAnalysis/ProbNoRain.Rdata");load("MarginalAnalysis/GPDfits.Rdata");load("DependenceAnalysis/fullspatfit.Rdata")
}else if(mix.boo==T & conv.boo == T){ load("Data/conv.Rdata");load("MarginalAnalysis/convProbNoRain.Rdata");load("MarginalAnalysis/convGPDfits.Rdata");load("DependenceAnalysis/convfullspatfit.Rdata")
}else if(mix.boo==T & conv.boo == F){load("Data/nonconv.Rdata");load("MarginalAnalysis/nonconvProbNoRain.Rdata");load("MarginalAnalysis/nonconvGPDfits.Rdata");load("DependenceAnalysis/nonconvfullspatfit.Rdata")}


##Loaded Objects:
#  
# For n observed fields with d sampling locations
#
#  Data: n x d matrix of hourly precipitation rate (mm/hour). We set all values <= 1e-5 to 0.
#  coords: d x 2 matrix of lon/lat coordinates
#  gpd_pred: a d x 3 vector of marginal parameter estimates. Each row corresponds to estimates at one site.
#  gpd_pred[,1]: d-vector of estimated GPD scale parameters
#  gpd_pred[,2]: d-vector of estimated GPD shape parameters
#  gpd_pred[,3]: d-vector of estimated GPD exceedance threshold
#  prob_zero: d-vector of estimated probability of no rain at each site
#  spat_pars: Vector of parameter estimates for the full spatial model, see spatial_fit.R

## Required for mapping
ncols=120
col1<-(sequential_hcl(ncols,"oslo"))
col1=col1[-c(2:6)]

lat_range=range(coords[,2])
lon_range=range(coords[,1])
plotmap     <- map('worldHires', xlim=lon_range,ylim=lat_range,interior=F,plot=FALSE)


mappanel <- function(x,y,...) { # I think this is the fn that plots the map lines within the call to contour
  panel.contourplot(x,y,...,xaxt="n",axes=F)
  llines(plotmap$x,plotmap$y,col="black")
  #map.grid(plotmap)
}


r=0.073 # r changes the resolution of map. Each coordinate is rounded to the nearest r for plotting


c.vec = qlaplace(prob_zero) #Censoring thresholds

#Choose conditioning site
Cond.ind=sample(1:nrow(coords),1)

#Set exceedance threshold u. 
if(mix.boo==F) u=qlaplace(0.98) #We take u as the 98% Laplace quantile in Richards et al. (2022a)
if(mix.boo==T & conv.boo == T) u=qlaplace(0.96) #We take u as the 96% Laplace quantile for convective rainfall in Richards et al. (2022b)
if(mix.boo==T & conv.boo == F) u=qlaplace(0.99) #We take u as the 96% Laplace quantile for nonconvective rainfall in Richards et al. (2022b)

#For simulations, require v >= u. We set v = u
v=u


simpar=spat_pars[-c(length(spat_pars)-1,length(spat_pars))] #Remove anisotropy parameters and perform transformation

tcoord=t(apply(coords,1,function(x){
  anisotransform(x,theta=spat_pars[length(spat_pars)-1],L=spat_pars[length(spat_pars)])
}))

#Get n.sims realisations
n.sims=50
if(mix.boo==F) Sim=rOneCondsite(par=simpar,v=v,c.vec=c.vec,coords=tcoord,n.sims=n.sims,Cond.ind=Cond.ind)
if(mix.boo==T & conv.boo == T) Sim=rOneCondsite_conv(par=simpar,v=v,c.vec=c.vec,coords=tcoord,n.sims=n.sims,Cond.ind=Cond.ind)
if(mix.boo==T & conv.boo == F) Sim=rOneCondsite_nonconv(par=simpar,v=v,c.vec=c.vec,coords=tcoord,n.sims=n.sims,Cond.ind=Cond.ind)


#Plot the t.ind realisation on Laplace scale

t.ind <- which(Sim==max(Sim),arr.ind=T)[1] # We take the most extreme event
c    <- 0
c    <- contourplot(Sim[t.ind,]~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(tck=c(0,0),draw=F), 
                    panel=mappanel, aspect=1, region=T, main= "",
                    contour=F, pretty=F, cuts=length(col1)-1,   col.regions=col1 )

print(c)

Sim_U=plaplace(Sim) #Transform to Uniform margins

#Back-transform to original margins

lambda=0.005 #Exceedance (above q) probability

unif=function(x) rank(x)/(length(x)+1)
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


#Plot the t.ind realisation on original scale

c    <- 0

plt=Sim_Orig[t.ind,]
plt[plt==0]=NA
c    <- contourplot(plt~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(tck=c(0,0),draw=F), 
                    panel=mappanel, aspect=1, region=T, main= "",
                    contour=F, pretty=F, cuts=length(col1)-1,   col.regions=col1 )

print(c)

