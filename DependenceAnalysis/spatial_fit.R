
source("src/DL_funcs.R")

# Import all required packages
source("RequiredPackages.R")

# Are you using the mixture model as defined in Richards et al. (2022b)?If so, set mix.boo==T. 
# Do you want the convective or non-convective fit? For the former, set conv.boo==T and for the latter, set conv.boo==F.
mix.boo <- T; conv.boo <-T #If using the model of Richards et al. (2022a), set mix.boo==F.

##Load required Rdata
if(mix.boo==F){load("MarginalAnalysis/Laplace_Data.Rdata")
}else if(mix.boo==T & conv.boo == T){ load("MarginalAnalysis/convLaplace_Data.Rdata")
}else if(mix.boo==T & conv.boo == F){load("MarginalAnalysis/nonconvLaplace_Data.Rdata") }


##Loaded objects:
# For n observed fields with d sampling locations
#
#  Data: n x d matrix of hourly precipitation rate (mm/hour). We set all values <= 1e-5 to 0.
#  coords: d x 2 matrix of lon/lat coordinates
#  Dat_Lap: n x d matrix of data on Laplace scale See marg_transform.R
#  c.vec: d-vector of censoring thresholds on Laplace scale


d_s=5000 #number of triples uses in pseudo-likelihood

if(mix.boo==F) h_max=28 #Maximum distance between conditioning sites and other sites within triple
if(mix.boo==T & conv.boo == T) h_max= 35
if(mix.boo==T & conv.boo == F) h_max=250
#Set exceedance threshold u. 
if(mix.boo==F) u=qlaplace(0.98) #We take u as the 98% Laplace quantile in Richards et al. (2022a)
if(mix.boo==T & conv.boo == T) u=qlaplace(0.96) #We take u as the 96% Laplace quantile for convective rainfall in Richards et al. (2022b)
if(mix.boo==T & conv.boo == F) u=qlaplace(0.99) #We take u as the 96% Laplace quantile for nonconvective rainfall in Richards et al. (2022b)
#Censor data
for(i in 1:dim(Dat_Lap)[2]){
  Dat_Lap[Dat_Lap[,i]<= c.vec[i],i]=NA
}

h=rdist.earth(coords,miles=F)
diag(h)=0


Exceed.Inds=c()

cond.site.inds=sample(1:nrow(coords),(d_s+100),replace=T) #Sample d_s+100 conditioning sites. Extra are sampled to ensure there are no duplicates
    
 for(j in 1:length(cond.site.inds)){
    
    dis.inds=which(h[cond.site.inds[j],] > 0.1 & h[cond.site.inds[j],]  <= h_max) #Find pairs of sites within h_max of conditioning site
    
    Pair.Inds=sample(dis.inds,2,replace=F)
    
    Exceed.Inds=rbind(Exceed.Inds,c(cond.site.inds[j],Pair.Inds))
  }

#We now remove duplicates and take only the first d_s triples 
new.Exceed.Inds=c()
uni.cond.inds=unique(Exceed.Inds[,1])
for(i in 1:length(uni.cond.inds)){
  temp=matrix(Exceed.Inds[which(Exceed.Inds[,1]==uni.cond.inds[i]),2:3],ncol=2)
  
  temp=matrix(temp[order(temp[,1]),],ncol=2)
  
  dup.inds=which(duplicated(temp))
  if(length(dup.inds)>0){
    temp=temp[-dup.inds,]
  }
  new.Exceed.Inds=rbind(new.Exceed.Inds,cbind(rep(uni.cond.inds[i],length=dim(temp)[1]),temp))
}

Exceed.Inds=new.Exceed.Inds[1:d_s,]


#Subset Laplace data and find data given exceedance above u at conditioning site

n.pairs=dim(Exceed.Inds)[1]
temp<-Exceed.pairs<-vector("list",n.pairs)
for( i in 1:n.pairs){Exceed.pairs[[i]]=vector("list",3)}

for(i in 1:n.pairs){
  Exceed.Dat=Dat_Lap[,c(Exceed.Inds[i,])]
  
  Exceed.Dat=Exceed.Dat[!is.na(Exceed.Dat[,1]),]
  Exceed.Dat=Exceed.Dat[Exceed.Dat[,1]>u,]
  
  temp[[i]]=Exceed.Dat
}

# We split the data into the 3 cases; no censoring, full censoring, one location censored
for(i in 1:n.pairs){
  ind1=which(is.finite(temp[[i]][,2])& is.finite(temp[[i]][,3]))
  Exceed.pairs[[i]][[1]]=temp[[i]][ind1,]
  ind2=which(is.na(temp[[i]][,2])& is.na(temp[[i]][,3]))
  Exceed.pairs[[i]][[2]]=temp[[i]][ind2,]
  Exceed.pairs[[i]][[3]]=temp[[i]][-c(ind1,ind2),]
  
  
}

rm(Dat_Lap)

source("src/spatial_fit_funcs.R")

#We now fit the full spatial model described in the paper. This can be particularly computationally intensive and so will require
#parallelised optimisation. The final fit for the paper took roughly 30 hours with 40 cpus.
ncores <- detectCores() # ncores <- 40
cl=makeCluster(ncores)
setDefaultCluster(cl=cl)
invisible(clusterEvalQ (cl , library("mvnfast")))
invisible(clusterEvalQ (cl , library("mvtnorm")))
invisible(clusterEvalQ (cl , library("mnormt")))
invisible(clusterEvalQ (cl , library("fields")))
clusterExport(cl , "dmvdlaplace")
clusterExport(cl , "ddlaplace")
clusterExport(cl , "pdlaplace")
clusterExport(cl , "anisotransform")


#Optimisation is conducted using box constraints. Whilst the lower parameter bounds are defined by the model and should not be changed,
#most of the upper bounds have been set to finite values rather than infinite ones. Our upper bounds were determined by fitting the model with lower d_s.

# Initial parameters, constraints and the nll function are specific to the paper and process


if(mix.boo==F){
#Please refer to spatmodnll_AI in spatial_fit_funcs.R for details of the full model and its parameters. We have set
#kapppa_beta_3=kappa_delta_4=1 (KB3=KD4=1) as in the paper, but this can be changed within the spatmodnll function.

lowpar=c(1e-4,  0, 
         1e-4,0, 
         -1,1e-5, 1e-4,
         1e-4, 0, 
         1e-4,1e-4,
         1e-4,1e-5,1e-4, 
         -pi/2,1e-4)

uppar=c(30,  3, 
        50, 3 ,
        1, 3,  150  ,
        50,  2,
        200,2,
        2,5,250,
        0,10)


#init.par: initial parameters taken from final estimates given in paper
init.par=c(1.95,  0.73, 
           38.6, 1.02 ,
           0.65, 0.28,  140  ,
           34.2,  0.89,
           58.7,0.53,
           0.43,0.46,142,
           -0.18,0.93)

opt=optimParallel(par=init.par,
                  lower=lowpar,
                  upper=uppar,
                  fn=spatmodnll_AI,Exceed.all=Exceed.pairs,Exceed.Inds=Exceed.Inds,coords=coords,
                  c.vec=c.vec,parallel=list(loginfo=TRUE))
print(opt)

stopCluster(cl)

spat_pars=opt$par
save(spat_pars,file="DependenceAnalysis/fullspatfit.Rdata")
}

## The intrinsically asymptotically-dependent spatial model from Richards et al. (2022a) can be similarly fitted by minimising spatmodnll_AD instead, see spatial_fit_funcs.R.


# lowpar=c( -1,1e-5, 1e-4, 
#          1e-4, 0,0, 
#          1e-4,1e-4,
#          1e-4,1e-5,1e-4, -1,
#          -pi/2,1e-4)
# uppar=c(1, 4,  400  
#         ,50,  2,5,
#         1000,2
#         ,2,5,1200,1,
#         0,10)
# 
# #init.par: initial parameters taken from final estimates given in paper
# init.par=c(-0.08, 0.86, 249.8, 
#           12.68, 0.63, 2.30,
#           655.03, 0.36, 
#           0.008, 1.17, 753.21, -0.03, 
#           -0.19, 0.95)
# opt=optimParallel(par=init.par,
#                       lower=lowpar,
#                       upper=uppar,
#                       fn=spatmodnll_AD,Exceed.all=Exceed.pairs,Exceed.Inds=Exceed.Inds,coords=coords,
#                       c.vec=c.vec,parallel=list(loginfo=TRUE))
# print(opt)
# 
# stopCluster(cl)
# 
# spat_pars=opt$par
# save(spat_pars,file="DependenceAnalysis/ADfit.Rdata")


## For mixture component models

if(mix.boo==T & conv.boo==T){
  #Please refer to spatmodnll_conv in spatial_fit_funcs.R for details of the full model and its parameters.
  
  lowpar=c(1e-4,  0, 0,
           1e-4,0, 0,
           -1,1e-5, 1e-4,
           1e-4, 0, 
           1e-4,1e-4,
           1e-4,1e-5,1e-4, -10,
           -pi/2,1e-4)
  
  uppar=c(30,  3, 100,
          50, 3 ,1,
          1, 3,  150  ,
          50,  2,
          200,2,
          2,5,100,1,
          0,10)

  
  #init.par: initial parameters taken from final estimates given in paper
  init.par=c(1.6, 0.64, 0,
             32.49,  0.73, 1,
             0.7, 0.28, 47.56,
             20.87, 0.76,
             104.81,   0.41,
             0.77, 0.32, 57.9, 1,
             -0.15, 0.97)

  opt=optimParallel(par=init.par,
                    lower=lowpar,
                    upper=uppar,
                    fn=spatmodnll_conv,Exceed.all=Exceed.pairs,Exceed.Inds=Exceed.Inds,coords=coords,
                    c.vec=c.vec,parallel=list(loginfo=TRUE))
  print(opt)
  
  stopCluster(cl)
  
  spat_pars=opt$par
  save(spat_pars,file="DependenceAnalysis/convfullspatfit.Rdata")
}else if(mix.boo==T & conv.boo==F){
  #Please refer to spatmodnll_nonconv in spatial_fit_funcs.R for details of the full model and its parameters.
  
  lowpar=c(1e-4,  0, 0,
           0,1e-5, 1e-5,
           -2,1e-5, 1e-4,
           1e-4, 1e-4, 1e-4,
           1e-4,1e-4,
           1e-5,1e-5,1e-4, -10,
           -pi/2,0.1)
  
  uppar=c(300,  2, 80,
          1, 2.5 ,500,
          1, 2.5,  100  ,
          5000/100,  2, 3,
          25,2,
          0.1*100,3.5,150,1,
          0,3)
  
  
  #init.par: initial parameters taken from final estimates given in paper
  
  #Note that init.par[13] and init.par[15] are scaled inside the nll function to avoid numerical instability during optimisiation.
  #This is particularlly important for init.par[15], which can become quite small.
  init.par=c(170.07, 0.871, 6.806,
             0.125, 1.357, 255.997, 
             -0.055, 0.395,   15.503, 
             16.794, 0.725, 2.113,
             1496.357/100, 0.364, 
             0.002*100, 1.652,   90.57, 0.285, 
             -0.322, 0.997)
  
  opt=optimParallel(par=init.par,
                    lower=lowpar,
                    upper=uppar,
                    fn=spatmodnll_nonconv,Exceed.all=Exceed.pairs,Exceed.Inds=Exceed.Inds,coords=coords,
                    c.vec=c.vec,parallel=list(loginfo=TRUE))
  
  
  print(opt)
  
  stopCluster(cl)
  
  spat_pars=opt$par
  spat_pars[13]=spat_pars[13]*100
  spat_pars[15]=spat_pars[15]/100
  
  save(spat_pars,file="DependenceAnalysis/nonconvfullspatfit.Rdata")
}

