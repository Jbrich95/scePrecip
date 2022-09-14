# Import all required packages
source("RequiredPackages.R")

# Are you using the mixture model as defined in Richards et al. (2022b)? If so, set mix.boo==T. 
# Do you want the convective or non-convective fit? For the former, set conv.boo==T and for the latter, set conv.boo==F.
mix.boo <- T; conv.boo <-T #If using the model of Richards et al. (2022a), set mix.boo==F.

##Load required Rdata
if(mix.boo==F){ load("Data/Data.Rdata")
}else if(mix.boo==T & conv.boo == T){ load("Data/conv.Rdata")
}else if(mix.boo==T & conv.boo == F){ load("Data/nonconv.Rdata")}
##Inputs:
#  
# For n observed fields with d sampling locations
#
#  Data: n x d matrix of hourly precipitation rate (mm/hour). We set all values <= 1e-5 to 0.
#  coords: d x 2 matrix of lon/lat coordinates

# If mix.boo <-T, require
#  elev: d vector of elevation values, which will be using in all GAM formulae




##Create map for plots
ncols=120
col1 <- viridis(ncols)

lat_range=range(coords[,2])+c(-0.2,0.2)
lon_range=range(coords[,1])+c(-0.1,0.1)
plotmap     <- map('worldHires', xlim=lon_range,ylim=lat_range,interior=F,plot=FALSE)

mappanel <- function(x,y,...) {
  panel.contourplot(x,y,...)
  llines(plotmap$x,plotmap$y,col="black")
}

#Set exceedance probability. This is constant across all grid boxes.

lambda=0.005

#Evaluate site-wise quantiles
if(mix.boo==F){ # if mix.boo = F, fits TPS through sitewise point estimate of quantiles
  
  quants<-rep(0, nrow(coords))
for(i in 1:nrow(coords)) quants[i]=quantile(Data[,i],prob=1-lambda)


#Fit thin-plate spline to pointwise estimates
fit<-Tps(coords,quants)


#Predict over entire spatial domain
pred.quant=predict(fit)


}else if(mix.boo=T){ # if mix.boo = F, uses qgam package
  
  n.sub=500 # sub-sample n.sub locations
  
  # We sub-sample sites with probability approximately equal to that of the empirical elevation distribution.
  # Create 10 bins of elevation and draw sites from the bins 
  elev.bins=seq(0,max(elev)+0.2,length=10)
  elev.bin.inds=vector("list",length(elev.bins)-1)
  for(i in 1:length(elev.bin.inds)) elev.bin.inds[[i]]=which(elev>=elev.bins[i] & elev < elev.bins[i+1])
  
  
  sub.inds=c()
  while(length(sub.inds)<n.sub){
    which_bin=sample(1:length(elev.bin.inds),1)
    if(length(elev.bin.inds[[which_bin]])>0){
      loc_ind=sample(elev.bin.inds[[which_bin]],1)
      sub.inds=c(sub.inds,loc_ind)
      elev.bin.inds[[which_bin]]=elev.bin.inds[[which_bin]][-which(elev.bin.inds[[which_bin]]==loc_ind)]
    }
    
  }
  
  
  sub.inds=sort(sub.inds)
  sub.data=Data[,sub.inds]

  #Plot subset of coordinates
  plot(coords,xlab="",ylab="")
  points(coords[sub.inds,],col="red",pch=20)

  #Create data.frame
  pr<-lon.vec<-lat.vec<-elev.vec<-c()
  for(i in 1:n.sub){
    temp=Data[,sub.inds[i]]

    
    pr=c(pr,temp)
    elev.vec=c(elev.vec,rep(elev[sub.inds[i]],length(temp)))
    lon.vec=c(lon.vec,rep(coords[sub.inds[i],1],length(temp)))
    lat.vec=c(lat.vec,rep(coords[sub.inds[i],2],length(temp)))
  }
  
  df=data.frame(pr,lon.vec,lat.vec,elev.vec)
  names(df)=c("pr","lon","lat","elev")
  
  # Fit quantile regression GAM to estimate 1-lambda quantile
  fit<-qgam(pr ~ 1+s(lon, lat, k = 4)+s(elev, k=4),data=df,qu=1-lambda)
  
  df=data.frame(cbind( 
    "lon" = coords[,1]),
    "lat"= coords[,2],
    "elev"=elev
  )
  
  #Predict exceedance thresholds
  pred.quant=predict(fit,newdata=df)
  
  
}

#Plot spatial map of threshold

r=0.073 # r changes the resolution of map. Each coordinate is rounded to the nearest r for plotting

c    <- 0
c    <- contourplot(as.matrix(pred.quant) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    panel=mappanel, aspect=1, region=T, main= "",  scales=list(tck=c(0,0),draw=F),
                    contour=F, pretty=F, cuts=ncols-1,   col.regions=col1 )

print(c)


## GAM fits for GPD tails 

#Put into dataframe
no.exceed=rep(0,nrow(coords))
Exceedances<-lon.vec<-lat.vec<-elev.vec<-c()
for(i in 1:nrow(coords)){
  temp=Data[,i]
  temp=temp[temp>pred.quant[i]]-pred.quant[i]
  no.exceed[i]=length(temp)
  
  Exceedances=c(Exceedances,temp)
  lon.vec=c(lon.vec,rep(coords[i,1],no.exceed[i]))
  lat.vec=c(lat.vec,rep(coords[i,2],no.exceed[i]))
  
  if(mix.boo==T)   elev.vec=c(elev.vec,rep(elev[i],no.exceed[i]))

}

if(mix.boo==F){  
df=data.frame(Exceedances,lon.vec,lat.vec)   
names(df)=c("Exceedances","lon","lat")
}else if(mix.boo==T){
  df=data.frame(Exceedances,lon.vec,lat.vec,elev.vec)   
  names(df)=c("Exceedances","lon","lat","elev")
  
}


#Fit GPD GAM to exceedances. 

if(mix.boo==F){
#Here we use k = 10 knots to speed up computation. In the paper, we set k = 300, but this takes some time to run.
fmla_gpd <- list(Exceedances ~ 1+s(lon, lat, k = 10),~1+s(lon,lat,k=10))
}else if(mix.boo==T){
  fmla_gpd <- list(Exceedances ~ 1+s(lon, lat, k = 4)+s(elev,k=4),~1)  #Note that shape is fixed over domain
  
}
gpd_gam_fit <- evgam(fmla_gpd, df, family = "gpd")
summary(gpd_gam_fit)

#Predict over entire domain
if(mix.boo==F){
gpd_pred=predict(gpd_gam_fit,newdata = data.frame("lon"=coords[,1],"lat"=coords[,2]),type="response")
}else if(mix.boo==T){
  gpd_pred=predict(gpd_gam_fit,newdata = data.frame("lon"=coords[,1],"lat"=coords[,2],"elev"=elev),type="response")
  
}
#Add estimated threshold to prediction matrix
gpd_pred[,3]=pred.quant
colnames(gpd_pred)[3]="Threshold"

if(mix.boo==F) save(gpd_pred,file=paste0("MarginalAnalysis/GPDfits.Rdata"))
if(mix.boo==T & conv.boo==T) save(gpd_pred,file=paste0("MarginalAnalysis/convGPDfits.Rdata"))
if(mix.boo==T & conv.boo==F) save(gpd_pred,file=paste0("MarginalAnalysis/nonconvGPDfits.Rdata"))


#Plot scale
c    <- 0
c    <- contourplot(as.matrix(gpd_pred[,1]) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(tck=c(0,0),draw=F),                
                    panel=mappanel, aspect=1, region=T, main= "",
                    contour=F, pretty=F, at=seq(min(gpd_pred[,1]),max(gpd_pred[,1]),length=ncols-1),   col.regions=col1 )

print(c)


#Plot shape if mix.boo==F. Else, the shape is fixed over the domain
if(mix.boo==F){
c    <- 0
c    <- contourplot(as.matrix(gpd_pred[,2]) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(tck=c(0,0),draw=F),            
                    panel=mappanel, aspect=1, region=T, main= "",
                    contour=F, pretty=F, at=seq(min(gpd_pred[,2]),max(gpd_pred[,2]),length=ncols-1),   col.regions=col1 )

print(c)
}

GPDthresh=gpd_pred[,3]
GPDshape=gpd_pred[,2]
GPDscale=gpd_pred[,1]

#Plot 1,10,20,50 year return level estimates

for(return.period in c(1,10,20,50)){
  prob=1-1/(return.period*dim(Data)[1]/20)
  
  quants=matrix(ncol=dim(coords)[1],nrow=1)
  for(i in 1:length(quants)){
    quants[,i]=GPDthresh[i]+qgpd((prob-(1-lambda))/lambda,loc=GPDthresh[i],scale=GPDscale[i],shape=GPDshape[i])
    
  }

  title=paste0(return.period, " year return-level estimate")
  c    <- 0
  c    <- contourplot(as.matrix(quants)  ~  (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                      ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                      scales=list(tck=c(0,0),draw=F),           
                      panel=mappanel, aspect=1, region=T, main= title,
                      contour=F, pretty=F, cuts=ncols-1,  col.regions=col1 )
  
  print(c)   
}


##We now fit the logistic GAM to estimate the probability that there is no rain at a site. Can use all sites for fitting.

#Calculate number of times with rain/norain at each site
dat=matrix(NA,ncol=2,nrow=dim(coords)[1])
for(i in 1:dim(coords)[1]){
  dat[i,2]=sum(Data[,i]>0)
  dat[i,1]=sum(Data[,i]==0)
}

if(mix.boo==F){
df_prob=data.frame(dat,coords[,1],coords[,2])
names(df_prob)=c("rain","norain","lon","lat")
}else if(mix.boo==T){
  df_prob=data.frame(dat,coords[,1],coords[,2],elev)
  names(df_prob)=c("rain","norain","lon","lat","elev")
  
}
#Fit model
if(mix.boo==F){
#As in the paper, we use k = 300 knots
logistic_fit<-gam(cbind(rain,norain)~1+s(lat,lon,k=300),family=binomial,data=df_prob)
}else if(mix.boo==T){
  #As in the paper, we use only k = 10 knots
  logistic_fit<-gam(cbind(rain,norain)~1+s(lat,lon,k=4)+s(elev,k=4),family=binomial,data=df_prob)
  
}
#Predictions
if(mix.boo==F){
prob_zero=exp(predict(logistic_fit, data.frame("lon"=coords[,1],"lat"=coords[,2])))/
  (1+exp(predict(logistic_fit, data.frame("lon"=coords[,1],"lat"=coords[,2]))))
}else if(mix.boo==T){
  prob_zero=exp(predict(logistic_fit, data.frame("lon"=coords[,1],"lat"=coords[,2],"elev"=elev)))/
    (1+exp(predict(logistic_fit, data.frame("lon"=coords[,1],"lat"=coords[,2],"elev"=elev))))
}

c    <- 0
c    <- contourplot(as.matrix(prob_zero) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(tck=c(0,0),draw=F),           
                    panel=mappanel, aspect=1, region=T, main= "",
                    contour=F, pretty=F, cuts=ncols-1,  col.regions=col1 )

print(c)

#Save estimated probabilities

if(mix.boo==F) save(prob_zero,file="MarginalAnalysis/ProbNoRain.Rdata")
if(mix.boo==T & conv.boo==T) save(prob_zero,file="MarginalAnalysis/convProbNoRain.Rdata")
if(mix.boo==T & conv.boo==F) save(prob_zero,file="MarginalAnalysis/nonconvProbNoRain.Rdata")

