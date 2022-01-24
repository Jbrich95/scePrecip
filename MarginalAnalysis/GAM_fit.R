
#Import all required packages
source("RequiredPackages.R")


##Objects:
#  
# For n observed fields with d sampling locations
#
#  Data: n x d matrix of hourly precipitation rate (mm/hour). We set all values <= 1e-5 to 0.
#  coords: d x 2 matrix of lon/lat coordinates



##Create map for plots
ncols=120
col1 <- viridis(ncols)
contour=FALSE
region=TRUE
pretty=FALSE


lat_range=range(coords[,2])+c(-0.2,0.2)
lon_range=range(coords[,1])+c(-0.1,0.1)
plotmap     <- map('worldHires', xlim=lon_range,ylim=lat_range,interior=F,plot=FALSE)

lat_scale <- pretty(as.vector(lat_range))
lon_scale <- pretty(as.vector(lon_range))

mappanel <- function(x,y,...) {
  panel.contourplot(x,y,...)
  llines(plotmap$x,plotmap$y,col="black")
}

#Set exceedance probability. This is constant across all grid boxes.

lambda=0.005

#Evaluate site-wise quantiles

quants<-rep(0, nrow(coords))
for(i in 1:nrow(coords)) quants[i]=quantile(Data[,i],prob=1-lambda)


#Fit thin-plate spline to pointwise estimates
fit<-Tps(coords,quants)


#Predict over entire spatial domain
pred.quant=predict(fit)


#Plot spatial map of threshold

r=0.073 # r changes the resolution of map. Each coordinate is rounded to the nearest r for plotting

c    <- 0
c    <- contourplot(as.matrix(pred.quant) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    panel=mappanel, aspect="iso", region=region, main= "",  scales=list(tck=c(0,0),draw=F),
                    contour=contour, pretty=pretty, cuts=ncols-1,   col.regions=col1 )

print(c)


##EVGAM fits for GPD tails 

#Put into dataframe
no.exceed=rep(0,nrow(coords))
Exceedances<-lon.vec<-lat.vec<-c()
for(i in 1:nrow(coords)){
  temp=Data[,i]
  temp=temp[temp>pred.quant[i]]-pred.quant[i]
  no.exceed[i]=length(temp)
  
  Exceedances=c(Exceedances,temp)
  lon.vec=c(lon.vec,rep(coords[i,1],no.exceed[i]))
  lat.vec=c(lat.vec,rep(coords[i,2],no.exceed[i]))
}

df=data.frame(Exceedances,lon.vec,lat.vec)
names(df)=c("Exceedances","lon","lat")


#Fit GPD GAM to exceedances. Note that shape parameter is fixed over all sites.
#Here we use k = 10 knots to speed up computation. In the paper, we set k = 300, but this takes some time to run.

fmla_gpd <- list(Exceedances ~ 1+s(lon, lat, k = 10),~1+s(lon,lat,k=10))

gpd_gam_fit <- evgam(fmla_gpd, df, family = "gpd")
summary(gpd_gam_fit)

#Predict over entire domain
gpd_pred=predict(gpd_gam_fit,newdata = data.frame("lon"=coords[,1],"lat"=coords[,2]),type="response")

#Add estimated threshold to prediction matrix
gpd_pred[,3]=pred.quant
colnames(gpd_pred)[3]="Threshold"

save(gpd_pred,file=paste0("MarginalAnalysis/GPDfits.Rdata"))


#Plot scale
c    <- 0
c    <- contourplot(as.matrix(gpd_pred[,1]) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(tck=c(0,0),draw=F),                    panel=mappanel, aspect="iso", region=region, main= "",
                    contour=contour, pretty=pretty, at=seq(min(gpd_pred[,1]),max(gpd_pred[,1]),length=ncols-1),   col.regions=col1 )

print(c)


#Plot shape
c    <- 0
c    <- contourplot(as.matrix(gpd_pred[,2]) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(tck=c(0,0),draw=F),                    panel=mappanel, aspect="iso", region=region, main= "",
                    contour=contour, pretty=pretty, at=seq(min(gpd_pred[,2]),max(gpd_pred[,2]),length=ncols-1),   col.regions=col1 )

print(c)

GPDthresh=gpd_pred[,3]
GPDshape=gpd_pred[,2]
GPDscale=gpd_pred[,1]

#Plot 1,10,20,50 year return level estimates

for(return.period in c(1,10,20,50)){
  prob=1-1/(return.period*dim(Data)[1]/20)
  
  quants=matrix(ncol=dim(coords)[1],nrow=1)
  for(i in 1:length(quants)){
    quants[,i]=GPDthresh[i]+qgpd((prob-lambda)/(1-lambda),loc=GPDthresh[i],scale=GPDscale[i],shape=GPDshape[i])
    
  }

  title=paste0(return.period, " year return-level estimate")
  c    <- 0
  c    <- contourplot(as.matrix(quants)  ~  (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                      ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                      scales=list(tck=c(0,0),draw=F),           
                      panel=mappanel, aspect="iso", region=region, main= title,
                      contour=contour, pretty=pretty, cuts=ncols-1,  col.regions=col1 )
  
  print(c)   
}


##We now fit the logistic GAM to estimate the probability that there is no rain at a site. Can use all sites for fitting.

#Calculate number of times with rain/norain at each site
dat=matrix(NA,ncol=2,nrow=dim(coords)[1])
for(i in 1:dim(coords)[1]){
  dat[i,2]=sum(Data[,i]>0)
  dat[i,1]=sum(Data[,i]==0)
}

df_prob=data.frame(dat,coords[,1],coords[,2])
names(df_prob)=c("rain","norain","lon","lat")

#Fit  
#As in the paper, we use k = 300 knots
logistic_fit<-gam(cbind(rain,norain)~1+s(lat,lon,k=300),family=binomial,data=df_prob)

#Predictions
prob_zero=exp(predict(logistic_fit, data.frame("lon"=coords[,1],"lat"=coords[,2])))/(1+exp(predict(logistic_fit, data.frame("lon"=coords[,1],"lat"=coords[,2]))))


c    <- 0
c    <- contourplot(as.matrix(prob) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(tck=c(0,0),draw=F),           
                    panel=mappanel, aspect="iso", region=region, main= "",
                    contour=contour, pretty=pretty, cuts=ncols-1,  col.regions=col1 )

print(c)

#Save estimated probabilities
save(prob_zero,file="MarginalAnalysis/ProbNoRain.Rdata")
