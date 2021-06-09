
#Import all required packages
source("RequiredPackages.R")

#Load Data, coordinates and elevation
load("Data/Data.Rdata")
load("Data/Elevation.Rdata")
load("Data/Coordinates.Rdata")

##Objects:
#  
# For n observed fields with d sampling locations
#
#  Data: n x d matrix of hourly precipitation rate (mm/hour). We have set all values <= 1e-5 to 0.
#  coords: d x 2 matrix of lon/lat coordinates
#  elev: d vector of elevation (m)


#Put everything in data.frame, required for predictions
df=data.frame(cbind( 
  "lon" = coords[,1]),
  "lat"= coords[,2],
  "elev"=elev
)

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

#Set exceedance probability. This is constant across all girdboxes.

lambda=0.995

#We use a subet of sites for estimating the exceedance threshold. These are sampled uniformly at random.
n.sub_thresh=500

sub.inds=sample(1:dim(coords)[1],n.sub_thresh)
sub.data_thresh=Data[,sub.inds]

sub.coords_thresh=coords[sub.inds,]
plot(coords,xlab="",ylab="")
points(sub.coords_thresh, col = "red",pch=20,cex=1.5)


##Data and covariates must be in a data frame to use qgam

pr<-lon.vec<-lat.vec<-elev.vec<-c()
n=dim(Data)[1]
for(i in 1:n.sub_thresh){
  
  pr=c(pr,sub.data_thresh[,i])
  elev.vec=c(elev.vec,rep(elev[sub.inds[i]],n))
  lon.vec=c(lon.vec,rep(sub.coords_thresh[i,1],n))
  lat.vec=c(lat.vec,rep(sub.coords_thresh[i,2],n))
}

df_sub=data.frame(pr,lon.vec,lat.vec,elev.vec)
names(df_sub)=c("pr","lon","lat","elev")

#Fit qgam model for exceedance threshold. We use k =4 knots for both splines.

##WARNING: Takes a long time to run!

#qgam_fit<-qgam(pr ~ 1+s(lon, lat, k = 4)+s(elev, k=4),data=df_sub,qu=lambda)

#Predict over entire spatial domain
pred.quant=predict(qgam_fit,newdata=df)


#Plot spatial map of threshold

r=0.03 # r changes the resolution of map. Each coordinate is rounded to the nearst r for plotting

c    <- 0
c    <- contourplot(as.matrix(pred.quant) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(x=list(at=lon_scale), y=list(at=lat_scale) ),
                    panel=mappanel, aspect="iso", region=region, main= "",
                    contour=contour, pretty=pretty, cuts=ncols-1,   col.regions=col1 )

print(c)


##EVGAM fits for GPD tails 
# Computationally, this is a lot easier than qgam fit, so using all sites is feasible. We use a subet here.

n.sub_ev=4000

sub.inds=sort(sample(1:dim(coords)[1],n.sub_ev))

sub.data_ev=Data[,sub.inds]

sub.coords_ev=coords[sub.inds,]

#Put into dataframe
no.exceed=rep(0,n.sub_ev)
Exceedances<-lon.vec<-lat.vec<-elev.vec<-c()
for(i in 1:n.sub_ev){
  temp=sub.data_ev[,i]
  temp=temp[temp>pred.quant[sub.inds[i]]]-pred.quant[sub.inds[i]]
  no.exceed[i]=length(temp)
  
  Exceedances=c(Exceedances,temp)
  lon.vec=c(lon.vec,rep(sub.coords_ev[i,1],no.exceed[i]))
  lat.vec=c(lat.vec,rep(sub.coords_ev[i,2],no.exceed[i]))
  elev.vec=c(elev.vec,rep(elev[sub.inds[i]],no.exceed[i]))
}

df_sub=data.frame(Exceedances,lon.vec,lat.vec,elev.vec)
names(df_sub)=c("Exceedances","lon","lat","elev")


#Fit GPD GAM. Again, k = 4 knots. Note that shape parameter is fixed over all sites.

##WARNING: Takes a long time to run!#

fmla_gpd <- list(Exceedances ~ 1+s(lon, lat, k = 4)+s(elev,k=4),~1)

#gpd_gam_fit <- evgam(fmla_gpd, df_sub, family = "gpd")
summary(gpd_gam_fit)

#Predict over entire domain
gpd_pred=predict(gpd_gam_fit,df,type="response")

#Add estimated threshold to prediction matrix
gpd_pred[,3]=pred.quant
colnames(gpd_pred)[3]="Threshold"

save(gpd_pred,file=paste0("MarginalAnalysis/GPDfits.Rdata"))


#Plot scale
c    <- 0
c    <- contourplot(as.matrix(gpd_pred[,1]) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(x=list(at=lon_scale), y=list(at=lat_scale) ),
                    panel=mappanel, aspect="iso", region=region, main= "",
                    contour=contour, pretty=pretty, at=seq(min(gpd_pred[,1]),max(gpd_pred[,1]),length=ncols-1),   col.regions=col1 )

print(c)


#Shape estimate
print(paste0("Shape esimate: ",round(gpd_pred[1,2],3)))


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
                      scales=list(x=list(at=lon_scale), y=list(at=lat_scale) ),
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

df_prob=data.frame(dat,coords[,1],coords[,2],elev)
names(df_prob)=c("rain","norain","lon","lat","elev")

#Fit  

logistic_fit<-gam(cbind(rain,norain)~1+s(lat,lon,k=4)+s(elev,k=4),family=binomial,data=df_prob)

#Predictions
prob=exp(predict(logistic_fit,df))/(1+exp(predict(logistic_fit,df)))


c    <- 0
c    <- contourplot(as.matrix(prob) ~ (ceiling(as.matrix(coords[,1])/r)*r)*(ceiling(as.matrix(coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    contour=F,aspect="iso",region=region,
                    panel=mappanel, main= "",
                    at=seq(min(prob),max(prob),length=ncols-1),   col.regions=col1 )

print(c)

#Save estimated probabilities
save(prob,file="MarginalAnalysis/ProbNoRain.Rdata")
