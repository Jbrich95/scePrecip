rm(list=ls())
library(evd)
library(ismev)
library(maps)
library(mapdata)
library(mapproj)
library(lattice)
library(evgam)
library(fields)
library(viridis)
library(qgam)
jobid<-commandArgs(trailing=T)
jobid<- as.numeric(jobid)
print(jobid)
load("All_Original_Data_Large_Masked4.Rdata")

ncols=120
col1 <- viridis(ncols)
contour=FALSE
region=TRUE
pretty=FALSE


Orig.coords=coords
Orig.Data=Data.Mat_all

rm(Data.Mat_all)
lat_range=range(Orig.coords[,2])+c(-0.2,0.2)
lon_range=range(Orig.coords[,1])+c(-0.1,0.1)
plotmap     <- map('worldHires', xlim=lon_range,ylim=lat_range,interior=F,plot=FALSE)


lat_scale <- pretty(as.vector(lat_range))
lon_scale <- pretty(as.vector(lon_range))

mappanel <- function(x,y,...) { # I think this is the fn that plots the map lines within the call to contour
  panel.contourplot(x,y,...)
  llines(plotmap$x,plotmap$y,col="black")
  #map.grid(plotmap)
}

load("Conv_alt.Rdata")
r=0.04
Data.sub=data.frame(cbind( 
  "lon" = Orig.coords[,1]),
  "lat"= Orig.coords[,2],
  "alt"=alt
)

alt_bins=seq(0,max(alt)+0.2,length=6)

if(jobid==1){
  q=0.98
  
}else if(jobid==2){
  q=0.9995
  
}else if(jobid==3){
  q=0.9925
  
}else if(jobid==4){
  q=0.99
}else if(jobid==5){
  q=0.999
  
}else if(jobid==6){
  q=0.995
  
}else if(jobid==7){
  q=0.9975
}
print(q)
n.sub=250
print(n.sub)
alt_bin_inds=vector("list",length(alt_bins)-1)

for(i in 1:length(alt_bin_inds)){
  alt_bin_inds[[i]]=which(alt>=alt_bins[i] & alt < alt_bins[i+1])
}

sub.inds=c()
while(length(sub.inds)<n.sub){
  which_bin=sample(1:length(alt_bin_inds),1)
  if(length(alt_bin_inds[[which_bin]])>0){
    loc_ind=sample(alt_bin_inds[[which_bin]],1)
    sub.inds=c(sub.inds,loc_ind)
    alt_bin_inds[[which_bin]]=alt_bin_inds[[which_bin]][-which(alt_bin_inds[[which_bin]]==loc_ind)]
  }
  
}


sub.inds=sort(sub.inds)
sub.data=Orig.Data[,sub.inds]

sub.coords=Orig.coords[sub.inds,]
plot(Orig.coords)



points(sub.coords,col="red",pch=20)
##Throw away some data and keep scaling.factor proportion

scaling.factor=0.4
pr<-lon.vec<-lat.vec<-alt.vec<-c()
for(i in 1:n.sub){
  temp=Orig.Data[,sub.inds[i]]
  temp=sort(temp,decreasing = T)[1:(length(temp)*scaling.factor)]
  
  
  pr=c(pr,temp)
  alt.vec=c(alt.vec,rep(alt[sub.inds[i]],length(temp)))
  lon.vec=c(lon.vec,rep(Orig.coords[sub.inds[i],1],length(temp)))
  lat.vec=c(lat.vec,rep(Orig.coords[sub.inds[i],2],length(temp)))
}

Data=data.frame(pr,lon.vec,lat.vec,alt.vec)
names(Data)=c("pr","lon","lat","alt")


threshfit<-qgam(pr ~ 1+s(lon, lat, k = 4)+s(alt, k=4),data=Data,qu=1-(1-q)/scaling.factor)


pred.quant=predict(threshfit,newdata=Data.sub)




jpeg(filename=paste0("Smooth_GAM/Smooth_qgam_quant_Bothq",q*10000,"_4.jpeg"),width=600,height=600)
title1 <- ''
c    <- 0
c    <- contourplot(as.matrix(pred.quant) ~ (ceiling(as.matrix(Orig.coords[,1])/r)*r)*(ceiling(as.matrix(Orig.coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(x=list(at=lon_scale), y=list(at=lat_scale) ),
                    panel=mappanel, aspect="iso", region=region, main= title1,
                    contour=contour, pretty=pretty, cuts=ncols-1,   col.regions=col1 )

print(c)
dev.off()


n.sub=4000
print(n.sub)
n.sub=dim(Orig.coords)[1]

sub.inds=sort(sample(1:dim(Orig.coords)[1],n.sub))

sub.data=Orig.Data[,sub.inds]

sub.coords=Orig.coords[sub.inds,]

no.exceed=rep(0,n.sub)
Exceeds<-lon.vec<-lat.vec<-alt.vec<-c()
for(i in 1:n.sub){
  temp=Orig.Data[,sub.inds[i]]
  temp=temp[temp>pred.quant[sub.inds[i]]]-pred.quant[sub.inds[i]]
  no.exceed[i]=length(temp)
  
  Exceeds=c(Exceeds,temp)
  lon.vec=c(lon.vec,rep(Orig.coords[sub.inds[i],1],no.exceed[i]))
  lat.vec=c(lat.vec,rep(Orig.coords[sub.inds[i],2],no.exceed[i]))
  alt.vec=c(alt.vec,rep(alt[sub.inds[i]],no.exceed[i]))
}

Data=data.frame(Exceeds,lon.vec,lat.vec,alt.vec)
names(Data)=c("Exceedances","lon","lat","alt")


#Fit GPD gam
fmla_gpd <- list(Exceedances ~ 1+s(lon, lat, k = 4)+s(alt,k=4),~1)


m_gpd <- evgam(fmla_gpd, Data, family = "gpd")
summary(m_gpd)
# plot(m_gpd)

#Plot
Data.sub=data.frame(cbind( 
  "lon" = Orig.coords[,1]),
  "lat"= Orig.coords[,2],
  "alt"=alt
)
gpd_pred=predict(m_gpd,Data.sub,type="response")
gpd_pred[,3]=pred.quant
colnames(gpd_pred)[3]="Threshold"

save(gpd_pred,file=paste0("Smooth_GAM/Both_qgam_GAM_GPDfitsq",q*10000,"_4.Rdata"))


jpeg(filename=paste0("Smooth_GAM/Smooth_qgam_Scales_Bothq",q*10000,"_4.jpeg"),width=600,height=600)
title1 <- ''
c    <- 0
c    <- contourplot(as.matrix(gpd_pred[,1]) ~ (ceiling(as.matrix(Orig.coords[,1])/r)*r)*(ceiling(as.matrix(Orig.coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(x=list(at=lon_scale), y=list(at=lat_scale) ),
                    panel=mappanel, aspect="iso", region=region, main= title1,
                    contour=contour, pretty=pretty, at=seq(min(gpd_pred[,1]),max(gpd_pred[,1]),length=ncols-1),   col.regions=col1 )

print(c)
dev.off()



jpeg(filename=paste0("Smooth_GAM/Smooth_qgam_Shapes_Bothq",q*10000,"_4.jpeg"),width=600,height=600)

title1 <- ''
c    <- 0
c    <- contourplot(as.matrix(gpd_pred[,2])  ~  (ceiling(as.matrix(Orig.coords[,1])/r)*r)*(ceiling(as.matrix(Orig.coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    scales=list(x=list(at=lon_scale), y=list(at=lat_scale) ),
                    panel=mappanel, aspect="iso", region=region, main= title1,
                    contour=contour, pretty=pretty, cuts=ncols-1,  col.regions=col1 )

print(c)    
dev.off()



GPDthresh=gpd_pred[,3]
GPDshape=gpd_pred[,2]
GPDscale=gpd_pred[,1]
unif=function(x){rank(x)/(length(x)+1)}

Orig.Data_U=Orig.Data
exceed.inds=vector("list",dim(Orig.Data)[2])

for(i in 1:dim(Orig.Data)[2]){
  U=GPDthresh[i]
  Orig.Data_U[,i]=unif(Orig.Data[,i])
  exceed.inds[[i]]=which(Orig.Data[,i]>U)
}

Exp.Data=c()
for(i in 1:dim(Orig.Data)[2]){
  U=GPDthresh[i]
  temp=Orig.Data[,i][exceed.inds[[i]]]
  Exp.Data=c(pgpd(temp,loc=U,scale=GPDscale[i],shape=GPDshape[i]),Exp.Data)
}

Exp.Data=qexp(Exp.Data)
probs=seq(0.0001,1-1/length(Exp.Data),length=1000)

qe=qexp(probs)
png(filename=paste0("Smooth_GAM/GAM_qgam_diag_pool_Bothq",q*10000,"_4.png"),width=800,height=600)

plot(qe,quantile(c(Exp.Data),probs), xlab="Fitted",main="Both",
     ylab= "",cex.main=1.75, cex.lab=1.75, cex.axis=1.75, xlim=c(0,15),ylim=c(0,15))
mtext(side=2, "Empirical",cex = 1.75,line=2.6)
abline(a=0,b=1,col="red",lwd=2)
points(qe,quantile(c(Exp.Data),probs))

dev.off()


for(return.period in c(1,10,20,50)){
  prob=1-1/(return.period*dim(Orig.Data)[1]/20)
  
  quants=matrix(ncol=dim(Orig.Data)[2],nrow=1)
  for(i in 1:length(quants)){
    quants[,i]=GPDthresh[i]+qgpd((prob-q)/(1-q),loc=GPDthresh[i],scale=GPDscale[i],shape=GPDshape[i])
    
  }
  
  png(filename=paste0("Smooth_GAM/RL_qgam_",return.period,"_Bothq",q*10000,"_4.png"),width=800,height=600)
  title1 <- ''
  c    <- 0
  c    <- contourplot(as.matrix(quants)  ~  (ceiling(as.matrix(Orig.coords[,1])/r)*r)*(ceiling(as.matrix(Orig.coords[,2])/r)*r),
                      ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                      scales=list(x=list(at=lon_scale), y=list(at=lat_scale) ),
                      panel=mappanel, aspect="iso", region=region, main= title1,
                      contour=contour, pretty=pretty, cuts=ncols-1,  col.regions=col1 )
  
  print(c)   
  dev.off()   
}


#Fit logistic

cstar=1e-5
sub.inds=sort(sample(1:dim(Orig.coords)[1],n.sub))
t.coord1.vec=Orig.coords[sub.inds,1]
t.coord2.vec=Orig.coords[sub.inds,2]
alt.vec=alt[sub.inds]
dat=matrix(NA,ncol=2,nrow=n.sub)
for(i in 1:n.sub){
  dat[i,2]=sum(Orig.Data[,sub.inds[i]]>cstar)
  dat[i,1]=sum(Orig.Data[,sub.inds[i]]<= cstar)
}
Data=data.frame(dat,t.coord1.vec,t.coord2.vec,alt.vec)
names(Data)=c("rain","norain","lon","lat","alt")
fit<-gam(cbind(rain,norain)~1+s(lat,lon,k=4)+s(alt,k=4),family=binomial,data=Data)

prob=exp(predict(fit,Data.sub))/(1+exp(predict(fit,Data.sub)))
jpeg(filename="Smooth_GAM/Smooth_probnorain_Both_4.jpg",width=600,height=600)
title1 <- ''
c    <- 0
c    <- contourplot(as.matrix(prob) ~ (ceiling(as.matrix(Orig.coords[,1])/r)*r)*(ceiling(as.matrix(Orig.coords[,2])/r)*r),
                    ylab="",xlab="", xlim=lon_range, ylim=lat_range,
                    contour=F,aspect="iso",region=region,
                    panel=mappanel, main= title1,
                    at=seq(min(prob),max(prob),length=ncols-1),   col.regions=col1 )

print(c)
dev.off()

save(prob,file="Smooth_GAM/SmoothProbNoRain_Both_4.Rdata")
