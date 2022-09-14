# Applies the convection identified algorithm described in Section 2.1 of Richards et al. (2022b)

#set hyperparameters
g_u = 1.0 # (mm/h) upper threshold
g_l = 0.01 # (mm/h) lower threshold
n_g = 9 # size of neighbourhood
p.star = 0.2 #proportion threshold

lon=lonlat.grid[,,1]
lat=lonlat.grid[,,2]
nlons=dim(Data.grid)[2]
nlats=dim(Data.grid)[1]

conv.precip=array(NA,dim(Data.grid))

# set gradhigh to 1 if high spatial gradient
# set gradlow to 1 if low or high (but non-zero) spatial gradient

for(t in 1:dim(Data.grid)[3]){
  gradlow <- gradhigh <- array(0,dim=c(nlats+1,nlons+1))
  gradlow[nlats+1,nlons+1]=NA
  gradhigh[nlats+1,nlons+1]=NA
  tempDat=Data.grid[,,t]
  for(i in 1:(nlats-1)){
    for(j in 1:nlons){  
      if(!is.na(tempDat[i+1,j]) & !is.na(tempDat[i,j])){
        if( tempDat[i+1,j] > tempDat[i,j]+g_u){
          gradhigh[i,j]=1.0
          gradhigh[i+1,j]=1.0
        }
        
        if(tempDat[i+1,j] > tempDat[i,j]+g_l){
          gradlow[i,j]=1.0
          gradlow[i+1,j]=1.0  
        }
        
        
        
      }
    }
    
    
  }
  
  for(i in 2:nlats){
    for(j in 1:nlons){    
      if(!is.na(tempDat[i-1,j])&!is.na(tempDat[i,j])){
        if( tempDat[i-1,j] > tempDat[i,j]+g_u){
          gradhigh[i,j]=1.0
          gradhigh[i-1,j]=1.0}
        
        if(tempDat[i-1,j] > tempDat[i,j]+g_l){
          gradlow[i,j]=1.0
          gradlow[i-1,j]=1.0  
        }
      }
    }
  }
  
  
  for(i in 1:nlats){
    for(j in 1:(nlons-1)){   
      if(!is.na(tempDat[i,j+1])&!is.na(tempDat[i,j])){
        if( tempDat[i,j+1] > tempDat[i,j]+g_u){
          gradhigh[i,j]=1.0
          gradhigh[i,j+1]=1.0}
        
        if(tempDat[i,j+1] > tempDat[i,j]+g_l){
          gradlow[i,j]=1.0
          gradlow[i,j+1]=1.0  
        }
      }
    }
  }
  for(i in 1:nlats){
    for(j in 2:nlons){      
      if(!is.na(tempDat[i,j-1])&!is.na(tempDat[i,j])){
        if( tempDat[i,j-1] > tempDat[i,j]+g_u){
          gradhigh[i,j]=1.0
          gradhigh[i,j-1]=1.0}
        
        if(tempDat[i,j-1] > tempDat[i,j]+g_l){
          gradlow[i,j]=1.0
          gradlow[i,j-1]=1.0  
        }
        
      }
    }}
  
  
  nlen = (n_g-1)/2
  for(i in (nlen+1):(nlats-nlen)){
    for(j in (nlen+1):(nlons-nlen)){
      
      ghdata = gradhigh[(i-nlen):(i+nlen),(j-nlen):(j+nlen)]
      gldata = gradlow[(i-nlen):(i+nlen),(j-nlen):(j+nlen)]
      if(sum(gldata) >= 0.5){
        nfrac = sum(ghdata)/sum(gldata)
        if(nfrac >= p.star){
          conv.precip[i,j,t] = Data.grid[i,j,t]
        }
      }
      if(any(is.na(lon[(i-nlen):(i+nlen),(j-nlen):(j+nlen)]))){
        conv.precip[i,j,t]=NA
      }
      
      if( sum(ghdata) > sum(gldata)){
        print("ERROR in totals")
        
      }
    }
  }
  
  if(t%%100==0){print(t)}
}


conv.precip=conv.precip[(nlen+1):(nlats-nlen),(nlen+1):(nlons-nlen),] #Remove edges
lat=lat[(nlen+1):(nlats-nlen),(nlen+1):(nlons-nlen)] 
lon=lon[(nlen+1):(nlats-nlen),(nlen+1):(nlons-nlen)]


#Re-shape to matrices
tmp<-matrix(ncol=prod(dim(conv.precip)[1:2]),nrow=dim(conv.precip)[3])
tmp2<-matrix(nrow=prod(dim(lat)[1:2]),ncol=2)
int<-1
for(i in 1:dim(conv.precip)[1]){
  for(j in 1:dim(conv.precip)[2]){
    tmp[,int]=conv.precip[i,j,]
    tmp2[int,1]=lon[i,j]; tmp2[int,2]=lat[i,j]
    int=int+1
  }
}

coords=tmp2

# Note that tmp currently includes only data identified as convective, with NAs otherwise. 
# We take the whole observed field as convective if any convective values are found within a field

conv.t.inds<-which(apply(tmp,1,sum,na.rm=T)>0)

Data<-tmp[conv.t.inds,]
save(Data,coords,file="Data/conv.Rdata")


Data<-tmp[-conv.t.inds,]
save(Data,coords,file="Data/nonconv.Rdata")