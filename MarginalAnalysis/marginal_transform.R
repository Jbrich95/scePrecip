#Import all required packages
source("RequiredPackages.R")


##Objects:
#  
# For n observed fields with d sampling locations
#
#  Data: n x d matrix of hourly precipitation rate (mm/hour). We set all values <= 1e-5 to 0.
#  coords: d x 2 matrix of lon/lat coordinates

##Load GAM fits
load("MarginalAnalysis/ProbNoRain.Rdata")
load("MarginalAnalysis/GPDfits.Rdata")



##Objects:
#  
# For d sampling locations
#  gpd_pred: a d x 3 vector of marginal parameter estimates. Each row corresponds to estimates at one site.
#  gpd_pred[,1]: d-vector of estimated GPD scale parameters
#  gpd_pred[,2]: d-vector of estimated GPD shape parameters
#  gpd_pred[,3]: d-vector of estimated GPD exceedance threshold
#  prob_zero: d-vector of estimated probability of no rain at each site


# Spatial-censoring thresholds on Laplace margins
c.vec=qlaplace(prob_zero)


# Transform all data to uniform margins. GPD tails above exceedance thresholds and a point mass at zero given by the prob vector. 
# For bulk of distribution, use empirical distribution of non-zero values ONLY.

unif=function(x){rank(x)/(length(x)+1)}

Data_U=Data
for(i in 1:dim(Data)[2]){
  U=as.numeric(gpd_pred[i,3])
  mass=prob_zero[i]
  temp=Data[,i]
  
  #Point mass for zero values
  zero.inds=which(temp==0)
  Data_U[zero.inds,i]=mass
  
  exceed.inds=which(temp >U)
  pos_emp_dist=unif(temp[-zero.inds]) # Empirical distribution of non-zero values

  Data_U[-zero.inds,i]=(1-mass)*pos_emp_dist+mass
  
  #Use GPD CDF above threshold U
  Data_U[exceed.inds,i]=(length(Data_U[,i])-length(exceed.inds))/(length(Data_U[,i]))+
    (length(exceed.inds))/(length(Data_U[,i]))*pgpd(Data[exceed.inds,i],loc=U,scale=gpd_pred[i,1],shape=gpd_pred[i,2])
}

# Transform to Laplace margins
Dat_Lap=qlaplace(Data_U)

# Save with censoring threshold
save(Dat_Lap,c.vec,file="MarginalAnalysis/Laplace_Data.Rdata")



