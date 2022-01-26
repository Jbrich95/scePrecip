#Import all required packages
source("RequiredPackages.R")

##Load required Rdata
load("Data/Data.Rdata")
load("MarginalAnalysis/ProbNoRain.Rdata")
load("MarginalAnalysis/GPDfits.Rdata")


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


# Spatial-censoring thresholds on Laplace margins
c.vec=qlaplace(prob_zero)


# Transform all data to uniform margins. GPD tails above exceedance thresholds gpd_pred[,3] and a point mass at zero given by the prob_zero vector. 
# For bulk of distribution, use empirical distribution of non-zero values ONLY.

unif=function(x){rank(x)/(length(x)+1)}

lambda=0.005 #Exceedance probability for estimated q

Data_U=Data #Transform to Uniform margins
for(i in 1:dim(Data)[2]){
  q=as.numeric(gpd_pred[i,3])
  p=prob_zero[i]
  temp=Data[,i]
  
  #Point mass for zero values
  zero.inds=which(temp==0)
  Data_U[zero.inds,i]=p
  
  exceed.inds=which(temp >q)
  pos_emp_dist=unif(temp[-zero.inds]) # Empirical distribution of non-zero values

  Data_U[-zero.inds,i]=(1-p-lambda)*pos_emp_dist/mean(Data[-zero.inds,i]<q)+p
  
  #Use GPD CDF above threshold q

  Data_U[exceed.inds,i]=1-lambda*pgpd(Data[exceed.inds,i],loc=q,scale=gpd_pred[i,1],shape=gpd_pred[i,2],lower.tail = F)
}

# Transform to Laplace margins
Dat_Lap=qlaplace(Data_U)

# Save with censoring threshold
save(Dat_Lap,c.vec,file="MarginalAnalysis/Laplace_Data.Rdata")



