#Import all required packages
source("RequiredPackages.R")

##Load required Rdata
load("Data/Data.Rdata")
load("AggregateAnalysis/R_A.Rdata")
##Loaded Objects:
#  #  
# For n observed fields with d sampling locations
#
#  Data: n x d matrix of hourly precipitation rate (mm/hour). We set all values <= 1e-5 to 0.
#  coords: d x 2 matrix of lon/lat coordinates
#  R_A: vector of simulated R_\mathcal{A} for a single region A
#  agg.inds: vector giving the indices of rows in coords that contain the coordinates within \mathcal{A}


#Calculate empirical R_A
R_emp=rowMeans(Data[,agg.inds])

#Produce Q-Q plot to assess fit - Figure 4 in paper
p.max=1-1/dim(Data)[1]
probs=seq(0.8,p.max,length=5000)
qrdata=quantile(R_emp,prob=probs)
qrsim=quantile(R_A,prob=probs)
plot(qrdata,qrsim,ylim=range(qrdata,qrsim),xlim=range(qrdata,qrsim),main="",xlab="",ylab="")
mtext(side=2, "Model",cex = 1.5, line=2.6)
mtext(side=1, "Empirical",cex = 1.5, line=2.6)
abline(a=0,b=1,col="red")




#Estimate return level curves for R_A as the average over {Y(s):s \in \mathcal{A}} (not sum as in paper)

return.period=seq(-0.3,4,length=1000)
return.period=10^return.period

nyears=20 #Set number of years of data
prob=1-1/return.period/(nyears*dim(Data)[1])

# Fit GPDs to estimate return level curves

q_gpd=0.95 #exceedance threshold for GPD fit
fit_sim=gpd.fit(R_A,threshold = quantile(R_A,q_gpd)) #This fit uses the simulated R_A
fit_emp=gpd.fit(R_emp,threshold = quantile(R_emp,q_gpd)) #This uses the empirical R_A

RL_sim=qgpd((prob-q_gpd)/(1-q_gpd),loc=quantile(R_A,q_gpd),scale=fit_sim$mle[1],shape=fit_sim$mle[2])
RL_emp=qgpd((prob-q_gpd)/(1-q_gpd),loc=quantile(R_emp,q_gpd),scale=fit_emp$mle[1],shape=fit_emp$mle[2])

plot(return.period,RL_emp,col=2,lwd=3,ylim=range(0,RL_emp,RL_sim),type="l",xlab="Return-Period (years)",ylab="Return-Level (mm/hour)",log="x")
points(return.period,RL_sim,col=1,lwd=3,type="l")


