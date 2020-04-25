# library(truncnorm)
# library(gplots)
# library(sparcl)
# library(mclust)
# library(edgeR)
# library(mvtnorm)
# library(MCMCpack)
##make sure all packages above has been imported
#####################
##The following shows how to get result for one run of simulation.
######################
library(snbClust)
load('BRCA_data.RData')
empirical_dist<-apply(data$data,1,mean)
sim_disp<-data$disp
####eff_a is the effect size to tune
sim.data<-Sim.Independent(ngenes=150,eff_a=1,percent_DE=1,sim_disp,empirical_dist)
disp<-1/estimateDisp(sim.data)$tagwise.dispersion
est_lib<-calcNormFactors(sim.data)

############get initial for snbClust
y1<-cpm(sim.data,prior.count=0.25,log=T)
data.sd<-t(apply(y1,1,function(x) (x-mean(x))))
result<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)

#########result is also the result for sparse kmeans
center1<-apply(sim.data[,which(result[[20]]$Cs==1)],1,mean)
center2<-apply(sim.data[,which(result[[20]]$Cs==2)],1,mean)
center3<-apply(sim.data[,which(result[[20]]$Cs==3)],1,mean)

center<-cbind(center1,center2,center3)

center[which(center==0)]<-0.1
tune<-0
##########snbClust
model_nb<-snbClust(data=sim.data,lib=est_lib,k=3,phi=disp,c_center=center,penalize=TRUE,tune=tune,max_iter_beta = 500)

###############get initial for sgClust
result<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)
center1<-apply(data.sd[,which(result[[20]]$Cs==1)],1,mean)
center2<-apply(data.sd[,which(result[[20]]$Cs==2)],1,mean)
center3<-apply(data.sd[,which(result[[20]]$Cs==3)],1,mean)
center_gauss<-cbind(center1,center2,center3)

model_gauss<-sgClust(data=data.sd,c_center=center_gauss,lambda=0,K=3)
##############sparse kmeans
model.skmeans<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)



