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
sim.data<-Sim.Corr(alpha=0.5,sim_disp,empirical_dist)
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

##########snbClust
tune<-seq(0,12,0.75)
model_nb<-lapply(1:length(tune),function(i){
  res<-snbClust(data=sim.data,lib=est_lib,k=3,phi=disp,c_center=center,penalize=TRUE,tune=tune[i],max_iter_beta = 500)
  return(res)
})


###############get initial for sgClust
result<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)
center1<-apply(data.sd[,which(result[[20]]$Cs==1)],1,mean)
center2<-apply(data.sd[,which(result[[20]]$Cs==2)],1,mean)
center3<-apply(data.sd[,which(result[[20]]$Cs==3)],1,mean)
center_gauss<-cbind(center1,center2,center3)
tuning_param_gauss<-seq(0,4,0.5)
model_gauss<-lapply(1:length(tuning_param_gauss),function(i){
  res<-sgClust(data=data.sd,c_center=center_gauss,lambda=tuning_param_gauss[i],K=3)
  return(res)
})
##############sparse kmeans
model.skmeans<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)



