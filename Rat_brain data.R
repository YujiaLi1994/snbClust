# library(truncnorm)
# library(gplots)
# library(sparcl)
# library(mclust)
# library(edgeR)
# library(mvtnorm)
# library(MCMCpack)
##make sure all packages above has been imported
#####################
library(snbClust)
load('dat_mice_ctrl.RData')
data<-dat_mice_ctrl$m_ctrl
lib<-dat_mice_ctrl$lib
disp<-1/dat_mice_ctrl$disp
set.seed(2000)


###########snbClust
y1<-cpm(data,prior.count=0.25,log=T)
data.sd<-t(apply(y1,1,function(x) (x-mean(x))))
result<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)
center1<-apply(data[,which(result[[20]]$Cs==1)],1,mean)
center2<-apply(data[,which(result[[20]]$Cs==2)],1,mean)
center3<-apply(data[,which(result[[20]]$Cs==3)],1,mean)
center<-cbind(center1,center2,center3)
center[which(center==0)]<-0.1
tune<-c(seq(0,160,1),seq(160,240,0.5))
mc.cores=1####parallel computing
model_nb<-mclapply(1:length(tune),function(i){
  res<-snbClust(data=data,lib=lib,k=3,phi=disp,c_center=center,penalize=TRUE,tune=tune[i],max_iter_beta = 500)
  return(res)
},mc.cores=mc.cores)
#############sgClust

result<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)
center1<-apply(data.sd[,which(result[[20]]$Cs==1)],1,mean)
center2<-apply(data.sd[,which(result[[20]]$Cs==2)],1,mean)
center3<-apply(data.sd[,which(result[[20]]$Cs==3)],1,mean)
center_gauss<-cbind(center1,center2,center3)
tuning_param_gauss<-c(seq(0,16.44,0.1),seq(16.44,16.52,0.01))
mc.cores=1
model_gauss<-mclapply(1:length(tuning_param_gauss),function(i){
  res<-sgClust(data=data.sd,c_center=center_gauss,lambda=tuning_param_gauss[i],K=3)
  return(res)
},mc.cores = mc.cores)
##############sparse kmeans
model.skmeans<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)
