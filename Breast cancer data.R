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
load('BRCA_data.RData')
lib<-data$lib
disp<-1/data$disp
data<-data$data

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
tune<-seq(0,1800,20)
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
tuning_param_gauss<-seq(0,175,5)
mc.cores=1
model_gauss<-mclapply(1:length(tuning_param_gauss),function(i){
  res<-sgClust(data=data.sd,c_center=center_gauss,lambda=tuning_param_gauss[i],K=3)
  return(res)
},mc.cores = mc.cores)
##############sparse kmeans
model.skmeans<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)
