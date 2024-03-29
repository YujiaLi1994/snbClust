# Project Title

The package implements sparse negative binomial mixture models (snbClust) proposed in the following paper: 
Li Y, Rahman T, Ma T, et al. A sparse negative binomial mixture model for clustering RNA-seq count data[J]. Biostatistics, 2021. (https://doi.org/10.1093/biostatistics/kxab025). Function of sparse Gaussian mixture models is also included in the package. 

Two real applications evaluated in the paper is not included in the package but uploaded in this github project (one rat brain application and one TCGA breast cancer application), along with the code the reproduce the result in this paper.
### Installing

For installation, we recommend to unzip the tar.gz file first and then use devtools::install() to install the package, which can make sure to install all the depends. Make sure R package truncnorm, gplots, sparcl, mclust, edgeR, mvtnorm, MCMCpack and Brobdingnag have all been properly imported.


## The following shows how to get result for one run of simulation.

## Simulate the data
```
library(snbClust)

load('BRCA_data.RData') ####Use real data to guide the simulatiom

empirical_dist<-apply(data$data,1,mean)

sim_disp<-data$disp

#eff_a is the effect size to tune, generally eff_a large than 
0.8 will be strongly enough to give good separation.

sim.data<-Sim.Independent(ngenes=1000,eff_a=1,percent_DE=0.15,sim_disp,empirical_dist) 
#####simulation data according to simulation 2 in the paper.

disp<-1/estimateDisp(sim.data)$tagwise.dispersion

est_lib<-calcNormFactors(sim.data)
```
## SnbClust
```
y1<-cpm(sim.data,prior.count=0.25,log=T)

data.sd<-t(apply(y1,1,function(x) (x-mean(x))))

result<-KMeansSparseCluster(t(data.sd),K=3,nstart=150)

center1<-apply(sim.data[,which(result[[20]]$Cs==1)],1,mean)

center2<-apply(sim.data[,which(result[[20]]$Cs==2)],1,mean)

center3<-apply(sim.data[,which(result[[20]]$Cs==3)],1,mean)

center<-cbind(center1,center2,center3)

center[which(center==0)]<-0.1

tune<-seq(0,12,0.75)

model_nb<-lapply(1:length(tune),function(i){
  res<-snbClust(data=sim.data,lib=est_lib,k=3,phi=disp,
  c_center=center,penalize=TRUE,tune=tune[i],max_iter_beta = 500)
  return(res)
})
```

## SgClust
```
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
```




