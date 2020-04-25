##' @title Simulation 1
##' @param alpha gene correlation parameter, larger alpha creates larger correlation
##' @param sim_disp gene dispersion from breast cancer data
##' @param empirical_dist average gene expression from breast cancer data
##' @return data with correlated genes for simulation 3.
##' @references Rahman T, Li Y, Ma T, et al. A sparse negative binomial mixture model for clustering RNA-seq count data[J]. arXiv preprint arXiv:1912.02399, 2019.
##' @export
##'
Sim.Corr<-function(alpha=0.5,sim_disp,empirical_dist){
  #empirical_dist<-apply(data$data,1,mean)
  quantile<-quantile(empirical_dist,0.70)
  empirical_dist<-empirical_dist[empirical_dist<quantile]
  #sim_disp<-data$disp
  sim_disp<-sim_disp[empirical_dist<quantile]
#empirical_dist<-empirical_dist[empirical_dist>500]
#ari_em<-c()
#ari_kmeans<-c()
#ari_da<-c()
#for(ari_ind in 1:20){

  min_cor<--0.1
  max_cor<-0.1
  ngenes=1000
  percent_DE=0.15
  select_ind<-sample(c(1:length(empirical_dist)),size = ngenes,replace = TRUE)
  
  #lfc=1
  #phi=1
  #sample1=15
  #sample2=15
  #sample3=15
  #empirical_dist<-empirical_dist[empirical_dist>600]select_ind<-sample(c(1:length(empirical_dist)),size = ngenes,replace = TRUE)
  select_mu<-empirical_dist[select_ind]
  select_disp<-sim_disp[select_ind]
  de<-rep(0,ngenes)
  de[1:as.integer(percent_DE*ngenes)]<-1
  #lib1<-rtruncnorm(15,a=1,mean=0,sd=1)
  #lib2<-rtruncnorm(15,a=1,mean=0,sd=1)
  #lib3<-rtruncnorm(15,a=1,mean=0,sd=1)
  lib1<-runif(15,min=0.90,1.10)
  lib2<-runif(15,min=0.90,1.10)
  lib3<-runif(15,min=0.90,1.10)
  #select_mu<-rep(1000,ngenes)
  lib<-c(lib1,lib2,lib3)
  lfc<-rtruncnorm(ngenes,a=0.5,mean=1,sd=1)
  comb<-matrix(c(rep(c(-1,0,1),as.integer(sum(de)/3)),rep(c(0,1,1),as.integer(sum(de)/3)),rep(c(1,-1,0),sum(de)-2*as.integer(sum(de)/3)),rep(c(0,0,0),ngenes-sum(de))),ncol=3,byrow = TRUE)
  power_comb<-(lfc*comb)
  log_mu_matrix<-log(select_mu,base=2)+power_comb
  log_matrix<-t(apply(log_mu_matrix,1,function(x) rep(x,each=15)))
  
  ind_module_matrix<-matrix(c(1:10,51:60,101:110,151:160,801:810),ncol=5)
  for(ind_mod in 1:5){
    log_mu<-log_mu_matrix[ind_module_matrix[,ind_mod],]
    
    #alpha<-a[a_ind]
    param_inwish<-(1-alpha)*diag(1,nrow=10,ncol=10)+alpha*matrix(1,nrow=10,ncol=10)
    sigma_prime<-riwish(60,param_inwish)
    sigma<-cov2cor(sigma_prime)
    
    log_matrix[ind_module_matrix[,ind_mod],1:15]<- t(rmvnorm(15,log_mu[1:10,1],sigma=sigma))
    
    param_inwish<-(1-alpha)*diag(1,nrow=10,ncol=10)+alpha*matrix(1,nrow=10,ncol=10)
    sigma_prime<-riwish(60,param_inwish)
    sigma<-cov2cor(sigma_prime)
    
    log_matrix[ind_module_matrix[,ind_mod],16:30]<- t(rmvnorm(15,log_mu[,2],sigma=sigma))
    
    param_inwish<-(1-alpha)*diag(1,nrow=10,ncol=10)+alpha*matrix(1,nrow=10,ncol=10)
    sigma_prime<-riwish(60,param_inwish)
    sigma<-cov2cor(sigma_prime)
    
    log_matrix[ind_module_matrix[,ind_mod],31:45]<- t(rmvnorm(15,log_mu[,3],sigma=sigma))
    
    #heatmap.2(mu_matrix/apply(mu_matrix,1,sum))
    
  }
  lib<-c(lib1,lib2,lib3)
  
  sim_matrix<-matrix(,nrow=ngenes,ncol=45)
  
  
  for(ind_gene in 1:ngenes){
    
    for(sample_ind in 1:45){
      
      sim_matrix[ind_gene,sample_ind]<-rnbinom(1,mu=lib[sample_ind]*2^(log_matrix[ind_gene,sample_ind]),size=select_disp[ind_gene])
      
      
    }
    
  }
  
  
  data<-sim_matrix
  return(data)
}
