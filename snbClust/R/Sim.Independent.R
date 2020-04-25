##' @title Simulation 1
##' @param ngenes total number of genes
##' @param eff_a effect size
##' @param percent_DE percentage of DE genes
##' @param sim_disp gene dispersion from breast cancer data
##' @param empirical_dist average gene expression from breast cancer data
##' @return data with independent genes for simulation 1 and 2
##' @references Rahman T, Li Y, Ma T, et al. A sparse negative binomial mixture model for clustering RNA-seq count data[J]. arXiv preprint arXiv:1912.02399, 2019.
##' @export
##'
Sim.Independent<-function(ngenes=150,eff_a=1,percent_DE=1,sim_disp,empirical_dist){
  quantile<-quantile(empirical_dist,0.70)
  empirical_dist<-empirical_dist[empirical_dist<quantile]
  #sim_disp<-data$disp
  sim_disp<-sim_disp[empirical_dist<quantile]
  select_ind<-sample(c(1:length(empirical_dist)),size = ngenes,replace = TRUE)
  select_mu<-empirical_dist[select_ind]
  select_disp<-sim_disp[select_ind]
  de<-rep(0,ngenes)
  de[1:as.integer(percent_DE*ngenes)]<-1
  lib1<-runif(15,min=1,1)
  lib2<-runif(15,min=1,1)
  lib3<-runif(15,min=1,1)
  lib<-c(lib1,lib2,lib3)
  lfc<-rtruncnorm(ngenes,a=eff_a/2,mean=eff_a,sd=1)
  comb<-matrix(c(rep(c(-1,0,1),as.integer(sum(de)/3)),rep(c(0,1,1),as.integer(sum(de)/3)),rep(c(1,-1,0),sum(de)-2*as.integer(sum(de)/3)),rep(c(0,0,0),ngenes-sum(de))),ncol=3,byrow = TRUE)
  power_comb<-2^(lfc*comb)
  mu_matrix<-select_mu*power_comb
  sim_matrix1<-sim_matrix2<-sim_matrix3<-matrix(,ncol=15,nrow=ngenes)
  for(i in 1:ngenes){

    sim_matrix1[i,]<-rnbinom(15,mu=lib1*mu_matrix[i,1],size=select_disp[i])
    sim_matrix2[i,]<-rnbinom(15,mu=lib2*mu_matrix[i,2],size=select_disp[i])
    sim_matrix3[i,]<-rnbinom(15,mu=lib3*mu_matrix[i,3],size=select_disp[i])
  }
  sim_matrix<-cbind(sim_matrix1,sim_matrix2,sim_matrix3)
  data<-sim_matrix
  return(data)
}
