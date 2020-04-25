##'sparse negative binomial mixture model for clustering RNA-seq count data with fused lasso penalty
##' @title sparse negative binomial mixture model for clustering RNA-seq count data with fused lasso penalty.
##' @param data data matrix, where rows represent genes and columns represents samples.
##' @param phi dispersion parameter for each gene.
##' @param lib library size normalization factor for each sample.
##' @param k number of clusters.
##' @param center Initial centers for each clusters. We suggest using
##' other clustering algorithms such as sparse K-means to speed up.
##' @param lambda penalty parameter needs to be tuned or specified.
##' @param penalty.option penalty type, the default is MCP.
##' @param tune a parameter controlling the ADMM algorithm
##' @param ita a parameter in MCP penalty
##' @return A list of two components:
##' \itemize{
##' \item{result: }{a list containing mean matrix (beta),
##' likelihood, penalized log likelihood and cluster assignment probability (z)}
##' \item{BIC: }{ Extended BIC of the result}
##'
##' }
##' @references Rahman T, Li Y, Ma T, et al. A sparse negative binomial mixture model for clustering RNA-seq count data[J]. arXiv preprint arXiv:1912.02399, 2019.
##' @export


em_glm_Yujia<-function(data,phi=1,lib,c_center=NULL,k=3,lambda,penalty.option="MCP",tune=1,ita=3,max_iter_IRLS=20,max_iter_ADMM=50,max_iter_EM=10,threshold=0.005,BIC.ita=1){
  start=Sys.time()
  result<-update(data,phi=phi,lib=lib,c_center=c_center,k=k,lambda=lambda,penalty.option=penalty.option,tune=tune,ita=ita,max_iter_IRLS=max_iter_IRLS,max_iter_ADMM=max_iter_ADMM)
  #beta_list[[1]]<-result$beta
  #z_list[[1]]<-result$z
  #ind<-2
  beta_pre<-result$beta
  log_lik_pre<-result$log_lik_pen
  no_iter_EM<-1
  result_list<-NULL
  repeat{
    #print(result$z)
    #print(result$beta)

    #print(no_iter_EM)
    result<-update(data,phi=phi,lib=lib,c_center=exp(beta_pre),k=k,lambda=lambda,penalty.option=penalty.option,tune=tune,ita=ita,max_iter_IRLS=max_iter_IRLS,max_iter_ADMM=max_iter_ADMM,store_lang = result$lang,store_delta = result$delta)
    #beta_list[[ind]]<-result$beta
    #z_list[[ind]]<-result$z
    #ind<-ind+1
    beta_new<-result$beta
    log_lik_new<-result$log_lik_pen
    #result_list[[no_iter_EM]]<-list(lik=log_lik_pre,log_lik=result$log_lik,beta=beta_pre)
    print(log_lik_new)


    # if(class(try(if(abs(updated_log_lik-log_lik)<10^(-7)) break))=='try-error'){
    #
    #   output<-list('z_init'=z_init,'beta_init'=beta_init,'z_matrix'=z_matrix,'z'=z,'beta_matrix'=beta_list,'beta'=beta,'pi'=pi)
    #   return(output)
    # }
    #if(abs(log_lik_new-log_lik_pre)<10^(-2)){break}
    if(abs(log_lik_new-log_lik_pre)<1e-2){break}
    if(no_iter_EM>max_iter_EM){break}
    no_iter_EM<-no_iter_EM+1
    log_lik_pre<-log_lik_new
    beta_pre<-beta_new
  }
  beta<-beta_new
  log_lik<-result$log_lik

  d<-0
  for(i in 1:nrow(data)){
    d<-d+length(unique(round(beta[i,],4)))
  }
  d1<-1
  for(j in 1:nrow(data)){
    d1<-d1*choose(5,length(unique(round(beta[j,],4))))
  }
  n<-ncol(data)
  #BIC<--2*log_lik+log(n)*(k+sum(d)-1)
  BIC<--2*log_lik+log(n)*(k+d-1)+2*BIC.ita*log(d1)
  end=Sys.time()
  time=difftime(end,start,units=c("secs"))
  list('result'=result,'BIC'=BIC)
}
