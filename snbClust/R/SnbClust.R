##'sparse negative binomial mixture model for clustering RNA-seq count data
##' @title sparse negative binomial mixture model for clustering RNA-seq count data.
##' @param data data matrix, where rows represent genes and columns represents samples.
##' @param phi dispersion parameter for each gene.
##' @param lib library size normalization factor for each sample.
##' @param k number of clusters.
##' @param center Initial centers for each clusters. We suggest using
##' other clustering algorithms such as sparse K-means to speed up.
##' @param tune penalty parameter needs to be tuned or specified.
##' @return A list of two components:
##' \itemize{
##' \item{result: }{a list containing mean matrix (beta),
##' likelihood, cluster assignment probability (z) and prior probability for each cluster (pi)}
##' \item{BIC: }{ Extended BIC of the result}
##'
##' }
##' @references Rahman T, Li Y, Ma T, et al. A sparse negative binomial mixture model for clustering RNA-seq count data[J]. arXiv preprint arXiv:1912.02399, 2019.
##' @export

snbClust<-function(data,phi=rep(1,dim(data)[1]),lib,c_center=NULL,k=2,penalize=FALSE,tune=0,no_init=1,max_iter_beta=200,max_iter=100,lambda=100,m_node=F){
  if(m_node==T){
    library(snowfall)



  }

  result_list<-list()
  log_lik_vector<-c()
  un_pen_lik_vector<-c()
  time<-c()
  z_matrix<-beta_list<-list()
  for(init in 1:no_init){
    start<-Sys.time()
    no_iter<-1
    error<-0
    #Column refers to samples and row refers to genes
    pi<-c()
    pi[1:k]<-1/k
    #  lib<-apply(data,2,sum)
    if(length(c_center)==0){
      mu<-matrix(rpois(dim(data)[1]*k,lambda =lambda),nrow=dim(data)[1],ncol=k)
      mu[which(mu==0)]<-100
    } else{
      if(is.vector(c_center)==T){c_center<-as.matrix(c_center,ncol=1)}
      mu=c_center

    }
    z<-matrix(,nrow=dim(data)[2],ncol=k)

    n<-dim(data)[2]

    mult_pdf<-matrix(,nrow=n,ncol=k)

    for(i in 1:n){

      for(j in 1:k){
        mult_pdf[i,j]<-as.numeric(mult_density1(data[,i],mu=lib[i]*mu[,j],phi=phi[j]))
      }
    }
    for(i in 1:n){
      d<-brob(mult_pdf[i,])*pi
      z[i,]<-as.numeric(d/sum(d))
    }
    # print(z)
    z_init<-z
    ###Qstep
    pi<-apply(z,2,sum)/n
    #####estimating global beta
    genes<-dim(data)[1]
    global_beta<-as.vector(estimate_beta_vector(y=data,lib=lib,phi=phi,pr_weight=rep(1,dim(data)[2]),max_iter=max_iter_beta))


    ###estimatig betas

    beta<-matrix(,nrow=genes,ncol=k)

    for(j in 1:k){

      z_node<-z[,j]

      beta[,j]<-as.vector(estimate_beta_vector(y=data,lib=lib,phi=phi,pr_weight=z_node,max_iter=max_iter_beta))


    }


    if(penalize==T){
      penalized_beta<-matrix(,nrow=genes,ncol=k)
      for(i in 1:genes){
        penalized_beta[i,]<-pen_func(data[i,],x=rep(1,n),lib=lib,pr_weights=z,beta_vector=beta[i,],global_beta=global_beta[i],phi=phi[i],tune=tune)
      }
      beta<-penalized_beta

    }
    beta_init<-beta
    ##updated log-multivariate_pdf by sample and cluster(colums) sample(row)
    mult_pdf<-matrix(,nrow=n,ncol=k)
    for(i in 1:n){
      up_mu<-lib[i]*exp(beta)
      for(j in 1:k){
        mult_pdf[i,j]<-as.numeric(mult_density1(data[,i],mu=up_mu[,j],phi=phi))
      }
    }
    log_lik<-0
    for(i in 1:n){
      for(j in 1:k){
        log_lik<-log_lik+z[i,j]*(log(pi[j])+mult_pdf[i,j])
      }
    }
    log_lik<-log_lik-tune*sum(abs(beta))
    updated_log_lik<-log_lik
    repeat{
      ###Updating Z matrix
      for(i in 1:n){
        d<-brob(mult_pdf[i,])*pi
        z[i,]<-as.numeric(d/sum(d))
      }
      #Updating pi
      pi<-apply(z,2,sum)/n
      #updating beta
      for(j in 1:k){

        z_node<-z[,j]



        beta[,j]<-as.vector(estimate_beta_vector(y=data,lib=lib,phi=phi,pr_weight=z_node,max_iter=max_iter_beta))



      }
      if(penalize==T){
        penalized_beta<-matrix(,nrow=genes,ncol=k)
        for(i in 1:genes){
          penalized_beta[i,]<-pen_func(data[i,],x=rep(1,n),lib=lib,pr_weights=z,beta_vector=beta[i,],global_beta=global_beta[i],phi=phi[i],tune=tune)

        }
        beta<-penalized_beta
      }
      ##updated log-multivariate_pdf by sample and cluster(colums) sample(row)
      mult_pdf<-matrix(,nrow=n,ncol=k)
      for(i in 1:n){
        up_mu<-lib[i]*exp(beta)
        for(j in 1:k){
          mult_pdf[i,j]<-as.numeric(mult_density1(data[,i],mu=up_mu[,j],phi=phi))
        }
      }
      log_lik<-0
      for(i in 1:n){
        for(j in 1:k){
          log_lik<-log_lik+z[i,j]*(log(pi[j])+mult_pdf[i,j])
        }
      }
      un_penalize_lik<-log_lik
      log_lik<-log_lik-tune*sum(abs(beta))
      #    print(log_lik)
      no_iter<-no_iter+1
      if(class(try(if(abs(updated_log_lik-log_lik)<10^(-7)) break))=='try-error'){

        output<-list('z_init'=z_init,'beta_init'=beta_init,'z_matrix'=z_matrix,'z'=z,'beta_matrix'=beta_list,'beta'=beta,'pi'=pi)
        return(output)
      }
      if(abs(updated_log_lik-log_lik)<10^(-7)) break;
      if(no_iter>max_iter) break;
      updated_log_lik<-log_lik
      z_matrix[[no_iter-1]]<-z
      beta_list[[no_iter-1]]<-beta
      print(abs(updated_log_lik-log_lik))
      print(no_iter)
    }
    end<-Sys.time()
    time[init]<-end-start
    result_list[[init]]<-  list('beta'=beta,'log_lik'=log_lik,'z'=z,'pi'=pi)
    log_lik_vector[init]<-log_lik
    un_pen_lik_vector[init]<-un_penalize_lik
    print(init)
  }
  max_lik<-2*un_pen_lik_vector[which.max(log_lik_vector)]


  optimal_result<-result_list[[which.max(log_lik_vector)]]
  max_beta<-optimal_result$beta
  s<-apply(max_beta,1,function(x) length(unique(x)))
  #  d<-sum(s[which(s!=1)])
  d<-sum(s)
  BIC<--2*max_lik+log(n)*(k+d-1)

  list('result'=optimal_result,'BIC'=BIC)
}
