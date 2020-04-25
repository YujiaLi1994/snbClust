##'sparse Gaussian mixture model
##' @title sparse Gaussian mixture model .
##' @param data data matrix, where rows represent genes and columns represents samples.
##' @param K number of clusters.
##' @param c_center Initial centers for each clusters. We suggest using
##' other clustering algorithms such as sparse K-means to speed up.
##' @param lambda penalty parameter needs to be tuned or specified.
##' @return A list of two components:
##' \itemize{
##' \item{result: }{a list containing mean matrix (mu),variable estimation for each feature (sigma)
##' likelihood, cluster assignment probability (z) and prior probability for each cluster (pi)}
##' \item{BIC: }{ Extended BIC of the result}
##'
##' }
##' @references Pan W, Shen X. Penalized model-based clustering with application to variable selection[J]. Journal of Machine Learning Research, 2007, 8(May): 1145-1164.
##'
##'
##' @export

sgClust<-function(data,lambda,c_center=NULL,v_int=NULL,mu=0,Sigma=1,pi_int=NULL,K=2,no_init=1,max_iter=200){
  result_list<-list()
  log_lik_vector<-c()
  unpen_lik_vector<-c()
  time<-c()
  #Column refers to samples and row refers to genes
  p<-dim(data)[1] # number of variables
  n<-dim(data)[2] # number of subjects

  for(init in 1:no_init){
    start<-Sys.time()
    #========== E-step: initialization===========#

    if(length(c_center)==0){
      # assume the means for each cluster follows the distribution N(0,1)
      mu_int<-matrix(rnorm(p*K,mu,Sigma),nrow=p,ncol=K)
      #mu[which(mu==0)]<-100
    } else{
      mu_int=c_center
    }

    v_int<-ifelse(length(v_int)!=p,rep(Sigma,p),v_int) # initial value of variance of each variables
    pi_int<-ifelse(length(pi_int)!=K,rep(1/K,K),pi_int) # initial value of pi
    z_int<-matrix(,nrow=n,ncol=K) # initial value of prob in cluster k for each subject

    mult_pdf<-matrix(,nrow=n,ncol=K)
    for(i in 1:n){
      for(j in 1:K){
        mult_pdf[i,j]<-as.numeric(mult_density(data[,i],mu=mu_int[,j],sigma=v_int))
      }
    }

    for(i in 1:n){
      d<-brob(mult_pdf[i,])*pi_int
      z_int[i,]<-as.numeric(d/sum(d))
    }

    pi<-pi_int
    v<-v_int
    mu<-mu_int
    z<-z_int

    #========= M step ==========#
    iter=1
    log_lik<-1 # initialize
    log_lik_up<-0 #initialize
    while(abs(log_lik_up-log_lik)>10^(-7) & iter <=max_iter){
      # log likelihood before updating
      log_lik<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf[x,])))-lambda*sum(abs(mu))
      #update pi
      pi_up<-apply(z,2,sum)/n
      #pi_up<-ifelse(pi_up==0,1*10^(-16),pi_up)
      #update v
      update_v<-function(x){
        sig<-sum(sapply(1:K,function(y) sum(z[,y]*(data[x,]-mu[x,y])^2)))/n
        return(sig)
      }
      v_up<-sapply(1:p,function(x) update_v(x))

      #update mu
      mu_tu<-matrix(,nrow=p,ncol=K)
      mu_up<-matrix(,nrow=p,ncol=K)

      for(i in 1:K){
        temp<-sapply(1:n,function(x) z[x,i]*data[,x])
        mu_tu[,i]<-apply(temp,1,sum)/sum(z[,i])
        mu_up[,i]<-ifelse(lambda<=abs(apply(temp,1,sum)/v),sign(mu_tu[,i])*(abs(mu_tu[,i])-(lambda/sum(z[,i]))*v),0)
      }


      #update z
      z_up<-matrix(,nrow=n,ncol=K)

      for(i in 1:n){
        for(j in 1:K){
          mult_pdf[i,j]<-as.numeric(mult_density(data[,i],mu=mu_up[,j],sigma=v_up))
        }
      }

      for(i in 1:n){
        d<-brob(mult_pdf[i,])*pi_up
        z_up[i,]<-as.numeric(d/sum(d))
      }

      # update parameters values
      pi<-pi_up
      v<-v_up
      mu<-mu_up
      z<-z_up
      #log likelihood after updating
      log_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf[x,])))-lambda*sum(abs(mu))
      unpen_lik_up<-sum(sapply(1:n,function(x) z[x,]*(log(pi)+mult_pdf[x,])))
      iter=iter+1
    }

    end<-Sys.time()
    time[init]<-end-start
    result_list[[init]] <- list('mu'=mu,'Sigma'=v,'log_lik'=log_lik,'z'=z,'pi'=pi,'no_init'=init)
    log_lik_vector[init]<-log_lik_up
    unpen_lik_vector[init]<-unpen_lik_up
    print(init)
  }
  max_lik<-unpen_lik_vector[which.max(log_lik_vector)]

  # the optimal initials is the one with the maximum log likelihood
  optimal_result<-result_list[[which.max(log_lik_vector)]]
  max_mu<-optimal_result$mu
  s<-apply(max_mu,1,function(x) sum(x==0))
  #d<-sum(s[which(s!=1)])
  q<-sum(s)
  #BIC
  if(lambda!=0) {
    BIC<--2*max_lik+log(n+K+p+K*p-q)
  } else if(lambda==0) {
    BIC<--2*max_lik+log(n+K+p+K*p)
  }

  list('result'=optimal_result,'BIC'=BIC)
}
