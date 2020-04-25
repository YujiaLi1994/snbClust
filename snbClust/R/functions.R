estimate_beta_vector<-function(y,lib,phi,init_beta=10,pr_weight,max_iter=500){
  n<-dim(y)[2]
  p<-dim(y)[1]

  init_beta<-rep(init_beta,p)

  s_log<-matrix(rep(log(lib),p),byrow = F,nrow=n)
  phi_mat<-matrix(rep(phi,n),nrow=n,byrow=T)

  iter<-1
  pr_weight<-t(as.matrix(pr_weight,nrow=n))
  repeat{
    bet_mat<-matrix(rep(init_beta,n),nrow=n,byrow=T)

    lnmu<-s_log+bet_mat
    mu<-exp(lnmu)
    w<-phi_mat/(1+phi_mat/mu)
    tau<-lnmu+(t(y)/mu)-1
    beta<-pr_weight%*%(w*(tau-s_log))/(pr_weight%*%w)
    if(max(abs(beta-init_beta))<10^(-7)){break}
    init_beta<-beta
    if(iter>max_iter){break}
    iter<-iter+1
    print(max(abs(beta-init_beta)))

  }
  return(beta)
}







pen_func<-function(y,x,lib,pr_weights,beta_vector,global_beta,phi=1,tune=0){

  mean_beta<-global_beta
  penal_beta_vector<-c()
  for(i in 1:length(beta_vector)){
    mu<-lib*exp(x*beta_vector[i])
    resid<-(y-mu)/mu
    w<-mu/(1+(mu/phi))
    q_num<-pr_weights[,i]*w*x*(resid+x*beta_vector[i])
    q_denom<-pr_weights[,i]*w*(x^2)
    q<-sum(q_num)/sum(q_denom)
    penalty<-tune/sum(q_denom)

    if(beta_vector[i]>mean_beta){
      penal_beta_vector[i]<-ifelse(q>mean_beta+penalty,q-penalty,mean_beta)
    }
    else{penal_beta_vector[i]<-ifelse(q<mean_beta-penalty,q+penalty,mean_beta)}
  }

  return(penal_beta_vector)
}
mult_density1<-function(x,mu,phi=1){

  sum<-sum(dnbinom(x,mu=mu,size=phi,log=TRUE))
  return(sum)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

em_nb_mult2<-function(dat,k=2,phi=3,c_center=FALSE,gamma=0,no_init=1,ht=F,log_penal=F,diff_cutoff=10^(-7)){
  nan_val<-0
  # Reading the data as matrix( columns=genes, rows=sample)
  y<-as.matrix(dat)
  #Initial values of the component probabilities
  p<-dim(y)[2]
  pi<-rep(1/k,k)
  final_lik<-c()
  unpen_final_lik<-c()
  param_set<-list()

  for(m in 1:no_init){
    if(c_center==FALSE){
      #Initial value of the mean centers
      lambda<-matrix(,nrow=k,ncol=p)
      mu_y<-apply(y,2,mean)
      for(n in 1:length(mu_y)){

        lambda[,n]<-rpois(k,lambda=mu_y[n])

      }
      lambda[which(lambda==0)]<-10
    }
    else{lambda=c_center}
    # Computing value of unobserved values
    z<-matrix(,nrow=dim(y)[1],ncol=k)
    for(i in 1:dim(y)[1]){
      s<-c()
      for(j in 1:k){
        log_sum<-0
        for(l in 1:p){
          log_sum<-log_sum+dnbinom(y[i,l],mu=lambda[j,l],size=phi,log=T)
        }
        s[j]<-log_sum
      }
      s<-pi*brob(s)
      z[i,]<-as.numeric(s/sum(s))
      z[i,][is.na(z[i,])]<-.5
    }
    #calculating the log_liklihood
    log_lik<-0
    for(i in 1:dim(y)[1]){
      for(j in 1:k){
        joint_pdf<-0
        for(l in 1:p){
          joint_pdf<-joint_pdf+dnbinom(y[i,l],mu=lambda[j,l],size=phi,log=T)
        }
        log_lik<-log_lik+z[i,j]*(log(pi[j])+joint_pdf)

      }
    }
    up_log_lik<-log_lik
    log_lik<-0

    #####Looping until convergence
    while(abs(up_log_lik-log_lik)>diff_cutoff){
      log_lik<-up_log_lik
      #E-step
      z<-matrix(,nrow=dim(y)[1],ncol=k)
      for(i in 1:dim(y)[1]){
        s<-c()
        for(j in 1:k){
          log_sum<-0
          for(l in 1:p){
            log_sum<-log_sum+dnbinom(y[i,l],mu=lambda[j,l],size=phi,log=T)
          }
          s[j]<-log_sum
        }
        s<-pi*brob(s)
        z[i,]<-as.numeric(s/sum(s))
        z[i,][is.na(z[i,])]<-.5
      }
      #M-step
      for(j in 1:k){
        pi[j]<-sum(z[,j])/dim(y)[1]
        for(l in 1:p){
          lambda[j,l]<-sum(z[,j]*y[,l])/sum(z[,j])
        }
      }

      if(anyNA(lambda)){
        lambda<-param_set[[m-1]]$lambda
        z<-param_set[[m-1]]$z
        pi<-param_set[[m-1]]$prob
        nan_val<-m
      }
      var_clust_mean<-apply(lambda,2,mean)

      for(l in 1:p){
        if(log_penal==F){
          for(j in 1:k){
            penalty<-gamma*(lambda[j,l]*(lambda[j,l]+phi)*(k-1))/(phi*k*sum(z[,j]))
            if(lambda[j,l]>var_clust_mean[l]){
              pen_val<-lambda[j,l]-penalty
              if(ht==F){
                lambda[j,l]<-ifelse(pen_val>var_clust_mean[l],pen_val,var_clust_mean[l])}
              else{ lambda[j,l]<-ifelse(pen_val>var_clust_mean[l],lambda[j,l],var_clust_mean[l])}

            }
            else if(lambda[j,l]<var_clust_mean[l]){
              pen_val<-lambda[j,l]+penalty
              if(ht==F){
                lambda[j,l]<-ifelse(pen_val<var_clust_mean[l],pen_val,var_clust_mean[l])}
              else{ lambda[j,l]<-ifelse(pen_val<var_clust_mean[l],lambda[j,l],var_clust_mean[l])}



            }
          }
        }

        else{

          for(j in 1:k){
            penalty<-gamma*(lambda[j,l]*(lambda[j,l]+phi)/(phi*sum(z[,j])))*(1/lambda[j,l]-1/(k*var_clust_mean[l]))

            pen_val<-lambda[j,l]-(penalty*sign(lambda[j,l]-var_clust_mean[l]))
            if(ht==F){
              lambda[j,l]<-ifelse(abs(penalty)<abs(lambda[j,l]-var_clust_mean[l]),pen_val,var_clust_mean[l])
            }
            else{ lambda[j,l]<-ifelse(abs(penalty)<abs(lambda[j,l]-var_clust_mean[l]),lambda[j,l],var_clust_mean[l])
            }
          }
        }
      }
      ###Up-log_lik
      up_log_lik<-0
      for(i in 1:dim(y)[1]){
        for(j in 1:k){
          joint_pdf<-0
          for(l in 1:p){
            joint_pdf<-joint_pdf+dnbinom(y[i,l],mu=lambda[j,l],size=phi,log=T)
            #  print(dnbinom(y[i,l],mu=lambda[j,l],size=phi,log=T))
            # print(lambda[j,l])

          }
          up_log_lik<-up_log_lik+z[i,j]*(log(pi[j])+joint_pdf)
        }
      }
      #    print('up_log_lik')
      #    print(up_log_lik)
      un_pen_log_lik<-up_log_lik
      up_log_lik<-up_log_lik-gamma*sum(log(lambda)-log(var_clust_mean))
      #print(lambda)
      #   print('sum_log_lambda')
      #    print(sum(log(lambda)))
      # print('sum_penalty')
      # sum(penalty)
      print(abs(up_log_lik-log_lik))
      if(is.na(up_log_lik-log_lik)==TRUE){

        output<-list('z_init'='not_converge')
        return(output)
      }
      print(m)
    }
    unpen_final_lik[m]<-un_pen_log_lik
    final_lik[m]<-up_log_lik
    # print(lambda)
    param_set[[m]]<-list('lambda'=lambda,'z'=z,'prob'=pi,'lik'=up_log_lik,'un_pen_lik'=unpen_final_lik)
  }
  max<-max(final_lik)
  mode<-Mode(final_lik)
  parameter<-param_set[[which(final_lik==max)[1]]]

  list('param_set'=param_set,'paramter'=parameter)
}

mult_density<-function(x,mu,sigma){

  sum<-sum(dnorm(x,mean=mu,sd=sqrt(sigma),log=TRUE))
  return(sum)
}

cross_sum<-function(x,y){
  matrix1<-matrix(0,length(x),length(y))
  for(i in 1:dim(matrix1)[1]){
    for(j in 1:dim(matrix1)[2]){
      matrix1[i,j]<-x[i]+y[j]
    }
  }
  return(matrix1)
}

update_tau<-function(lib,beta,data,mu){
  tau_update<-matrix(0,length(lib),length(beta))
  for(i in 1:dim(tau_update)[1]){
    for(j in 1:dim(tau_update)[2]){
      tau_update[i,j]<-lib[i]+beta[j]+(data[i]-mu[i,j])/mu[i,j]
    }
  }
  return(tau_update)
}

soft_thresholding<-function(z,t){
  if((abs(z)-t)<0){
    return(0)
  }else{
    return(sign(z)*(abs(z)-t))
  }
}
update_keci<-function(tune=tune,penalty.option=penalty.option,keci,ita=ita,lambda=lambda){
  output<-rep(NA,length(keci))
  if(penalty.option=="lasso"){
    for(i in 1:length(output)){
      output[i]<-soft_thresholding(keci[i],lambda/tune)
    }
    return(output)
  }
  if(penalty.option=="MCP"){
    for(i in 1:length(output)){
      if(abs(keci[i])>lambda*ita){
        output[i]<-keci[i]
      }else{
        output[i]<-soft_thresholding(keci[i],lambda/tune)/(1-1/(ita*tune))
      }

    }
    return(output)
  }
}

mult_density2<-function(x,mu,phi=1){

  prod<-prod(dnbinom(x,mu=mu,size=phi,log=FALSE))
  return(prod)
}
calculate_penalty<-function(x,lambda,ita){
  y<-x
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      if(x[i,j]<lambda*ita){
        y[i,j]<-lambda*abs(x[i,j])-abs(x[i,j])*abs(x[i,j])/(2*ita)
      }else{
        y[i,j]<-0.5*ita*lambda^2
      }

    }
  }
  return(sum(y))
}
#phi<-rep(1,dim(data)[1])
update<-function(data,phi=1,lib,c_center=NULL,k=3,lambda=3,penalty.option="MCP",tune=1,ita=3,max_iter_IRLS=200,max_iter_ADMM=100,store_lang=NULL,store_delta=NULL){
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

    mu=c_center

  }

  z<-matrix(0,nrow=dim(data)[2],ncol=k)

  n<-dim(data)[2]

  mult_pdf<-matrix(0,nrow=n,ncol=k)

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

  ###estimatig betas
  ngenes<-dim(data)[1]
  beta_init<-log(c_center)
  #beta<-matrix(beta_init,nrow=ngenes,ncol=k)
  beta<-beta_init
  ###############IRLS

  iter_IRLS<-0
  max_iter_IRLS<-max_iter_IRLS
  store_lang=NULL
  store_delta=NULL
  if(is.null(store_lang)){
    store_lang<-matrix(runif(choose(k,2)*nrow(data),0,1),nrow=nrow(data),ncol=choose(k,2))
  }
  #store_lang<-matrix(runif(choose(k,2)*nrow(data),0,1),nrow=nrow(data),ncol=choose(k,2))
  if(is.null(store_delta)){
    store_delta<-matrix(runif(choose(k,2)*nrow(data),0,1),nrow=nrow(data),ncol=choose(k,2))
  }
  #store_delta<-matrix(runif(choose(k,2)*nrow(data),0,1),nrow=nrow(data),ncol=choose(k,2))

  for(j in 1:ngenes){
    #print(j)
    beta_pre<-beta[j,]
    data1<-data[j,]

    repeat{
      mu_update<-cross_sum(log(lib),beta_pre)
      mu_update<-exp(mu_update)
      v_update<-mu_update/(1+mu_update/phi[j])
      tau_update<-update_tau(lib=log(lib),beta=beta_pre,data=data[j,],mu=mu_update)
      ############################ADMM
      beta_update<-beta[j,]
      iter_ADMM<-0
      max_iter_ADMM<-max_iter_ADMM
      # init_delta<-runif(choose(k,2),0,1)
      # init_lang<-runif(choose(k,2),0,1)
      # delta_ADMM<-init_delta
      # lang_ADMM<-init_lang
      delta_ADMM<-store_delta[j,]
      lang_ADMM<-store_lang[j,]
      beta_ADMM_pre<-beta[j,]


      repeat{


        q_ADMM<-matrix(0,length(lib),dim(beta)[2])
        for(i in 1:dim(q_ADMM)[1]){
          q_ADMM[i,]<-tau_update[i,]-log(lib)[i]
        }
        q_ADMM<-as.vector(t(q_ADMM))

        W_ADMM<-v_update*z_init
        W_ADMM<-diag(as.vector(t(W_ADMM)))
        A_ADMM<-NULL
        for(i in 1:n){
          A_ADMM<-rbind(A_ADMM,diag(k))
        }
        #library(ppls)
        #Penalty.matrix(m=4,order=1)
        index_matrix<-combn(k,2)
        B_ADMM<-matrix(0,choose(k,2),k)
        for(i in 1:dim(B_ADMM)[1]){
          B_ADMM[i,][index_matrix[1,i]]<-1
          B_ADMM[i,][index_matrix[2,i]]<-(-1)
        }
        R_ADMM<-delta_ADMM-lang_ADMM/tune
        # beta_ADMM_new<-solve(t(A_ADMM)%*%W_ADMM%*%A_ADMM+tune*t(B_ADMM)%*%B_ADMM)%*%
        #  t((t(q_ADMM)%*%W_ADMM%*%A_ADMM+tune*t(R_ADMM)%*%B_ADMM))
        beta_ADMM_new<-solve(t(A_ADMM)%*%W_ADMM%*%A_ADMM+tune*t(B_ADMM)%*%B_ADMM)%*%
          (t(A_ADMM)%*%W_ADMM%*%q_ADMM+tune*t(B_ADMM)%*%R_ADMM)
        #beta_ADMM_pre<-beta_ADMM_new
        keci_ADMM<-B_ADMM%*%beta_ADMM_new+lang_ADMM/tune
        delta_ADMM<-update_keci(tune=tune,penalty.option=penalty.option,keci_ADMM,ita=ita,lambda=lambda)
        lang_ADMM<-lang_ADMM+tune*(B_ADMM%*%beta_ADMM_new-delta_ADMM)



        iter_ADMM<-iter_ADMM+1
        if(iter_ADMM>max_iter_ADMM){break}
        if(sum(abs(beta_ADMM_new-beta_ADMM_pre))<10^(-10)){break}
        #if(sum((B_ADMM%*%beta_ADMM_new-delta_ADMM)^2)<10^(-10)){break}
        #print(paste("ADMM",sum((B_ADMM%*%beta_ADMM_new-delta_ADMM)^2)))
        #if(sum(abs(beta_ADMM_new-beta_ADMM_pre))<10^(-10)){break}
        # beta_ADMM_pre<-beta_ADMM_new
        # keci_ADMM<-B_ADMM%*%beta_ADMM_new+lang_ADMM/tune
        # delta_ADMM<-update_keci(tune=tune,penalty.option=penalty.option,keci_ADMM,ita=ita,lambda=lambda)
        # lang_ADMM<-lang_ADMM+tune*(B_ADMM%*%beta_ADMM_new-delta_ADMM)
        #


      }
      #############################end of ADMM
      beta_new<-beta_ADMM_new
      iter_IRLS<-iter_IRLS+1
      #print(paste("IRLS",iter_IRLS))
      if(iter_IRLS>max_iter_IRLS){break}
      if(sum(abs(beta_new-beta_pre))<10^(-8)){break}
      beta_pre<-beta_new


    }
    beta[j,]<-beta_new
    store_delta[j,]<-delta_ADMM
    store_delta[j,]<-lang_ADMM

  }



  mult_pdf<-matrix(,nrow=n,ncol=k)
  for(i in 1:n){
    up_mu<-lib[i]*exp(beta)
    for(j in 1:k){
      mult_pdf[i,j]<-as.numeric(mult_density1(data[,i],mu=up_mu[,j],phi=phi))
    }
  }
  log_lik<-0
  max_mult<-apply(mult_pdf,1,max)
  for(i in 1:nrow(mult_pdf)){
    #log_lik<-log_lik+sum(pi*exp(mult_pdf[i,]))
    log_lik<-log_lik+sum(pi*exp(mult_pdf[i,]-max_mult[i]))+max_mult[i]
  }
  # mult_pdf<-mult_pdf%*%pi
  #
  # log_lik<-sum(log(mult_pdf))
  # mult_pdf<-colSums(mult_pdf)
  # mult_pdf<-mult_pdf+pi
  # log_lik<-sum(z%*%mult_pdf)


  # for(i in 1:n){
  #   for(j in 1:k){
  #     log_lik<-log_lik+z[i,j]*(log(pi[j])+mult_pdf[i,j])
  #   }
  # }


  log_lik_pen<-log_lik-calculate_penalty(B_ADMM%*%t(beta),lambda=lambda,ita=ita)
  return(list(z=z,beta=beta,log_lik=log_lik,log_lik_pen=log_lik_pen))
}


