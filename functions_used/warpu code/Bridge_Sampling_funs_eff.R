# the newest version of the bridge sampling functions, 
# last updated on Dec 7th, 2013


Opt_BS_hd = function(x1,x2,q1,q2,r0=1,opt_niters=50){
    #estimating log(c1/c2)
    l1_star = exp(q2(x1,log=T)-q1(x1,log=T))
    l2 = exp(q1(x2,log=T)-q2(x2,log=T))
    s1 = 1/(1+dim(x2)[1]/dim(x1)[1])
    s2 = 1-s1
    for(i in 1:opt_niters){
        r0 = mean(l2/(s1*l2+s2*r0))/mean(l1_star/(s1+s2*r0*l1_star))
    }
    print('r0')
    print(r0)
    return(log(r0))
}

Geo_BS_hd = function(x1,x2,q1,q2){
    #estimating log(c1/c2)
    l1 = exp((q2(x1,log=T)-q1(x1,log=T))/2)
    l2 = exp((q1(x2,log=T)-q2(x2,log=T))/2)
    log_r0 = log(mean(l2)) - log(mean(l1))
    return(log_r0)
}

phi_mix = function(x,log=F){
    log_phi_mix = vapply(1:K_prop,FUN=function(k){log(prop_prop[k])+dmvnorm(x,mu_prop[k,],diag(sigma2_prop[k,]),log=T)
        },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_phi_mix = apply(log_phi_mix,1,max)
    if(!log){
        return(rowSums(exp(log_phi_mix-max_log_phi_mix))*exp(max_log_phi_mix))
    }
    if(log){
        return(log(rowSums(exp(log_phi_mix-max_log_phi_mix)))+max_log_phi_mix)
    }
}

rphi_mix = function(n){
    theta = sample(c(1:K_prop),n,replace=T,prob=prop_prop)
    result = matrix(NA,ncol=Dim,nrow=n)
    for(k in 1:K_prop){
        index = which(theta==k)
        if(length(index)>0){result[index,]=rmvnorm(length(index),mu_prop[k,],diag(sigma2_prop[k,]))}
    }
    return(result)
}

prop_fun = function(x,k,log=F){
    yyy = x%*%diag(sqrt(sigma2_prop[k,]))+
    matrix(mu_prop[k,],byrow=T,nrow=dim(x)[1],ncol=Dim)
    
    log_prop_denominator = vapply(1:K_prop,FUN=function(kk){
        log(prop_prop[kk])+dmvnorm(yyy,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T) - zzz
    },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_prop_denominator = apply(log_prop_denominator,1,max)
    
    if(!log){
        return(prop_prop[k]/exp(max_log_prop_denominator)/rowSums(exp(log_prop_denominator-max_log_prop_denominator)))
    }
    if(log){
        return(log(prop_prop[k])-max_log_prop_denominator-log(rowSums(exp(log_prop_denominator-max_log_prop_denominator))))
    }
}

tp1 = function(x,log=T){
    # n standard Normal evaluations
    zzz = dmvnorm(x,rep(0,Dim),diag(rep(1,Dim)),log=T)
    log_tp1 = array(NA,c(dim(x)[1],K_prop))
    # nK^2 Normal evaluations, nK target evaluations
    for(k in 1:K_prop) {
        yyy = x%*%diag(sqrt(sigma2_prop[k,]))+
        matrix(mu_prop[k,],byrow=T,nrow=dim(x)[1],ncol=Dim)
        # nK evaluations
        log_prop_denominator = vapply(1:K_prop,FUN=function(kk){
            log(prop_prop[kk])+dmvnorm(yyy,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T) - zzz
        },FUN.VALUE=rep(0,dim(x)[1]))
        max_log_prop_denominator = apply(log_prop_denominator,1,max)
        # n target evaluations
        log_tp1[,k] = log(prop_prop[k])-max_log_prop_denominator-log(rowSums(exp(log_prop_denominator-max_log_prop_denominator))) +
        p1(yyy,log=T)
    }
    max_log_tp1 = apply(log_tp1,1,max)
    print('max_log_tp1')
    print(max_log_tp1)
    print('rowSums(exp(log_tp1-max_log_tp1)))')
    print(rowSums(exp(log_tp1-max_log_tp1)))
    print('log(rowSums(exp(log_tp1-max_log_tp1)))+max_log_tp1')
    print(log(rowSums(exp(log_tp1-max_log_tp1)))+max_log_tp1)
    return(log(rowSums(exp(log_tp1-max_log_tp1)))+max_log_tp1)
}

transform_prop = function(x,k){
    zzz = dmvnorm(x,mu_prop[k,],diag(sigma2_prop[k,]),log=T)
    log_denominator = vapply(1:K_prop,FUN=function(kk){
        log(prop_prop[kk])+(dmvnorm(x,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T)-zzz)
    },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_denominator = apply(log_denominator,1,max)
    return(prop_prop[k]/rowSums(exp(log_denominator-max_log_denominator))/exp(max_log_denominator))
}

rtp1 = function(x){
    T_prop = array(NA,c(dim(x)[1],K_prop))
    # K(K+1)n evaluations
    for(k in 1:K_prop) {
        # n evalutions
        zzz = dmvnorm(x,mu_prop[k,],diag(sigma2_prop[k,]),log=T)
        # nK evaluations
        log_denominator = vapply(1:K_prop,FUN=function(kk){
            log(prop_prop[kk])+(dmvnorm(x,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T)-zzz)
        },FUN.VALUE=rep(0,dim(x)[1]))
        max_log_denominator = apply(log_denominator,1,max)
        T_prop[,k] = prop_prop[k]/rowSums(exp(log_denominator-max_log_denominator))/exp(max_log_denominator)
    }    
    # dim1 by K_prop matrix
    theta = vapply(1:dim(x)[1],FUN=function(i){sample(1:K_prop,1,prob=T_prop[i,])},FUN.VALUE=1)
    t_x = x
    for(k in 1:K_prop){
        index = which(theta==k)
        if(length(index)>0){t_x[index,]=(x[index,]-matrix(mu_prop[k,],byrow=T,nrow=length(index),ncol=Dim))%*%diag(1/sqrt(sigma2_prop[k,]))}
    }
    return(t_x)
}

p2_warp0 = function(x,log=F){
    dmvnorm(x,mean=rep(0,Dim),sigma=diag(rep(1,Dim)),log=log)}

p2_warp1 = function(x,log=F){dmvnorm(x,mean=mu_warp1,sigma=diag(rep(1,Dim)),log=log)}

p2_warp2 = function(x,log=F){dmvnorm(x,mu_warp2,sigma2_warp2,log=log)}

p1_warp3 = function(x,log=F){
    matrix_result = vapply(c(-1,1),FUN=function(signx){p1(signx*x%*%sigma_warp2+matrix(mu_warp2,byrow=T,nrow=dim(x)[1],ncol=Dim),log=T)-log(2)},
        FUN.VALUE=rep(0,dim(x)[1]))
    max_matrix_result = apply(matrix_result,1,max)
    
    if(!log){
        return(rowSums(exp(matrix_result-max_matrix_result))*exp(max_matrix_result)*det(sigma_warp2))
    }
    if(log){
        return(log(rowSums(exp(matrix_result-max_matrix_result)))+max_matrix_result+determinant(sigma_warp2,log=T)$modulus[1])
    }
}
    

# commit 24/09/2019
# tp1_variant
tp1_variant = function(x,log=T){
  # n standard Normal evaluations
  zzz = dmvnorm(x,rep(0,Dim),diag(rep(1,Dim)),log=T)
  log_tp1 = array(NA,c(dim(x)[1],K_prop))
  
  ## sample one component by probability pi_k
  com_k = sample(x = (1:K_prop), size = 1,prob = prop_prop)
  
  ## one component 
  yyy = x%*%diag(sqrt(sigma2_prop[com_k,]))+
    matrix(mu_prop[com_k,],byrow=T,nrow=dim(x)[1],ncol=Dim)
  # nK evaluations
  log_prop_denominator = vapply(1:K_prop,FUN=function(kk){
    log(prop_prop[kk])+dmvnorm(yyy,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T) - zzz
  },FUN.VALUE=rep(0,dim(x)[1]))
  max_log_prop_denominator = apply(log_prop_denominator,1,max)
  # n target evaluations
  log_tp1_k = -max_log_prop_denominator-log(rowSums(exp(log_prop_denominator-max_log_prop_denominator))) + p1(yyy, log= T)
  return(log_tp1_k)
}


# repeated sample com_approx times
tp1_variant_1 = function(x,log=T){
  # n standard Normal evaluations
  log_tp1_mean = 0
  zzz = dmvnorm(x,rep(0,Dim),diag(rep(1,Dim)),log=T)
  log_tp1 = array(NA,c(dim(x)[1],K_prop))
  for (ii in 1:com_approx) {
    ## sample one component by probability pi_k
    com_k = sample(x = (1:K_prop), size = 1,prob = prop_prop)
    
    ## one component 
    yyy = x%*%diag(sqrt(sigma2_prop[com_k,]))+
      matrix(mu_prop[com_k,],byrow=T,nrow=dim(x)[1],ncol=Dim)
    # nK evaluations
    log_prop_denominator = vapply(1:K_prop,FUN=function(kk){
      log(prop_prop[kk])+dmvnorm(yyy,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T) - zzz
    },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_prop_denominator = apply(log_prop_denominator,1,max)
    # n target evaluations
    log_tp1_k = -max_log_prop_denominator-log(rowSums(exp(log_prop_denominator-max_log_prop_denominator))) + p1(yyy, log= T)

    log_tp1_mean = log_tp1_mean + exp(log_tp1_k)
  }
  log_tp1_mean  = log(log_tp1_mean  / com_approx)
  return(log_tp1_mean)
}


# sample com_approx one time
tp1_variant_2 = function(x,log=T){
  # n standard Normal evaluations
  zzz = dmvnorm(x,rep(0,Dim),diag(rep(1,Dim)),log=T)
  log_tp1 = array(NA,c(dim(x)[1],K_prop))
  log_tp1_mean = 0
  ## sample one component by probability pi_k
  com_k_all = sample(x = (1:K_prop), size = com_approx,prob = prop_prop,replace = FALSE)
  ii = 0
  for (com_k in com_k_all) {
    ii = ii+1
    ## one component 
    yyy = x%*%diag(sqrt(sigma2_prop[com_k,]))+
      matrix(mu_prop[com_k,],byrow=T,nrow=dim(x)[1],ncol=Dim)
    # nK evaluations
    log_prop_denominator = vapply(1:K_prop,FUN=function(kk){
      log(prop_prop[kk])+dmvnorm(yyy,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T) - zzz
    },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_prop_denominator = apply(log_prop_denominator,1,max)
    # n target evaluations
    log_tp1_k = -max_log_prop_denominator-log(rowSums(exp(log_prop_denominator-max_log_prop_denominator))) + p1(yyy, log= T)
    log_tp1_mean = log_tp1_mean + exp(log_tp1_k)
  }
  log_tp1_mean  = log(log_tp1_mean  / com_approx)
  return(log_tp1_mean)
}

# MSE_tp1_variant_1
MSE_tp1_variant_1 = function(x,log=T){
  # n standard Normal evaluations
  log_tp1_mean = 0
  zzz = dmvnorm(x,rep(0,Dim),diag(rep(1,Dim)),log=T)
  log_tp1 = array(NA,c(dim(x)[1],K_prop))
  for (ii in 1:com_approx) {
    ## sample one component by probability pi_k
    com_k = sample(x = (1:K_prop), size = 1,prob = prop_prop)
    
    ## one component 
    yyy = x%*%diag(sqrt(sigma2_prop[com_k,]))+
      matrix(mu_prop[com_k,],byrow=T,nrow=dim(x)[1],ncol=Dim)
    # nK evaluations
    log_prop_denominator = vapply(1:K_prop,FUN=function(kk){
      log(prop_prop[kk])+dmvnorm(yyy,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T) - zzz
    },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_prop_denominator = apply(log_prop_denominator,1,max)
    # n target evaluations
    log_tp1_k = -max_log_prop_denominator-log(rowSums(exp(log_prop_denominator-max_log_prop_denominator))) + p1(yyy, log= T)
    
    log_tp1_mean = log_tp1_mean + exp(log_tp1_k)
  }
  log_tp1_mean  = log(log_tp1_mean  / com_approx)
  MSE_tp1 = mean( (exp(log_tp1_mean) - exp( tp1(x, log = T) ))^2 )
  return(MSE_tp1)
}

MSE_tp1_variant_2 = function(x,log=T){
  # n standard Normal evaluations
  zzz = dmvnorm(x,rep(0,Dim),diag(rep(1,Dim)),log=T)
  log_tp1 = array(NA,c(dim(x)[1],K_prop))
  log_tp1_mean = 0
  ## sample one component by probability pi_k
  com_k_all = sample(x = (1:K_prop), size = com_approx,prob = prop_prop,replace = FALSE)
  ii = 0
  for (com_k in com_k_all) {
    ii = ii+1
    ## one component 
    yyy = x%*%diag(sqrt(sigma2_prop[com_k,]))+
      matrix(mu_prop[com_k,],byrow=T,nrow=dim(x)[1],ncol=Dim)
    # nK evaluations
    log_prop_denominator = vapply(1:K_prop,FUN=function(kk){
      log(prop_prop[kk])+dmvnorm(yyy,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T) - zzz
    },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_prop_denominator = apply(log_prop_denominator,1,max)
    # n target evaluations
    log_tp1_k = -max_log_prop_denominator-log(rowSums(exp(log_prop_denominator-max_log_prop_denominator))) + p1(yyy, log= T)
    log_tp1_mean = log_tp1_mean + exp(log_tp1_k)
  }
  log_tp1_mean  = log(log_tp1_mean  / com_approx)
  MSE_tp1 = mean( (exp(log_tp1_mean) - exp(tp1(x,log=T)))^2 )
  return(MSE_tp1)
}


# MSE_original
MSE_tp1 = function(x,log=T){
  zzz = dmvnorm(x,rep(0,Dim),diag(rep(1,Dim)),log=T)
  log_tp1 = array(NA,c(dim(x)[1],K_prop))
  # nK^2 Normal evaluations, nK target evaluations
  for(k in 1:K_prop) {
    yyy = x%*%diag(sqrt(sigma2_prop[k,]))+
      matrix(mu_prop[k,],byrow=T,nrow=dim(x)[1],ncol=Dim)
    # nK evaluations
    log_prop_denominator = vapply(1:K_prop,FUN=function(kk){
      log(prop_prop[kk])+dmvnorm(yyy,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T) - zzz
    },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_prop_denominator = apply(log_prop_denominator,1,max)
    # n target evaluations
    log_tp1[,k] = log(prop_prop[k])-max_log_prop_denominator-log(rowSums(exp(log_prop_denominator-max_log_prop_denominator))) +
      p1(yyy,log=T)
  }
  max_log_tp1 = apply(log_tp1,1,max)
  log_tp1_mean = log(rowSums(exp(log_tp1-max_log_tp1)))+max_log_tp1

  MSE_tp1 = mean( (exp(log_tp1_mean) - exp(tp1(x,log = T)))^2 )
  return(MSE_tp1)
}



## stochastic opt Bridge stimation
# Stoc_Opt_BS_hd = function(x1,x2,q1,q2,r0=1,opt_niters=12){
#   #estimating log(c1/c2)
#   if(is.list(x1)){
#     K = length(x1)
#   }
#   r0_n = 0
#   r0_d = 0
#   num_n = 0
#   num_d = 0
#   for (kk in 1:K) {
#     num_n = num_n + dim(x1[[kk]])[1]
#     num_d = num_d + dim(x2)[1]
#     l1_star = exp(q2(x1[[kk]],log=T)-q1(x1[[kk]],log=T))
#     l2 = exp(q1[[kk]](x2,log=T)-q2(x2,log=T))
#     s1 = 1/(1+dim(x2)[1]/dim(x1)[1])
#     s2 = 1-s1
#     r0_n = sum(l2/(s1*l2+s2*r0)) + r0_n
#     r0_d = sum(l1_star/(s1+s2*r0*l1_star)) + r0_d
#   }
#   r0 = (r0_n/num_n) / (r0_d/num_d)
# 
#   print('r0')
#   print(r0)
#   return(log(r0))
# }



qtilde_k = function(x,com_k, log=T){
  # n standard Normal evaluations
  zzz = dmvnorm(x,rep(0,Dim),diag(rep(1,Dim)),log=T)
  log_tp1 = array(NA,c(dim(x)[1],K_prop))
  ## one component 
  yyy = x%*%diag(sqrt(sigma2_prop[com_k,]))+
    matrix(mu_prop[com_k,],byrow=T,nrow=dim(x)[1],ncol=Dim)
  # nK evaluations
  log_prop_denominator = vapply(1:K_prop,FUN=function(kk){
    log(prop_prop[kk])+dmvnorm(yyy,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T) - zzz
  },FUN.VALUE=rep(0,dim(x)[1]))
  max_log_prop_denominator = apply(log_prop_denominator,1,max)
  # n target evaluations
  log_tp1_k = -max_log_prop_denominator-log(rowSums(exp(log_prop_denominator-max_log_prop_denominator))) + p1(yyy, log= T)
  return(log_tp1_k)
  #return(log(rowSums(exp(log_tp1_k-max_log_tp1)))+max_log_tp1)
}



split_data_k = function(x_in){
  split_data = list()
  length(split_data) = K_prop
  n_x = dim(x_in)[1]
  Ind_com = rep(0, n_x)
  prob_matrix = matrix(0, nrow = n_x, ncol = K_prop)
  for (k in 1:K_prop) {
    prob_matrix[ , k ] = prop_prop[k] * exp(qtilde_k(x = x_in, com_k = k, log = T) - tp1(x = x_in, log = T))
  }
  for (i in 1:n_x) {
    num_use = sample(1:K_prop, size = 1, prob = prob_matrix[ i,])
    Ind_com[i] = num_use
  }
  for (k in 1:K_prop) {
    split_data[[k]] = x_in[ Ind_com == k,]
  }
  return(split_data)
}
# this is the stochastics Bridge estimation function with same normal data for different r_k
Opt_BS_hd_Sto = function(split_data){
  rk = rep(0, K_prop)
  for (k in 1:K_prop) {
    qk = function(x, log= T){
      qtilde_k(x,k, log=T)
    }
    rk[k] = (Opt_BS_hd(x1=split_data[[k]], x2=data_normal, q1 = qk ,q2=p2_warp0,r0=10))
  }
  r = sum(exp(rk) * prop_prop)
  return(r)
}

## for different number of normal distribution
## this is the stochastics Bridge estimation function
Opt_BS_hd_Sto_splitNorm_Error = function(split_data,split_norm){
  rk = rep(0, K_prop)
  for (k in 1:K_prop) {
    qk = function(x, log= T){
      qtilde_k(x,k, log=T)
    }
    if(is.null(split_norm[[k]])){
      rk[k] = 0
    }else{
      if( dim(split_data[[k]])[1] == 1 ){
        rk[k] = 0
      }else{
        rk[k] = (Opt_BS_hd(x1=split_data[[k]], x2=split_norm[[k]], q1 = qk ,q2=p2_warp0,r0=10))
        
      }
      
    }
  }
  r = sum(exp(rk) * prop_prop)
  return(r)
}
Opt_BS_hd_Sto_splitNorm = function(split_data,split_norm){
  rk = rep(0, K_prop)
  for (k in 1:K_prop) {
    qk = function(x, log= T){
      qtilde_k(x,k, log=T)
    }
    rk[k] = (Opt_BS_hd(x1=split_data[[k]], x2=split_norm[[k]], q1 = qk ,q2=p2_warp0,r0=10))
  }
  r = sum(exp(rk) * prop_prop)
  return(r)
}

## different pi_k
Opt_BS_hd_Sto_splitNorm_pik = function(split_data,split_norm){
  rk = rep(0, K_prop)
  pik = rep(0,K_prop)
  for (k in 1:K_prop) {
    qk = function(x, log= T){
      qtilde_k(x,k, log=T)
    }
    rk[k] = (Opt_BS_hd(x1=split_data[[k]], x2=split_norm[[k]], q1 = qk ,q2=p2_warp0,r0=10))
    pik[k] = dim(split_norm[[k]])[1]
  }
  pik = pik/sum(pik)
  r = sum(exp(rk) * pik)
  return(r)
}
  
## split start with w to tilde{w}_k
rtp1_split = function(x){
  split_tilde_w  = list()
  length(split_tilde_w) = K_prop
  T_prop = array(NA,c(dim(x)[1],K_prop))
  # K(K+1)n evaluations
  for(k in 1:K_prop) {
    # n evalutions
    zzz = dmvnorm(x,mu_prop[k,],diag(sigma2_prop[k,]),log=T)
    # nK evaluations
    log_denominator = vapply(1:K_prop,FUN=function(kk){
      log(prop_prop[kk])+(dmvnorm(x,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T)-zzz)
    },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_denominator = apply(log_denominator,1,max)
    T_prop[,k] = prop_prop[k]/rowSums(exp(log_denominator-max_log_denominator))/exp(max_log_denominator)
  }    
  # dim1 by K_prop matrix
  theta = vapply(1:dim(x)[1],FUN=function(i){sample(1:K_prop,1,prob=T_prop[i,])},FUN.VALUE=1)

  for(k in 1:K_prop){
    index = which(theta==k)
    if(length(index)>0){split_tilde_w[[k]]=(x[index,]-matrix(mu_prop[k,],byrow=T,nrow=length(index),ncol=Dim))%*%diag(1/sqrt(sigma2_prop[k,]))}
    
  }
  return(split_tilde_w)
}

rtp1_split = function(x){
  split_tilde_w  = list()
  length(split_tilde_w) = K_prop
  T_prop = array(NA,c(dim(x)[1],K_prop))
  # K(K+1)n evaluations
  for(k in 1:K_prop) {
    # n evalutions
    zzz = dmvnorm(x,mu_prop[k,],diag(sigma2_prop[k,]),log=T)
    # nK evaluations
    log_denominator = vapply(1:K_prop,FUN=function(kk){
      log(prop_prop[kk])+(dmvnorm(x,mu_prop[kk,],diag(sigma2_prop[kk,]),log=T)-zzz)
    },FUN.VALUE=rep(0,dim(x)[1]))
    max_log_denominator = apply(log_denominator,1,max)
    T_prop[,k] = prop_prop[k]/rowSums(exp(log_denominator-max_log_denominator))/exp(max_log_denominator)
  }    
  # dim1 by K_prop matrix
  theta = vapply(1:dim(x)[1],FUN=function(i){sample(1:K_prop,1,prob=T_prop[i,])},FUN.VALUE=1)
  
  for(k in 1:K_prop){
    index = which(theta==k)
    if(length(index)>0){split_tilde_w[[k]]=(x[index,]-matrix(mu_prop[k,],byrow=T,nrow=length(index),ncol=Dim))%*%diag(1/sqrt(sigma2_prop[k,]))}
    
  }
  return(split_tilde_w)
}
rtp1_split_equal = function(x){
  split_tilde_w  = list()
  length(split_tilde_w) = K_prop
  data_new = rtp1(x)
  theta = vapply(1:dim(x)[1],FUN=function(i){sample(1:K_prop,1)},FUN.VALUE=1)
  
  for(k in 1:K_prop){
    index = which(theta==k)
    if(length(index)>0){split_tilde_w[[k]]=data_new[index,]}
    
  }
  return(split_tilde_w)
}

# compare function - input boostrap_num, N2_vec- number of normal distribution, output - MSE, Bias^2, Var(hat_r), mean_r
compare_method <- function(boostrap_num, N2_vec, w_in){
  N2_len <- length(N2_vec)
  r1_mean <- rep(0, N2_len)
  r2_mean <- rep(0, N2_len)
  r1_Bias2 <- rep(0, N2_len)
  r2_Bias2 <- rep(0, N2_len)
  r1_Var <- rep(0, N2_len)
  r2_Var <- rep(0, N2_len)
  r1_MSE <- rep(0, N2_len)
  r2_MSE <- rep(0, N2_len)
  r1_boos <- rep(0, boostrap_num)
  r2_boos <- rep(0, boostrap_num)
  for ( i_in in 1:N2_len) {
    N2_in = N2_vec[i_in]
    # Normal sample for original samples
    data_normal_ori = rmvnorm(N2_in,mean=rep(0,Dim),sigma=diag(1,Dim))
    # Normal sample for our samples
    split_norm = list()
    length(split_norm) = K_prop
    for (kk in 1:K_prop) {
      split_norm[[kk]] = rmvnorm(N2_in,mean=rep(0,Dim),sigma=diag(1,Dim))
      
    }
    for (j_in in 1:boostrap_num) {
      data_normal = data_normal_ori
      r1_boos[j_in] = method_original(w_in = w_in, normal_in = data_normal_ori)
      r2_boos[j_in] = method_w_Diff_normal_in(w_in = w_in, normal_in = split_norm)
    }
    print('r1_boos')
    print(r1_boos)
    print('r2_boos')
    print(r2_boos)
    r1_mean[i_in] <- mean(r1_boos)
    r2_mean[i_in] <- mean(r2_boos)
    r1_Var[i_in] <- var(r1_boos)*(boostrap_num-1)/(boostrap_num)
    r2_Var[i_in] <- var(r2_boos)*(boostrap_num-1)/(boostrap_num)
    r1_Bias2[i_in] <- (mean(r1_boos) - 2*pi)^2
    r2_Bias2[i_in] <- (mean(r2_boos) - 2*pi)^2
    r1_MSE[i_in] <- mean((r1_boos -2*pi)^2)
    r2_MSE[i_in] <- mean((r2_boos - 2*pi)^2)
  }
  return(list(r1_mean = r1_mean, r2_mean= r2_mean, r1_Var= r1_Var, r2_Var = r2_Var, r1_Bias2 = r1_Bias2, r2_Bias2 = r2_Bias2, r1_MSE = r1_MSE, r2_MSE = r2_MSE))
}


# method original 
method_original = function(w_in, normal_in){
  data_trans = rtp1(w_in) 
  return(exp(Opt_BS_hd(x1=data_trans,x2= normal_in, q1=tp1,q2=p2_warp0,r0=10)))
}
# using different normal distribution
method_w_Diff_normal_in = function(w_in, normal_in){
  tildeW_k  = rtp1_split(w_in)
  return(Opt_BS_hd_Sto_splitNorm(split_data = tildeW_k,split_norm = normal_in))
}

method_w_Diff_normal_in_pik = function(w_in, normal_in){
  tildeW_k  = rtp1_split(w_in)
  return(Opt_BS_hd_Sto_splitNorm_pik(split_data = tildeW_k,split_norm = normal_in))
}

