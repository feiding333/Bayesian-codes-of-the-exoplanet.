# This is the newest version, corrected on July 14th, 2014
# EM algorithm in 1 dimension
#library(mvtnorm)

# In total seems to cost reps*n*(K+2K(iters-1)) Normal evaluations
# (ignores variance update which also involves matrix algebra)
EM_rep = function(x,K=1,reps=20,niters=650){
    if(K==1){n=dim(data)[1];return(list(mu=colSums(data)/n,sigma2=var(x)*(n-1)/n,weight=1))}
    log_lik_mix_combine = 1:(reps*2)
    EM_lists = list()
    for(rep in 1:reps){
        #set.seed(rep) Dave commented out 12th Mar 2018
        EM_lists[[rep]]=EM_hd(x,K,niters,tol=1e-8)
        log_lik_mix_combine[rep]=tail(EM_lists[[rep]]$log_lik_mix,1)
    }
    print("half EM is done!")
    for(rep in (1+reps):(2*reps)){
        #set.seed(rep) Dave commented out 12th Mar 2018
        EM_lists[[rep]]=EM_hd(x,K,niters,even_initial=F,tol=1e-8)
        log_lik_mix_combine[rep]=tail(EM_lists[[rep]]$log_lik_mix,1)
    }
    pick = which(log_lik_mix_combine==max(log_lik_mix_combine,na.rm=T))[1]
    mu = EM_lists[[pick]]$mu
    sigma2 = EM_lists[[pick]]$sigma2
    weight = EM_lists[[pick]]$weight
    
    return(list(EM_lists = EM_lists, mu=mu, sigma2=sigma2, weight=weight, log_lik_mix_all = log_lik_mix_combine))
}

# Seems to cost Kn+2Kn(iters-1) Normal evaluations
EM_hd = function(x,K=1,niters=650,even_initial=T,tol=1e-6){
    n = dim(x)[1]; Dim = dim(x)[2]
    if(is.vector(x)){# if x is 1-D, then, we just use EM_1d funciton
        return(EM_1d(x,K,niters,even_initial))
    }
    if(K==1){# if K==1, the result is easy 
        return(list(mu=colSums(data)/dim(data)[1],sigma2=var(x)*(n-1)/n,weight=1))
    }
    
    # if K is greater than 1 and Dim is greater than 1, then we do the following
    bounds = vapply(1:Dim,FUN=function(i){as.numeric(quantile(x[,i],c(0.25,0.75)))},
    FUN.VALUE=c(0,0))
    S_X = bounds[2,]-bounds[1,]
    large_spread_dim = which(S_X==max(S_X))
    
    # record the log_lik_mix to make sure the EM algorithm is correct and the likelihood is increasing! 
    log_lik_mix = rep(NA,niters)
    
    # set the initial values
    if(even_initial){
        bounds = as.numeric(quantile(x[,large_spread_dim],c(0,c(1:K)*0.95/K)+0.025))
        mu_init_index = vapply(1:K,FUN=function(i){
            sample(which(x[,large_spread_dim]<=bounds[i+1]&x[,large_spread_dim]>bounds[i]),1)
        },FUN.VALUE=1)
    }
    if(!even_initial){
        mu_init_index = sort(sample(1:n,K))
    }
    mu_curr = x[mu_init_index,]  # K by Dim
    prop_curr = rep(1,K)/K
    sigma2_curr = matrix(S_X^2*1.5,ncol=Dim,nrow=K,byrow = T)
    
    # Kn Normal evaluations 
    mixture_log = vapply(1:K,FUN=function(k){
       dmvnorm(x,mu_curr[k,],diag(sigma2_curr[k,]),log=T)+log(prop_curr[k])},FUN.VALUE=rep(1,n))
    max_mixture_log = apply(mixture_log,1,max)
    log_lik_mix[1] = mean(log(rowSums(exp(mixture_log-max_mixture_log))))+mean(max_mixture_log)

    # calculate 
    # 2Kn proposal evaluations per loop
    for(iter in 2:niters){
        # Kn proposal evaluations
        log_P = vapply(1:K,FUN=function(k){log(prop_curr[k])+dmvnorm(x,mu_curr[k,],diag(sigma2_curr[k,]),log=T)},
        FUN.VALUE=rep(1,n))
        max_log_P = apply(log_P,1,max)
        P = exp(log_P-max_log_P)
        P = P/rowSums(P)
        # EM updates
        PK = colSums(P)
        prop_curr = PK/sum(PK)
        mu_curr = t(P) %*% x / PK
        # What about this?
        sigma2_curr = t(vapply(1:K,FUN=function(k){((t(x)-mu_curr[k,])^2)%*%matrix(P[,k])+(S_X/K/2)^2*0.1},FUN.VALUE=rep(0,Dim)))/(PK+0.1)
        
        # calculate the likelihood
        # Kn proposal evalutions
        mixture_log = vapply(1:K,FUN=function(k){
           dmvnorm(x,mu_curr[k,],diag(sigma2_curr[k,]),log=T)+log(prop_curr[k])},FUN.VALUE=rep(1,n))
        max_mixture_log = apply(mixture_log,1,max)
        log_lik_mix[iter] = mean(log(rowSums(exp(mixture_log-max_mixture_log))))+mean(max_mixture_log)
        # stop the iterations if no improvement
        # Dave introduced tol 12th MAr 2018 -- want to reduce to 1e-6 to avoid large numebr of iters
        if(abs(log_lik_mix[iter]-log_lik_mix[iter-1])<abs(log_lik_mix[iter])*tol*Dim){break}
    }
    log_lik_mix = log_lik_mix[1:iter]
    
    order_prop = order(prop_curr)
    prop_curr = sort(prop_curr)
    mu_curr = mu_curr[order_prop,]
    sigma2_curr = sigma2_curr[order_prop,]
    
    return(list(mu=mu_curr,sigma2=sigma2_curr,weight=prop_curr,log_lik_mix=log_lik_mix/n,iter=iter))
}


# EM algorithm in 1 dimension
EM_1d = function(x,K=1,niters=650,even_initial=T){
    # if K==1, the result is easy
    if(K==1){
        valid_data = x[which(x<as.numeric(quantile(x,0.995)) & x>as.numeric(quantile(x,0.005)))]
        return(list(mu=mean(valid_data),sigma2=var(valid_data)))
    }
    
    # x is the vector of i.i.d samples from the target distribution
    n = length(x)
    bounds = as.numeric(quantile(x,c(0.25,0.75)))
    S_X = bounds[2] - bounds[1]
    
    log_lik_mix = rep(NA,niters)
    
    # initialize the parameters
    if(even_initial){
        bounds = as.numeric(quantile(x,c(0,c(1:K)*0.95/K)+0.025))
        mu_init_index = vapply(1:K,FUN=function(i){
            sample(which(x<=bounds[i+1]&x>bounds[i]),1)
        },FUN.VALUE=1)
        mu = x[mu_init_index]
    }
    if(!even_initial){
        mu = sample(x,K)
    }    
    
    prop = rep(1,K)/K
    sigma2 = rep(S_X^2*1.5,K)
    
    mixture_log = vapply(1:K,FUN=function(k){
       dnorm(x,mu[k],sd=sqrt(sigma2[k]),log=T)+log(prop[k])},FUN.VALUE=rep(1,n))
    max_mixture_log = apply(mixture_log,1,max)
    log_lik_mix[1] = mean(log(rowSums(exp(mixture_log-max_mixture_log))))+mean(max_mixture_log)
        
    # calculate 
    for(iter in 2:niters){
        log_P = vapply(1:K,FUN=function(k){log(prop[k])+dnorm(x,mu[k],sd=sqrt(sigma2[k]),log=T)},
        FUN.VALUE=rep(1,n))
        max_log_P = apply(log_P,1,max)
        P = exp(log_P-max_log_P)
        P = P/rowSums(P)
        # EM updates
        PK = colSums(P)
        prop = PK/sum(PK)
        mu = as.vector(t(P) %*% matrix(x))/PK
        sigma2 = (colSums(vapply(1:K,FUN=function(k){(x-mu[k])^2*P[,k]},FUN.VALUE=rep(1,n)))+(S_X/K/2)^2*0.1)/(PK+0.1)
        # calculate the likelihood
        mixture_log = vapply(1:K,FUN=function(k){
           dnorm(x,mu[k],sd=sqrt(sigma2[k]),log=T)+log(prop[k])},FUN.VALUE=rep(1,n))
        max_mixture_log = apply(mixture_log,1,max)
        log_lik_mix[iter] = mean(log(rowSums(exp(mixture_log-max_mixture_log))))+mean(max_mixture_log)
        # stop the iterations if no improvement
        if(abs(log_lik_mix[iter]-log_lik_mix[iter-1])<abs(log_lik_mix[iter])*1e-8){break}
    }
    log_lik_mix = log_lik_mix[1:iter]
    
    order_prop = order(prop)
    prop = sort(prop)
    mu = mu[order_prop]
    sigma2 = sigma2[order_prop]
    
    
    return(list(mu=mu,sigma2=sigma2,weight=prop,log_lik_mix=log_lik_mix))
}




