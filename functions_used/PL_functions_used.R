## this document contains the functions used to conduct PL
##
U=function(gam,v_w)
{
  -gam * log(q_density(v_w))
}
chains_parallel_tempering=function(pot=U,  Sigma_mat = (1* diag(Dim)), init = rep(1,Dim), MSD_flag = NULL,effective_sigma_flag = NULL)
{
  accept_rate = rep(0,length(temps))
  num_cold_change <- 0
  MSD <- 0
  w= t( matrix(data = rep(init,length(temps)), nrow = length(init), ncol = length(temps) ) )
  w_array = array(0, dim = c(length(temps), length(init), iters))
  pro_matrix <- matrix(0, nrow = length(temps), ncol = Dim)
  #effective_sigma <- 55*2.38^2/(temps* Dim)
  #effective_sigma <- seq(from  = 5, to = 1 , length.out = num_level)
  #effective_sigma =temps*2.38^2/(Dim)
  #effective_sigma <- rep(55,num_level)
  effective_sigma = 2.38^2/(temps*Dim)
  effective_sigma[[1]] = effective_sigma[[1]] * special_factor
  for (i in 1:iters) {
    print('i=')
    print(i)
    # can=x+rnorm(length(temps),0,tune)
    if(is.null(effective_sigma_flag)){
      candidate = w + mvrnorm(n = length(temps),mu = rep(0,length(init)),Sigma = Sigma_mat)
    }else{
      for (jj in 1:length(temps)) {
        pro_matrix[jj, ] <- mvrnorm(n = 1 ,mu = rep(0,length(init)),Sigma = effective_sigma[jj] *Sigma_mat)
      }
      candidate = w + pro_matrix
    }
    logA = mapply(pot,temps,lapply(seq_len(nrow(w)), function(i) w[i,])) - mapply(pot,temps,lapply(seq_len(nrow(candidate)), function(i) candidate[i,]))
    # logA=unlist(Map(pot,temps,x))-unlist(Map(pot,temps,can))
    accept=(log(runif(length(temps)))<logA)
    accept_rate = accept_rate + accept
    w[accept , ]=candidate[accept , ]
    # now the coupling update
    swap=sample(1:length(temps),2)
    if(!is.null(MSD_flag)){
      if( (sum(swap==1)) > 0){
        num_cold_change = num_cold_change +1
      }
    }
    logA=pot(temps[swap[1]],w[swap[1], ])+pot(temps[swap[2]],w[swap[2], ])-
      pot(temps[swap[1]],w[swap[2], ])-pot(temps[swap[2]],w[swap[1], ])
    if (log(runif(1))<logA){
      w[swap, ]=w[rev(swap), ]
      if(!is.null(MSD_flag)){
        if( (sum(swap==1)) > 0){
          MSD = MSD + (temps[swap[1]] - temps[swap[2]])^2
        }
      }
    }
    # end of the coupling update
    w_array[, , i]= w
  }
  ##colnames(xmat)=paste("gamma=",temps,sep="")
  if(!is.null(MSD_flag)){
    accept_rate = accept_rate/iters
    return(list(swapefficient = MSD/num_cold_change, accept_rate = accept_rate ))
  }else{
    return(w_array)
  }
}


## get the highest tempe
U_fix_gam=function(gam)
{
  function(x) U(gam,x)
}
##
# select the scale of the inverse temperature
select_scale = function(initi_rho,beta_min, Sigma_mat, pot = U,init = rep(1,Dim)){
  w= t( matrix(data = rep(init,2), nrow = length(init), ncol = 2 ) )
  beta_0 = 1
  beta= beta_0
  rho = initi_rho
  beta_store = c(beta)
  accept_rate_store = c()
  alpha_plus_store = c()
  alpha_plus_list = c()
  while(beta > beta_min){
    alpha_plus_each = c()
    accept_rate = rep(0,2)
    alpha_plus = 1
    n = 0
    while((n<iner_maximum) & ( abs(alpha_plus - 0.23) >  0.01 )) {
      n = n+1
      beta_pri = beta * (1+exp(rho))^(-1)
      ## one step MH algorithm
      pro_matrix <- matrix(0, nrow = 2, ncol = Dim)
      temps = c(beta,beta_pri)
      effective_sigma = 2.38^2/(temps*Dim)
      if(beta == beta_0){
        effective_sigma[[1]] = effective_sigma[[1]] * special_factor
      }
      for (jj in 1:2) {
        pro_matrix[jj, ] <- mvrnorm(n = 1 ,mu = rep(0,Dim),Sigma = effective_sigma[jj] *Sigma_mat)
      }
      candidate = w + pro_matrix
      logA = mapply(pot,temps,lapply(seq_len(nrow(w)), function(i) w[i,])) - mapply(pot,temps,lapply(seq_len(nrow(candidate)), function(i) candidate[i,]))
      # logA=unlist(Map(pot,temps,x))-unlist(Map(pot,temps,can))
      accept=(log(runif(length(temps)))<logA)
      accept_rate = accept_rate + accept
      w[accept , ]=candidate[accept , ]
      swap = c(1,2)
      B_n =pot(temps[swap[1]],w[swap[1], ])+pot(temps[swap[2]],w[swap[2], ])-
        pot(temps[swap[1]],w[swap[2], ])-pot(temps[swap[2]],w[swap[1], ])
      print('B_n = ')
      print(B_n)
      print('beta_pri')
      print(beta_pri)
      alpha_plus = min(1, exp(B_n))
      alpha_plus_each = c(alpha_plus_each,alpha_plus)
      rho_plus = rho + (1/n)*(alpha_plus - 0.23)
      rho = rho_plus
    }
    accept_rate = accept_rate/iner_maximum
    accept_rate_store = c(accept_rate_store,accept_rate)
    beta = beta_pri
    alpha_plus_store = c(alpha_plus_store,alpha_plus)
    alpha_plus_list = rbind(alpha_plus_list,alpha_plus_each)
    beta_store = c(beta_store,beta)
  }
  
  out_put = list()
  out_put$accept_rate_store = accept_rate_store
  out_put$alpha_plus_store = alpha_plus_store
  out_put$beta_store = beta_store
  out_put$alpha_plus_list= alpha_plus_list
  
  return(out_put)
}












###