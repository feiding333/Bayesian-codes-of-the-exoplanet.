load('./Results_stored/Data_Jones.Rdata')
library(sn)
# function H(x)/phi(x), the given energy function
p1_single <- function(x,log=F){
  if(is.matrix(x)){
    if(dim(x)[2] != 10){
      x = t(x)
      print('fine')
    }
  }else{
    if(length(x) == 10){
    }else{
      print('############')
    }
  }
  log_results = sapply(1:lweights,function(i){log(weights[i])+dmst(x,xis[[i]],Omegas[[i]],alphas[[i]],nus[i],log=T)})
  max_log_results = max(log_results)
  if(!log){
    return(sum(exp(log_results-max_log_results))*exp(max_log_results))
  }
  if(log){
    return(log(sum(exp(log_results-max_log_results)))+max_log_results)
  }
}
phix = function(x_vec){
  if ((sum(x_vec >= -20) == length(x_vec)) & (sum(x_vec <= 20) == length(x_vec)) ) {
    tmpphix = p1_single(x_vec)
    return(tmpphix)
  }else if(x_vec[1] >= 20 & x_vec[1] <= 21 & x_vec[2] >= -20 & x_vec[2] <= 20 & x_vec[3] >= -20 & x_vec[3] <= 20 & x_vec[4] >= -20 & x_vec[4] <= 20 & x_vec[5] >= -20 & x_vec[5] <= 20 & x_vec[6] >= -20 & x_vec[6] <= 20 & x_vec[7] >= -20 & x_vec[7] <= 20 & x_vec[8] >= -20 & x_vec[8] <= 20 & x_vec[9] >= -20 & x_vec[9] <= 20 & x_vec[10] >= -20 & x_vec[10] <= 20){
    tmpphix = 1/(40^9)
    return(tmpphix)
  }else{
    return(0)
  }
}
# function indx(x*, H), map from x* to the number of the energy level
get_indx = function(x_samp,phix){
  #interval = seq(-4.9,30,by = 0.1)
  interval = seq(from = 10.83, length.out = 350 ,by = 0.1)
  tmplgphix = -log(phix(x_samp))
  if (x_samp[1] >= 20 & x_samp[1] <= 21 & x_samp[2] >= -20 & x_samp[2] <= 20 & x_samp[3] >= -20 & x_samp[3] <= 20 & x_samp[4] >= -20 & x_samp[4] <= 20 & x_samp[5] >= -20 & x_samp[5] <= 20 & x_samp[6] >= -20 & x_samp[6] <= 20 & x_samp[7] >= -20 & x_samp[7] <= 20 & x_samp[8] >= -20 & x_samp[8] <= 20 & x_samp[9] >= -20 & x_samp[9] <= 20 & x_samp[10] >= -20 & x_samp[10] <= 20) {
    return(0)
  }else if(tmplgphix <= 10.83){
    return(1)
  }else if(tmplgphix > 45.73){
    return(351)
  }else{
    for (i in 1:349) {
      interval_star = interval[i]
      interval_end = interval[i]+0.1
      if (tmplgphix <= interval_end & tmplgphix > interval_star) {
        return(i+1)
      }
    }
  }
}
update_delta = function(delta_s){
  delta_s = sqrt(delta_s+1) - 1
  return(delta_s)
}
get_estimate = function(g_init,delta_1,delta_end,get_indx,phix,update_delta,n_1,MH_steps = NULL,stage_Max = 10,proposal_specific = 0){
  s = 0
  m = length(tar_pi)
  delta_s = delta_1
  n_s = n_1
  g = g_init
  theta = log(g)
  x_cur = runif(10,-10,10)
  tmp_ind_cur = get_indx(x_cur,phix)+1
  distance_z = c()
  x_final_stage = c()
  x_final_index = c()
  if(is.null(MH_steps)){steps_MH =1}else{steps_MH = MH_steps} 
  #while (delta_s > delta_end) {
  for (stage_repeat in 1:stage_Max) {
    delta_s_next = update_delta(delta_s)
    vist_frequence = rep(0,m)
    for ( k  in 1:(steps_MH*n_s)) {
      #gamma_t = t_0/(max(t_0,t))
      e_cur = rep(0,m)
      if(proposal_specific == 0){
        x_new = mvrnorm(n = 1,x_cur, 0.5*diag(rep(Dim)))
        flag1 = (x_new[1] >= -20 & x_new[1] <= 21 & x_new[2] >= -20 & x_new[2] <= 20 & x_new[3] >= -20 & x_new[3] <= 20 & x_new[4] >= -20 & x_new[4] <= 20 & x_new[5] >= -20 & x_new[5] <= 20 & x_new[6] >= -20 & x_new[6] <= 20 & x_new[7] >= -20 & x_new[7] <= 20 & x_new[8] >= -20 & x_new[8] <= 20 & x_new[9] >= -20 & x_new[9] <= 20 & x_new[10] >= -20 & x_new[10] <= 20)
      }else{
        flag1 = TRUE
      }
      if (flag1) {
        if(proposal_specific == 0){
          t_num = dmvnorm(x_cur,mean= x_new,sigma = 0.5*diag(rep(Dim)))
          t_den = dmvnorm(x_new,mean= x_cur,sigma = 0.5*diag(rep(Dim)))
          phi_num = phix(x_new)
          phi_den = phix(x_cur)
          tmp_ind_new = get_indx(x_new,phix)+1
          theta_num = theta[tmp_ind_cur]
          theta_den = theta[tmp_ind_new]
          #g_num = g[tmp_ind_cur]
          #g_den = g[tmp_ind_new]
          tmp_accept = (t_num*phi_num*exp(theta_num))/(t_den*phi_den*exp(theta_den))
          p_accept = min(tmp_accept,1)
          
        }else{
          norm_unif_index = sample(1:2,size = 1,prob = c(40/41,1/41))
          if(norm_unif_index == 1){
            x_new = rp1(1)
            if((sum(x_new >= -20) == length(x_new)) & (sum(x_new <= 20) == length(x_new)) ){
              x_new = x_new
            }else{
              x_new = rp1(1)
            }
          }else{
            x_new = c()
            for (unif_i in 1:length(x_cur)) {
              if(unif_i == 1){
                tmp_new_star = runif(n = 1,min = 20,max = 21)
              }else{
                tmp_new_star = runif(n = 1,min = -20,max = 20)
              }
              x_new = c(x_new,tmp_new_star)
            }
          }
          #
          t_num = q_pro_special(w_pre =x_new ,w_star = t(x_cur))
          t_den = q_pro_special(w_pre =t(x_cur) ,w_star = x_new)
          
          phi_num = phix(x_new)
          phi_den = phix(x_cur)
          tmp_ind_new = get_indx(x_new,phix)+1
          theta_num = theta[tmp_ind_cur]
          theta_den = theta[tmp_ind_new]
          tmp_accept = (t_num*phi_num*exp(theta_num))/(t_den*phi_den*exp(theta_den))
          p_accept = min(tmp_accept,1)
          
        }
        

        tmp_compare = runif(1,min = 0,max = 1)
        if (tmp_compare <= p_accept) {
          x_cur = x_new
          tmp_ind_cur = tmp_ind_new
          if(k%%steps_MH == 0){
            e_cur[tmp_ind_cur] = 1
            # caution
            #g = g*exp(gamma_t*(e_cur-tar_pi))
            
            #theta = theta + gamma_t*(e_cur- tar_pi)
            #theta = theta + delta_s*(e_cur- tar_pi)
            # causion
            g[tmp_ind_cur] = g[tmp_ind_cur] + delta_s * g[tmp_ind_cur]
            #g[tmp_ind_cur] = g[tmp_ind_cur] * exp(delta_s * g[tmp_ind_cur])
            vist_frequence[tmp_ind_cur] = vist_frequence[tmp_ind_cur]+1
            theta = log(g)
          }
        }else{
          if(k%%steps_MH == 0){
            e_cur[tmp_ind_cur] = 1
            #g = g*exp(gamma_t*(e_cur-tar_pi))
            #theta = theta + gamma_t*(e_cur- tar_pi)
            #theta = theta + delta_s*(e_cur- tar_pi)
            # causion
            g[tmp_ind_cur] = g[tmp_ind_cur] + delta_s * g[tmp_ind_cur]
            #g[tmp_ind_cur] = g[tmp_ind_cur] * exp(delta_s * g[tmp_ind_cur])
            
            vist_frequence[tmp_ind_cur] = vist_frequence[tmp_ind_cur]+1
            
            theta = log(g)
          }

        }
      }else{
        if(k%%steps_MH == 0){
          e_cur[tmp_ind_cur] = 1
          #g = g*exp(gamma_t*(e_cur-tar_pi))
          #theta = theta + gamma_t*(e_cur- tar_pi)
          #theta = theta + delta_s*(e_cur- tar_pi)
          # causion
          g[tmp_ind_cur] = g[tmp_ind_cur] + delta_s * g[tmp_ind_cur]
          #g[tmp_ind_cur] = g[tmp_ind_cur] * exp(delta_s * g[tmp_ind_cur])
          vist_frequence[tmp_ind_cur] = vist_frequence[tmp_ind_cur]+1
          theta = log(g)
        }
      }
      # causion plus
      if (max(theta) > 1e2) {
        theta = theta - 1e2
        g = exp(theta)
      }
      if (min(theta) < -1e2) {
        #print("***")
        #print(theta)
        theta = theta + 1e2
        g = exp(theta)
      }
      #if(delta_s_next <= delta_end){
      if(stage_repeat == stage_Max){
        x_final_stage = rbind(x_final_stage,x_cur)
        x_final_index = c(x_final_index, tmp_ind_cur)
        
      }
      
      #print('flag12')
      # g = g/sum(g)
      # print(g)
      # theta = log(g)
      
    }
    
    # if (max(abs(theta)) > 1e1) {
    #   theta = theta - 1e1
    #   g = exp(theta)
    # }
    # if (min(abs(theta)) < 1e-1) {
    #   #print(theta)
    #   theta = theta + 1e-1
    #   g = exp(theta)
    # }
    # g = g/sum(g)
    # theta = log(g)
    s = s+1
    n_s = 1.1*n_s
    #n_s = 1.1 * n_s
    #delta_s = sqrt(1+delta_s)-1
    delta_s = update_delta(delta_s)


    estimate_z = (sum(g)-g[1])/g[1]
    distance_z = c(distance_z,estimate_z)
    
  }
  result = list()
  result$distance_z = distance_z
  result$vist_frequence = vist_frequence
  result$g = g
  result$x_final_stage = x_final_stage
  result$x_final_index = x_final_index
  return(result)
}


q_pro_special = function(w_pre,w_star){
  if ((sum(w_star >= -20) == length(w_star)) & (sum(w_star <= 20) == length(w_star)) ){
    s = (40/41)*p1_single(w_star)
  }else{
    s = (1/41)*(1/(40^9))
  }
  return(s)
}



