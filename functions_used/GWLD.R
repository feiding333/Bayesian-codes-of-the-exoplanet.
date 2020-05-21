# function H(x)/phi(x), the given energy function
phix = function(x_vec){
  if (x_vec[1] <= 10 & x_vec[1] >= -10 & x_vec[2]>= -10 & x_vec[2]<= 10) {
    tmpphix = (1/3)*exp((-1/2)*sum((x_vec+5)^2)) + (2/3)*exp((-1/2)*sum((x_vec-5)^2))
    return(tmpphix)
  }else if(x_vec[1] <= 11 & x_vec[1] >= 10 & x_vec[2]>= -10 & x_vec[2]<= 10){
    tmpphix = 0.05
    return(tmpphix)
  }else{
    return(0)
  }
}
# function indx(x*, H), map from x* to the number of the energy level
get_indx = function(x_samp,phix){
  interval = seq(-4.9,30,by = 0.1)
  tmplgphix = -log(phix(x_samp))
  if (x_samp[1] >= 10 & x_samp[1] <= 11 & x_samp[2] >= -10 & x_samp[2] <= 10) {
    return(0)
  }else if(tmplgphix <= -4.9){
    return(1)
  }else if(tmplgphix > 30){
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
get_estimate = function(g_init,delta_1,delta_end,get_indx,phix,update_delta,n_1,MH_steps = NULL,stage_Max = 10){
  s = 0
  m = length(tar_pi)
  delta_s = delta_1
  n_s = n_1
  g = g_init
  theta = log(g)
  x_cur = runif(2,-10,10)
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
      x_new = mvrnorm(n = 1,x_cur, 20*diag(c(1,1)))
      e_cur = rep(0,m)
      if (x_new[1] >= -10 & x_new[1] <= 11 & x_new[2] >= -10 & x_new[2] <= 10) {
        t_num = dmvnorm(x_cur,mean= x_new,sigma = 20*diag(c(1,1)))
        t_den = dmvnorm(x_new,mean= x_cur,sigma = 20*diag(c(1,1)))
        phi_num = phix(x_new)
        phi_den = phix(x_cur)
        tmp_ind_new = get_indx(x_new,phix)+1
        theta_num = theta[tmp_ind_cur]
        theta_den = theta[tmp_ind_new]
        #g_num = g[tmp_ind_cur]
        #g_den = g[tmp_ind_new]
        tmp_accept = (t_num*phi_num*exp(theta_num))/(t_den*phi_den*exp(theta_den))
        p_accept = min(tmp_accept,1)
        print("p_accept=")
        print(p_accept)
        print("p_accept=")
        print(p_accept)
        print("flagtheta")

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
      if(delta_s_next <= delta_end){
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






