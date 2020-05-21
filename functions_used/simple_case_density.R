# simple case density
# q_density = function(x_vec){
#   if (x_vec[1] <= 10 & x_vec[1] >= -10 & x_vec[2]>= -10 & x_vec[2]<= 10) {
#     tmpphix = (1/3)*exp((-1/2)*sum((x_vec+5)^2)) + (2/3)*exp((-1/2)*sum((x_vec-5)^2))
#     return(tmpphix)
#   }else if(x_vec[1] <= 11 & x_vec[1] >= 10 & x_vec[2]>= -10 & x_vec[2]<= 10){
#     tmpphix = 0.05
#     return(tmpphix)
#   }else{
#     return(0)
#   }
# }

# only
# q_density = function(x_vec){
#   tmpphix = (1/3)*exp((-1/2)*sum((x_vec+5)^2)) + (2/3)*exp((-1/2)*sum((x_vec-5)^2))
#   return(tmpphix)
# }

# only multidimensional 
# q_density = function(x_vec){
#   tmpphix = exp((-1/2)*sum((x_vec+5)^2))
#   return(tmpphix)
# }
# multimode
# q_density = function(x_vec){
#   tmpphix = (1/15)*exp((-1/2)*sum((x_vec+1)^2)) + (2/15)*exp((-1/2)*sum((x_vec-2)^2))
#   +(3/15)*exp((-1/2)*sum((x_vec+3)^2))+(4/15)*exp((-1/2)*sum((x_vec-4)^2))+(5/15)*exp((-1/2)*sum((x_vec+5)^2))
#   return(tmpphix)
# }
# new version that mode is separated
q_density = function(x_vec){
  tmpphix = (1/15)*exp((-1/2)*sum((x_vec+11)^2)) + (2/15)*exp((-1/2)*sum((x_vec-12)^2)) + (3/15)*exp((-1/2)*sum((x_vec+8)^2))+(4/15)*exp((-1/2)*sum((x_vec-7)^2))+(5/15)*exp((-1/2)*sum((x_vec+2)^2))
  return(tmpphix)
}
q_density_est = function(x_vec,log = TRUE){
  result1 = rep(0,dim(x_vec)[1])
  for (i in (1:dim(x_vec)[1])) {
    x_vec1 = x_vec[i, ]
    tmpphix = (1/15)*exp((-1/2)*sum((x_vec1+11)^2)) + (2/15)*exp((-1/2)*sum((x_vec1-12)^2)) + (3/15)*exp((-1/2)*sum((x_vec1+8)^2))+(4/15)*exp((-1/2)*sum((x_vec1-7)^2))+(5/15)*exp((-1/2)*sum((x_vec1+2)^2))
    result1[i] = tmpphix
    }
  return(log(result1))
}


### ignore
## compare with MCMC with phi_mix as proposal. -proposal_density(y|x)
# Stand_MCMC = function(Num_sample, proposal_density, proposal_sample_fun, w_star = NULL,Dim_w,dens){
#   w_cur <- w_star
#   sample_get <- matrix(0,nrow = Num_sample, ncol = Dim_w)
#   for (i in 1:Num_sample) {
#     w_next <- proposal_sample_fun(w_cur,dens)
#     t_num = proposal_density(w_cur,dens)
#     t_den = proposal_density(w_next,dens)
#     qden_num = q_density(w_next)
#     qden_den = q_density(w_cur)
#     tmp_accept = (t_num*qden_num)/(t_den*qden_den)
#     p_accept = min(tmp_accept,1)
#     tmp_compare = runif(1,min = 0,max = 1)
#     if (tmp_compare <= p_accept) {
#       w_cur = w_next
#     }
#     sample_get[i, ] <- w_cur
#   }
#   return(sample_get)
# }

# proposal density function -- phi_mix
proposal_density = function(w,dens){
  tmpw <-as.vector(w)
  tmpw <- matrix(tmpw,nrow = 1)
  print('w')
  print(w)
  d2 = predict(dens, tmpw, what = c("dens"), logarithm = FALSE)
  phi_mix = d2[1]
  return(phi_mix)
}
# proposal sample function -- phi_mix
proposal_sample_fun = function(w = NULL,dens){
  component_num = dens$G
  pi_vec = dens$parameters$pro
  # mean_vec = dens$parameters$mean[variable_num,]
  # sigma_vec = c()
  # for (i in 1:component_num) {
  #   tmpsigma = dens$parameters$variance$sigma[,,i][variable_num,variable_num]
  #   sigma_vec = c(sigma_vec,tmpsigma)
  # }
  tmpsample1 = sample(x =1:component_num ,size = 1,replace = TRUE,prob = pi_vec)
  mean_cur = dens$parameters$mean[,tmpsample1]
  Sigma_cur = dens$parameters$variance
  Sigma_cur = Sigma_cur$sigma[,,tmpsample1]
  tmpsample2 = mvrnorm(n = 1,mean_cur,Sigma_cur)
  return(tmpsample2)
}
