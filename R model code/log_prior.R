log_prior <- function(for_prior){
  
  sigmaJ <- exp(for_post$log_sigmaJ)
  C <- for_post$C
  planet_paras <- for_post$planet_paras
  
  log_prior_value <- log_sigmaJprior(sigmaJ) + log_Cprior(C)
  for (planet_k in 1:num_planets){
    log_prior_value <- log_prior_value + log_planet_prior(planet_paras[[planet_k]],planet_k)
  }
  
  return(log_prior_value)
  
}