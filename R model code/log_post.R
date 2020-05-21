log_post <- function(for_post){
  
  # Inputs
  t <- for_post$t
  y <- for_post$y
  sd <- for_post$sd
  sigmaJ <- for_post$sigmaJ # Note changed to original scale to avoid alternativing normalizing constant
  C <- for_post$C
  cov_paras <- for_post$cov_paras
  planet_paras <- for_post$planet_paras
  
  mean_vec <- C
  if (num_planets>0){
    mean_vec <- numeric(length(y))
    #add
    mean_vec <- rep(C,length(y))
    ## add above
    for (planet_k in 1:num_planets){
      mean_vec <- mean_vec + v(t,planet_paras[[planet_k]])
    }
  }
  
  # Covaraince matrix
  k_full <- cov_make(t,t,cov_fun,cov_paras)
  obs_var <- diag(sd^2) 
  COV.obs.obs <- k_full + diag(rep(sigmaJ^2,length(y))) + obs_var
  
  if (prod(is.finite(COV.obs.obs)) == 1){ # Deal with numerical problems
    
    # Invert covaraince matrix
    chol_mat <- chol(COV.obs.obs) 
    COV.obs.obs.inverse  <- chol2inv(chol_mat)
    
    # Compute log-likelihood
    log_normalizing <- -0.5*determinant(COV.obs.obs, logarithm=TRUE)$modulus - length(t)*log(2*pi)/2
    quadratic <- -0.5*t(y-mean_vec)%*%COV.obs.obs.inverse%*%(y-mean_vec)
    value <- as.numeric(log_normalizing + quadratic)
    
    log_prior_value <- log_sigmaJprior(sigmaJ) + log_Cprior(C)
    if (num_planets>0){
      for (planet_k in 1:num_planets){
        log_prior_value <- log_prior_value + log_planet_prior(planet_paras[[planet_k]],planet_k)
      }
    }
    
    value <- value + log_prior_value
    
    if (is.finite(value)==FALSE){
      
      value <- -743.7469 # smallest value allowed by R (converted to log scale)
      
    }
    
  } else {
    
    value <- -743.7469 # smallest value allowed by R (converted to log scale)
    
  }
  
  return(value)
}
