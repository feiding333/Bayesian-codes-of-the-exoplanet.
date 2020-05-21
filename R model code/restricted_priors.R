# P prior - truncated Jeffrey's prior
max_num_planets <- 3
Pmin <- numeric(max_num_planets)
Pmax <- numeric(max_num_planets)
for (k in 1:max_num_planets){
  Pmin[k] <- prior_bounds[k,3] #1.25
  Pmax[k] <- prior_bounds[k,4] #10^4
}
log_Pprior <- function(x,k){
  if (x>=Pmin[k] & x<=Pmax[k]){
    value <- -log(x) - log(log(Pmax[k]/Pmin[k]))
  } else {
    value <- -Inf
  }
  return(value)
}

# K prior - truncated modified Jeffrey's prior
K0 <- 1
Kmax <- 999
log_Kprior <- function(x){
  if (x>0 & x<=Kmax){
    value <- -log(K0*(1+x/K0)) - log(log(1+Kmax/K0))
  } else {
    value <- -Inf
  }
  return(value)
}

# M0 prior - uniform
log_M0prior <- function(x){
  if (x>=0 & x<(2*pi)){
    value <- -log(2*pi)
  } else {
    value <- -Inf
  }
  return(value)
}

# e prior - truncated Rayleigh distribution
emax <- 1
sigma_e <- 0.2
log_eprior <- function(x){
  if (x>=0 & x<emax){
    value <- log(x) - 2*log(sigma_e) - x^2/(2*sigma_e^2) - log(1-exp(-emax^2/(2*sigma_e^2)))
  } else {
    value <- -Inf
  }
  return(value)
}

# w prior (argument of pericenter) - uniform
log_wprior <- function(x){
  if (x>=0 & x<(2*pi)){
    value <- -log(2*pi)
  } else {
    value <- -Inf
  }
  return(value)
}

# Planet prior
log_planet_prior <- function(keps,k){
  value <- log_Pprior(keps$tau,k) + log_Kprior(keps$K) + log_eprior(keps$e) + log_wprior(keps$w) + log_M0prior(keps$M0)
  return(value)
}

# sigmaJ prior - truncated modified Jeffrey's prior
sigmaJ0 <- 1
sigmaJmax <- 99
log_sigmaJprior <- function(x){
  if (x>0 & x<=sigmaJmax){
    value <- -log(sigmaJ0*(1+x/sigmaJ0)) - log(log(1+sigmaJmax/sigmaJ0))
  } else {
    value <- -Inf
  }
  return(value)
}

# C prior (offset) - uniform
Cmax <- 1000
log_Cprior <- function(x){
  if (x>=-Cmax & x<=Cmax){
    value <- -log(2*Cmax)
  } else {
    value <- -Inf
  }
  return(value)
}

