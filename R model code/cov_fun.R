cov_fun <- function(t1,t2,params){
  # Quasi-periodic kernel
  tau <- exp(params[1])
  alpha <- exp(params[2])
  lambda_e <- exp(params[3])
  lambda_p <- exp(params[4])
  t.diff <- t1-t2
  phi <- t.diff*pi/tau
  value <- alpha^2*exp(-(sin(phi))^2/(2*lambda_p^2) - t.diff^2/(2*lambda_e^2))
  return(value)
}