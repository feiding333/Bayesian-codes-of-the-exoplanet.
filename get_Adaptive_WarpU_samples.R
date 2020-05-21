rm(list = ls())
library(FNN)
library(mvtnorm)
library(mnormt)
library(numDeriv)
library(stats4)
library(MASS)
library(mvtnorm)
library(numDeriv)
library(ggplot2)
library(mclust)
library(transport)
library(e1071)
library(sn)
library(RiskPortfolios)
#library(bridgesampling)
# source the function
## skew_t correnspond to density number 123
#density_fun_used = 'skew_t'
density_fun_used = "real_exo"
restrict_prior_flag = FALSE
density_num = 123
set.seed(23)
## source the GWL algorithm
if(density_fun_used == "skew_t"){
  source('./functions_used/GWLD_skewT.R')
}else if(density_fun_used == "real_exo"){
  load('./set posterior function/dataset1_one_planet_posterior.RData')
  source("./R model code/log_like.R")
  source("./R model code/log_post.R")
  source("./R model code/cov_fun.R")
  source("./R model code/cov_make.R")
  load('./prior_bounds')
  source("./R model code/planet_model.R") # Added in sqrt to phi
  if(restrict_prior_flag == TRUE){
    source("./R model code/restricted_priors.R")
  }else{
    source('./R model code/full_priors.R')
  }
  
}else{
  source('./functions_used/GWLD.R')
  }
  
source('./functions_used/Functions_AdaptiveWarpU.R')
if(density_fun_used == "skew_t"){
  source('./functions_used/skew_t_density.R')
}



N = 20000
sigma_I_11_factor = 1
#
if(density_fun_used == "skew_t"){
  alpha = 0.5
  sigma2 = 2.38^2/(10)
  stage_num =20
  sample_range = rbind(rep(-20,10), rep(20,10))
  sigma_I = diag(10)
  com_num = 25
}else if(density_fun_used == "real_exo"){
  com_num = 10
  stage_num = 4
  #sigma2 = 10
  sigma2 = 2.38^2/(7) * 5
  sigma2 = 1
  alpha1 = 0.5
  ###
  ###
  if(restrict_prior_flag == FALSE){
    Pmin <- rep(1.25,3)
    Pmax <- rep(10^4,3)
  }else{
    Pmax = c(44.6684, 12.8825, 10.7152)
    Pmin = c(39.8107, 11.4815, 10.0000)
  }
  num_planets <- 1
  dataset_num <- 1
  paras <- matrix(NA,2,7)
  paras[1,] <- c(C=-1,tau=42.4,K=0.1,e=0.1,w=2.0,M0=3.0,sigmaJ=2)
  paras[2,] <- c(C=-1,tau=42.4,K=0.1,e=0.1,w=2.0,M0=3.0,sigmaJ=3)
  p1(paras)
  ###
  ###
  sample_range = matrix(c(-2.2292,1.2552,41.67,42.49,1.832,3.381,0.007472,0.538570,0.006595,4.860418,0.000676,4.834212,1.057,1.788),2,7)
  load('stanData')
  True_Data = data
  sigma_I_vec = c()
  for(i in 1:(dim(True_Data)[2])){
    sigma_I_vec = c(sigma_I_vec,var(True_Data[,i]))
  }
  sigma_I = diag(sigma_I_vec)
}
sigma_I  = covEstimation(True_Data, control = list(type = 'ewma', lambda = 0.9))
sigma_I  = diag(sigma_I)
sigma_I = diag(sigma_I)
#sigma_I = diag(7)
sigma_I[1,1] = sigma_I[1,1] * sigma_I_11_factor

#
inserversion_num = 23 # 12: warpU+gibbs_changes  13: for gibbs only, 10000: only for warpU, 25: MH(Norm+Unif)+warpU,24: MH(Norm+Unif) only # 23 adaptive warp U method


time_begin = proc.time()
#
# mu2 = c(-5,-5,-5,-5,-5,-5)
# mu2 = c(-5,-5)


# Sigma2 = diag(10)

#
density_num = 123



N_Rpeat = 1
results_sample_list = c()
for( rep in 1:N_Rpeat){
  #results_sample = adaptive_warpU(com_num = com_num,N = N,inserversion_num = inserversion_num, stage_num=stage_num, sample_range = sample_range,density_num = density_num,alpha = alpha,sigma2 = sigma2,sigma_I = sigma_I,True_com = True_com,mean_dis_true = mean_dis_true,pro_dis_true = pro_dis_true,store_sample_dim = 2,Pre_stage_form = 1)
  results_sample = adaptive_warpU(com_num = com_num,N = N,inserversion_num = inserversion_num, stage_num=stage_num, sample_range = sample_range,density_num = density_num,alpha = alpha,sigma2 = sigma2,sigma_I = sigma_I,star_w = True_Data[2,])
  results_sample_list = c(results_sample_list,list(results_sample))
}

#last_num = c(20000,30000,40000,50000,60000)

#distance_last_change = distance_last_stage (com_num = com_num,N = N,inserversion_num = inserversion_num, stage_num=stage_num, sample_range = sample_range,True_Data = True_Data,density_num = density_num,last_num = last_num ,results_sample= results_sample,alpha = alpha,sigma2 = sigma2 )


  
time_used = proc.time() - time_begin
save(list = ls(),file = paste0('./Results_stored/May19_5_covest_set_start_AdaptiveWarpU_Skewt_N2_com_num_',com_num,'stage_num_',stage_num,'inser_',inserversion_num,'density_',density_num,'N=',N,'N_Reat',N_Rpeat,'density_fun_used_',density_fun_used,'sigma2_',sigma2,'restricted_flag_',restrict_prior_flag,'sigma_I_11_factor',sigma_I_11_factor,'diagnal_TRUE'))


