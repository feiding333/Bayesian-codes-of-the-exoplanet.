rm(list = ls())
library(doRNG)
library(foreach)
library(doParallel)
##
# our method
#####################################################
library(mvtnorm)
library(mnormt)
library(numDeriv)
library(stats4)
library(MASS)
library(mvtnorm)
library(smfsb)
library(ggplot2)
library(LaplacesDemon)
library(ggsci)
library(microbenchmark)
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
library(RiskPortfolios)
#library(bridgesampling)
# source the function
# source('./functionused/Functionsexo_new.R')
# source('./functionused/simple_case_density.R')
## skew_t correnspond to density number 123
# density_fun_used = 'skew_t'
# density_num = 'skewt'
density_fun_used = "real_exo"
density_num = "real_exo"
density_num = 123
set.seed(23)
if(density_fun_used == "skew_t"){
  source('./functions_used/GWLD_skewT.R')
  source('./functions_used/skew_t_density.R')
}else if(density_fun_used == "real_exo"){
  load('./set posterior function/dataset1_one_planet_posterior.RData')
  source("./R model code/log_like.R")
  source("./R model code/log_post.R")
  source("./R model code/cov_fun.R")
  source("./R model code/cov_make.R")
  source("./R model code/planet_model.R") # Added in sqrt to phi
  source('./functions_used/Functions_AdaptiveWarpU.R')
}else{
  source('./functions_used/GWLD.R')
}
source('./functions_used/PL_functions_used.R')


#
if(density_fun_used == "skew_t"){
  Dim = 10
  num_level = 8
  iters=  N = 1e6
}else if(density_fun_used == "real_exo"){
  N = 5e4
  Dim = 7
  num_level = 7
  iters=5e4
}




# ### 
# ## calculate the sampling efficiency fo PT
# Sample_efficiency = c()
# for (num_level in 2:10) {
#   print('num_level')
#   print(num_level)
#   temps <- 10^( -2 * ((1:num_level) -1)/(num_level-1))
#   iters=2e4
#   Sample_efficiency_tmp = chains_parallel_tempering(MSD_flag = 1,effective_sigma_flag = 1)
#   Sample_efficiency = c(Sample_efficiency,Sample_efficiency_tmp)
# }
# op = par(mfrow = c(1,1))
# plot(x = 2:10, y = Sample_efficiency, type = 'o', xlab = 'Number of the Temperature LeveLs',ylab = 'MSD Jumps of Cold Chain (1/T)',main = 'Sampling Efficiency of Parallel Tempering')

###

## report that the num_level should be 6
# nworkers <- detectCores() 
# cl <- makeCluster(nworkers) 
# registerDoParallel(cl)
time_begin_PL = proc.time()
nrep = 1
set.seed(223)
PL_list = list()
length(PL_list) = nrep
#load('./Results_stored/temps_0.1')
#temps =  c(1, 0.99,0.95,0.9,0.85,0.8,0.7,0.6030716, 0.5316445, 0.4652398, 0.3989199, 0.2722091, 0.0717438)
load('./Results_stored/real_temps')
special_factor = 1
element_factor = 1
num_level = length(temps)
for (i_in in 1:nrep) {
  #####################################################
  
  ############################################################
  ## PL
  # op=par(mfrow=c(4,1))
  # temps=c(1,30,50,100)
  # temps = 1/temps
  # for (gamma_use in temps) {
  #   U_fix = U_fix_gam(gamma_use)
  #   y_plot = c()
  #   x_plot =seq(from = -15, to = 15, length.out = 1000)
  #   for (t in x_plot) {
  #     y_plot = c(y_plot, exp(-U_fix(t)))
  #   }
  #   plot(x_plot,y_plot, type = 'l')
  # 
  # }
  # par(op)
  if(density_fun_used == "skew_t"){
    True_data_simulated = rp1(1000000)
    Sigma_mat = covEstimation(True_data_simulated, control = list(type = 'ewma', lambda = 0.9))
    #Init_mean = meanEstimation(True_data_simulated, control = list(type = 'ewma',lambda = 0.9))
    Init_mean = True_data_simulated[1,]
  }else if(density_fun_used == "real_exo"){
    load('stanData')
    True_Data = data
    Sigma_mat = covEstimation(True_Data, control = list(type = 'ewma', lambda = 0.9))
    Init_mean = True_Data[2,]
    Sigma_mat = diag(Sigma_mat)
    Sigma_mat = diag(Sigma_mat)
    Sigma_mat[1,1] = element_factor * Sigma_mat[1,1]
  }
  

  # temps <- 10^( -2 * ((1:num_level) -1)/(num_level-1))
  # temps = (1- max( log10(10^( -1* ((1:num_level) -1)/(num_level-1)) + 1) )) + log10(10^( -1* ((1:num_level) -1)/(num_level-1)) + 1)
  # temps = (-8^(1:num_level))
  # temps = (temps - min(temps))/(max(temps) - min(temps))*(0.9) + 0.1
  # temps = sort(temps,decreasing = TRUE)
  sample_parallel = chains_parallel_tempering(effective_sigma_flag = 1, Sigma_mat = Sigma_mat, init = Init_mean)
  #sample_parallel = chains_parallel_tempering()
  variable_num = 2
  sample_parallel_one_dim = matrix( 0,nrow = iters, ncol = length(temps))
  for (i in 1:iters) {
    sample_parallel_one_dim[ i, ] = sample_parallel[ ,variable_num ,i]
  }
  colnames(sample_parallel_one_dim)=paste("temps=",1/temps,sep="")
  
  ###########################################################
  ###########################################################
  # Our_PL_list = list(distance_each_sample_wass_our=distance_each_sample_wass_our,distance_each_sample_wass_PL=distance_each_sample_wass_PL)
  PL_list[[i_in]] = sample_parallel
}
time_used_PL =proc.time()  - time_begin_PL
# stopCluster(cl)
save(list = c("PL_list","time_used_PL", "num_level","temps"), file =  paste0('./Results_stored/samples_May19_real_Temps_Repeat_num_',nrep,'density_fun_used_',density_fun_used,'num_level_',num_level,'N_',iters,'special_factor',special_factor,'element_factor',element_factor))


