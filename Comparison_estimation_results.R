## Comparison Simulation -- Stochastic Bridge Estimation
rm(list = ls())
# Libraries
library(mvtnorm)
library(mnormt)
library(numDeriv)
library(stats4)
library(MASS)
library(mvtnorm)
library(ggplot2)
## laod the Samples from GWLD algorithm.
#load('./Results_stored/GWL_WithIndex_density_fun_usedskew_trepeatNum_1n1_10000delta_end1e-05stage_Max15proposal_specific1')
#load('./Results_stored/gwl_104_cheat_prop_with_10000_storage_10d_example_run1.RData')
num_level = 1
#load('./Results_stored/PL_samples_Repeat_num_1density_fun_used_skew_tnum_level_13N_5e+05')
#load('./Results_stored/10_1init_1000max_dividPL_samples_Repeat_num_1density_fun_used_skew_tnum_level_13N_50000')
#load('./Results_stored/samples_April_4_firsta_Temps_Repeat_num_1density_fun_used_skew_tnum_level_12N_1e+06special_factor0.1')
#load('./Results_stored/samples_April_6_firsta_Temps_Repeat_num_1density_fun_used_skew_tnum_level_12N_1e+06special_factor0.15')
#load('Results_stored/April_11_PL_List_Final_small_num_sample_animation_com_meanLatest_New_Adaptive_large_sigma50_ourdistance_with_three_distance_N2_com_num_2stage_num_10inser_23density_2N=4000N_Reat1Pre_TRUElarger_than_3')
#load('Results_stored/samples_April_16_real_Temps_Repeat_num_1density_fun_used_real_exonum_level_7N_60000special_factor1')
#load('Results_stored/April_28_narrow_set_start_AdaptiveWarpU_Skewt_N2_com_num_15stage_num_10inser_23density_123N=5000N_Reat1density_fun_used_real_exosigma2_0.8092')
#load('./Results_stored/Initial_Exo_Real_Data_com_num_9stage_num_3inser_23density_2')
#load('Results_stored/May3_narrow_stage_3_20_specialfactor1_Comparisoncom_10Repeat_num1num_level7method_used_PLtime_num_used1density_fun_usedreal_exo')
#load('Results_stored/May5_covest_set_start_AdaptiveWarpU_Skewt_N2_com_num_10stage_num_4inser_23density_123N=10000N_Reat1density_fun_used_real_exosigma2_0.1restricted_flag_FALSE')
#load('Results_stored/standData')
#load('Results_stored/May6_wide_covest_specialfactor1_Comparisoncom_10Repeat_num1num_level1method_used_Adaptive_warpUtime_num_used1density_fun_usedreal_exosigma21')
#load('Results_stored/May19_3_covest_set_start_AdaptiveWarpU_Skewt_N2_com_num_10stage_num_4inser_23density_123N=5000N_Reat1density_fun_used_real_exosigma2_1restricted_flag_FALSEsigma_I_11_factor1diagnal_TRUE')
#load('Results_stored/samples_May19_real_Temps_Repeat_num_1density_fun_used_real_exonum_level_7N_50000special_factor1element_factor1')
load('Results_stored/May19_3_covest_set_start_AdaptiveWarpU_Skewt_N2_com_num_10stage_num_4inser_23density_123N=5000N_Reat1density_fun_used_real_exosigma2_1restricted_flag_FALSEsigma_I_11_factor1diagnal_TRUE')
#load('Results_stored/May19_5_covest_set_start_AdaptiveWarpU_Skewt_N2_com_num_10stage_num_4inser_23density_123N=50000N_Reat1density_fun_used_real_exosigma2_1restricted_flag_FALSEsigma_I_11_factor1diagnal_TRUE')
#load('./Results_stored/AdaptiveWarpU_Skewt_N2_com_num_25stage_num_20inser_23density_123N=30000N_Reat1Pre_TRUEsigma2density_fun_used_skew_t')
# provided data method
#method_used = 'gwl'
#method_used = 'PL'
time_num_used = 1
method_used = 'Adaptive_warpU'
## source the related functions
## skew_t correnspond to density number 123
#density_fun_used = 'skew_t'
#density_fun_used = 'mixture_Normal'
density_fun_used = "real_exo"
narrow_flag = FALSE
density_num = 123
set.seed(23)
## source the GWL algorithm
if(density_fun_used == "skew_t"){
  source('./functions_used/GWLD_skewT.R')
}else{
  source('./functions_used/GWLD.R')
}
source('./functions_used/Stochastics_functions.R')
source('./functions_used/warpu code/EMhd_diagonal.R')
source('./functions_used/warpu code/Bridge_Sampling_funs_eff.R')
if(density_fun_used == "skew_t"){
  source('./functions_used/skew_t_density.R')
}
#####
#total_sample_num = dim(est_result$x_final_stage)[1]
if(density_fun_used == "skew_t"){
  true_Z = 1
}else if(density_fun_used == 'mixture_Normal'){
  true_Z = 2*pi
}else if(density_fun_used == 'real_exo'){
  #true_Z = -193.68
  true_Z = -193.71
}
RepeatNum = 1
### suppose we have already got the GWL samples.
K_prop = com_num = 10
if(density_fun_used == "skew_t"){
  Dim = 10
}else if(density_fun_used == 'mixture_Normal'){
  Dim = 2
  p1 = function(x,log=F){
    nn = dim(x)[1]
    return_results = c()
    for (i in 1:nn) {
      return_results = c(return_results, q_density(x[i,]))
      
    }
   
    if(log){
      return(log(return_results))
    }else{
      return(return_results)
    }
  }
}else if(density_fun_used == 'real_exo'){
  Dim = 7
  load('./set posterior function/dataset1_one_planet_posterior.RData')
  if(narrow_flag == FALSE){
    source('./R model code/full_priors.R')
    Pmin <- rep(1.25,3)
    Pmax <- rep(10^4,3)
  }
  source("./R model code/log_like.R")
  source("./R model code/log_post.R")
  source("./R model code/cov_fun.R")
  source("./R model code/cov_make.R")
  source("./R model code/planet_model.R") # Added in sqrt to phi
  source('./functions_used/Functions_AdaptiveWarpU.R')
  true_Z = -193.71
  matrix_orig = c()
  matrix_Warp = c()
  matrix_Stoc = c()
}

p1
if(method_used == "PL"){
  num_level = dim(PL_list[[1]])[1]
  ## condition
  # sample_condition = t(PL_list[[1]][1,,])
  # max_condition = apply(sample_condition, 1, max)
  # min_condition = apply(sample_condition, 1, min)
  # condition_index = (max_condition<=20)&(min_condition>=-20)
  # sample_condition = sample_condition[condition_index,]
  # aa = (dim(sample_condition)[1]/2):(dim(sample_condition)[1])
  # index_PL_used = sample(aa,25937)
  # sample_results = sample_condition[ index_PL_used, ]
  ##
  aa = floor((dim(PL_list[[1]])[3]/3)) : ( dim(PL_list[[1]])[3])
  #index_PL_used = sort(1:dim(PL_list[[1]])[3],decreasing = TRUE)[1:25937]
  #index_PL_used  = sample(1:dim(PL_list[[1]])[3],25937)
  index_PL_used=  floor(seq(from = min(aa), to = max(aa), length.out = 5000))
  #index_PL_used  = sample(aa,25937)
  index_PL_used = 15715:31428
  sample_results = t(PL_list[[1]][ time_num_used , , index_PL_used  ])
  # sample_results = sample_results[1:25937, ]
  hist(sample_results[,6])
  plot(density(sample_results[,6]))
  
}else if(method_used == 'Adaptive_warpU'){
  Index_set = seq(from = 1, to= N, by = 1)
  #Index_set = sample(1:10000, size = 5000)
  sample_results = results_sample$w_all[Index_set,]
  #hist(sample_results[,6])
}

if(method_used == 'gwl'){
  mean_zhat = c()
  for (rept_num in 1:RepeatNum) {
    mean_zhat = c(mean_zhat,Repeat_store_zhat[[rept_num]][[11]])
  }
  mean_zhat = mean(mean_zhat)
  GWL_estimator = mean_zhat
  GWL_estimator_RMSE = sqrt((GWL_estimator - true_Z)^2)
  print('this is RMSE of GWL')
  GWL_estimator_RMSE
  
  length_record = c()
  for (i in 1:length(x_subregion_store)) {
    length_record = c(length_record, dim(x_subregion_store[[i]])[1])
    
  }
  sum(length_record)
  x_subregion_store
  #number_last_stage_target_evaluation = 2* (sum(length_record)  -  (dim(x_subregion_store[[1]])[1]-1)) +  (dim(x_subregion_store[[1]])[1]-1)
  number_last_stage_target_evaluation = n_store[[11]]
  
}else if(method_used  == 'PL'){
  number_last_stage_target_evaluation = dim(sample_results)[1]
  sampleNum_by_GWL = dim(sample_results)[1]
}else if (method_used == 'Adaptive_warpU'){
  number_last_stage_target_evaluation = dim(sample_results)[1]
  sampleNum_by_GWL = dim(sample_results)[1]
}
if(density_fun_used == 'real_exo'){
  #Cost_est = 49685
  Cost_est = 283000
  N1_in = number_last_stage_target_evaluation
  number_last_stage_target_evaluation = Cost_est
}


### get the samples from GWL results
#scale_grid = c(-1.5,-1,-0.5,0,0.5,1)
#number_grid = floor ( ((10^scale_grid) * number_last_stage_target_evaluation ) / com_num )
#sampleNum_by_GWL = floor(number_last_stage_target_evaluation/com_num)
#sample_by_GWL =  GetSamples_from_LastStage(sample_size = sampleNum_by_GWL, estimate_result_GWL = est_result,david_data = 1)
RMSE_orig_mat = c()
RMSE_WarpBs_mat = c()
RMSE_StochBs_mat = c()
time_com_begin = proc.time()
for (repNum in 1:RepeatNum) {
  set.seed(repNum+23)
  RMSE_orig = c()
  RMSE_WarpBs = c()
  RMSE_StochBs = c()
  Est_orig = c()
  Est_WarpBs = c()
  Est_StochBs = c()
  #number_grid = seq(from = 5000, to = 20000, by = 50)
  #scale_grid = c(-1,-0.5,0,0.5)
  scale_grid = c(-1,-0.5,0,0.5,1)
  scale_grid = c(0)
  number_grid = floor ( ((10^scale_grid) * number_last_stage_target_evaluation ) / (2*com_num-1) )
  if(method_used == 'gwl'){
    log_g_estimates = Repeat_store_log_g_estimate[[repNum]]
    zhat_store = Repeat_store_zhat[[repNum]]
    x_subregion_store = Repeat_store_x_subregion_store[[repNum]]
  }
  
  for (num_sample in (2*number_grid)) {
    sample_used = num_sample
    print('timebegin')
    time_beggin = proc.time()
    print(time_beggin)
    if(method_used == 'gwl'){
      all_GWLSample_used =  GetSamples_from_LastStage(sample_size = sample_used, estimate_result_GWL = est_result,david_data = 1)
    }else{
      #index_used = sample(x = (1:sampleNum_by_GWL), size = (((2*com_num-1) ) *sample_used) ,replace = TRUE)
      index_used = sample(x = (1:sampleNum_by_GWL), size = (sample_used) ,replace = TRUE)
      all_GWLSample_used = sample_results[index_used,]
    }
    print('timeused')
    print(proc.time()-time_beggin)
    
    # divide the data into two parts
    fit_mat = all_GWLSample_used
    simulation_sample_size = dim(fit_mat)[1] 
    n = simulation_sample_size
    first_half = sample(1:n,size=floor(n/2),replace=FALSE)# 1:(n/2)
    second_half = setdiff(1:n,first_half) #(n/2+1):n
    data = fit_mat
    data_p1_first_half = data[first_half,]
    data_p1_second_half = data[second_half,]
    
    # First run EM
    sample_size = min(K_prop*500,floor(simulation_sample_size/2))
    index1 = sample(x=first_half,size=sample_size)
    x = data[index1,]
    tem_em = EM_rep(x,K_prop,reps=1,niters=1000)
    total_iter_runs <- tem_em$EM_lists[[1]]$iter+tem_em$EM_lists[[2]]$iter
    proposal_evaluations <- sample_size*2*K_prop*(total_iter_runs-1)
    prop_prop = tem_em$weight
    valid_index = which(prop_prop!=0)
    
    # get the parameters from EM 
    prop_prop = prop_prop[valid_index]
    K_prop_original <- K_prop
    K_prop = length(prop_prop)
    mu_prop = matrix(matrix(tem_em$mu,K_prop_original,Dim)[valid_index,],K_prop,Dim)
    sigma2_prop = matrix(matrix(tem_em$sigma2,K_prop_original,Dim)[valid_index,],K_prop,Dim)
    
    # transfer the data from the target density
    #subsample_ind = sample(x = (1:dim(data_p1_second_half)[1]), size =  floor(dim(data_p1_second_half)[1]/(2*K_prop-1)) )
    if(density_fun_used == 'real_exo'){
      data_p1_second_half = sample_results
      data_tp1 = rtp1(sample_results)
    }else{
      data_tp1 = rtp1(data_p1_second_half)
    }
      
    proposal_evaluations <- proposal_evaluations+dim(data_p1_second_half)[1]*K_prop*(K_prop+1)
    
    ##########
    # specify the number of samples from normal
    if(density_fun_used == 'real_exo'){
      
      if(method_used == 'Adaptive_warpU'){
        # N2 = floor(n/2) * (2*com_num-1)
        N2 = floor(n/2*(2*com_num-1)/com_num)
      }else{
        N2 = floor( (n/2*(2*com_num-1) - (com_num-1)*N1_in) /com_num)
      }
      
    }else{
      if(method_used == 'Adaptive_warpU'){
        # N2 = floor(n/2) * (2*com_num-1)
        N2 = floor(n/2*(2*com_num-1)/com_num)
      }else{
        N2 = floor(n/2)
      }
    }
    
    
    #N2 = floor(N2/(2*K_prop-1))
    data_normal = rmvnorm(N2,mean=rep(0,Dim),sigma=diag(1,Dim))
    # apply the optimal bridge sampling on data_tp1 and data_normal
    # (n/2+N2)K target evaluations
    # (n/2+N2)K^2 Normal evaluations
    # 2(N2 + (n/2)) standard Normal evaluations
    proposal_evaluations <- proposal_evaluations+(dim(data_p1_second_half)[1]+N2)*K_prop^2 
    target_evaluations <- (dim(data_p1_second_half)[1]+N2)*K_prop
    standard_normal_evaluations <- (dim(data_p1_second_half)[1]+N2)*2
    
    ###
    ### get the estimator of original Bridge Sampling method
    #data_phimix =rphi_mix(floor(n/2))
    data_phimix =rphi_mix(floor(n/2*(2*com_num-1)))
    log_bs_Z_hat_original_Bs <- Opt_BS_hd(x1=data_p1_second_half,x2=data_phimix, q1=p1,q2=phi_mix,r0=10)
    print('this is the estimator of original Bridge method ')
    exp(log_bs_Z_hat_original_Bs)
    if(density_fun_used == 'real_exo'){
      RMSE_orig =c(RMSE_orig, ((log10(exp(log_bs_Z_hat_original_Bs)) - true_Z)^2))
      Est_orig = c(Est_orig, log10(exp(log_bs_Z_hat_original_Bs)))
    }else{
      RMSE_orig =c(RMSE_orig, ((exp(log_bs_Z_hat_original_Bs) - true_Z)^2))
    }
    
    ###
    ### get the estimator of warpU Bridge Sampling method
    log_bs_Z_hat_WarpU_Bs <- Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1,q2=p2_warp0,r0=10)
    print('this is the estimator of WarpU Bridge method ')
    exp(log_bs_Z_hat_WarpU_Bs)
    if(density_fun_used == 'real_exo'){
      RMSE_WarpBs =c(RMSE_WarpBs, ((log10(exp(log_bs_Z_hat_WarpU_Bs)) - true_Z)^2))
      Est_WarpBs = c(Est_WarpBs, log10(exp(log_bs_Z_hat_WarpU_Bs)))
    }else{
      RMSE_WarpBs =c(RMSE_WarpBs, ((exp(log_bs_Z_hat_WarpU_Bs) - true_Z)^2))
      
    }
    
    ###
    ### get the estimator of Stochastics Bridge Sampling method
    
    ## one way to use the normal data as the same size with tilde{w_k}
    # split normal distribution
    if(density_fun_used == 'real_exo'){
      
      Opt_Stochas_BS_splitNorm_same = function(w_in){
        tildeW_k  = rtp1_split(w_in)
        split_norm = list()
        length(split_norm) = K_prop
        for (kk in 1:K_prop) {
          if(ceiling(length(tildeW_k[[kk]])) > 0){
            split_norm[[kk]] = rmvnorm( (( floor(Cost_est/N1_in) )*ceiling(length(tildeW_k[[kk]]))),mean=rep(0,Dim),sigma=diag(1,Dim))
          }
          
        }
        return(Opt_BS_hd_Sto_splitNorm(split_data = tildeW_k,split_norm = split_norm))
      }
      
    }else{
      
      Opt_Stochas_BS_splitNorm_same = function(w_in){
        tildeW_k  = rtp1_split(w_in)
        split_norm = list()
        length(split_norm) = K_prop
        for (kk in 1:K_prop) {
          if(ceiling(length(tildeW_k[[kk]])) > 0){
            split_norm[[kk]] = rmvnorm( ((2*com_num-1)*ceiling(length(tildeW_k[[kk]]))),mean=rep(0,Dim),sigma=diag(1,Dim))
          }
          
        }
        return(Opt_BS_hd_Sto_splitNorm(split_data = tildeW_k,split_norm = split_norm))
      }
    }
    
    
    
    ## excute
    #bs_Z_hat_StochasticWarp_Bs = try(Opt_Stochas_BS_splitNorm_same(data_p1_second_half), silent = TRUE)
    
    for(try_num in 1:100){
      try({
        bs_Z_hat_StochasticWarp_Bs = Opt_Stochas_BS_splitNorm_same(data_p1_second_half)
        break #break/exit the for-loop
      }, silent = FALSE)
    }
    
    
    if (try_num == 100){
      if(density_fun_used == 'real_exo'){
        
        Opt_Stochas_BS_splitNorm_same = function(w_in){
          tildeW_k  = rtp1_split(w_in)
          split_norm = list()
          length(split_norm) = K_prop
          for (kk in 1:K_prop) {
            if(ceiling(length(tildeW_k[[kk]])) > 0){
              split_norm[[kk]] = rmvnorm( (( floor(Cost_est/N1_in) )*ceiling(length(tildeW_k[[kk]]))),mean=rep(0,Dim),sigma=diag(1,Dim))
            }
            
          }
          return(Opt_BS_hd_Sto_splitNorm(split_data = tildeW_k,split_norm = split_norm))
        }
        
      }else{
        Opt_Stochas_BS_splitNorm_same = function(w_in){
          tildeW_k  = rtp1_split(w_in)
          split_norm = list()
          length(split_norm) = K_prop
          for (kk in 1:K_prop) {
            if(ceiling(length(tildeW_k[[kk]])) > 0){
              split_norm[[kk]] = rmvnorm( ((2*com_num-1)*ceiling(length(tildeW_k[[kk]]))),mean=rep(0,Dim),sigma=diag(1,Dim))
            }
            
          }
          return(Opt_BS_hd_Sto_splitNorm_Error(split_data = tildeW_k,split_norm = split_norm))
        }
      }
      
      
      
      bs_Z_hat_StochasticWarp_Bs = Opt_Stochas_BS_splitNorm_same(data_p1_second_half)
      if(density_fun_used == 'real_exo'){
        RMSE_StochBs =c(RMSE_StochBs, ((log10((bs_Z_hat_StochasticWarp_Bs)) - true_Z)^2))  
        Est_StochBs = c(Est_StochBs,log10((bs_Z_hat_StochasticWarp_Bs)) )
      }else{
        RMSE_StochBs =c(RMSE_StochBs, (((bs_Z_hat_StochasticWarp_Bs) - true_Z)^2))  
        
      }
    }else{ 
      print('this is the estimator of stochastics WarpU Bridge method ')
      #exp(log_bs_Z_hat_StochasticWarp_Bs)
      if(density_fun_used == 'real_exo'){
        RMSE_StochBs =c(RMSE_StochBs, ((log10((bs_Z_hat_StochasticWarp_Bs)) - true_Z)^2))  
        Est_StochBs = c(Est_StochBs,log10((bs_Z_hat_StochasticWarp_Bs)) )
        
      }else{
        RMSE_StochBs =c(RMSE_StochBs, (((bs_Z_hat_StochasticWarp_Bs) - true_Z)^2))  
        
      }
    }
    
    
  }
  RMSE_orig_mat = rbind(RMSE_orig_mat, (RMSE_orig))
  RMSE_WarpBs_mat = rbind(RMSE_WarpBs_mat, (RMSE_WarpBs))
  RMSE_StochBs_mat = rbind(RMSE_StochBs_mat, (RMSE_StochBs))
  if(density_fun_used == 'real_exo'){
    matrix_orig = rbind(matrix_orig,Est_orig)
    matrix_Warp = rbind(matrix_Warp, Est_WarpBs)
    matrix_Stoc = rbind(matrix_Stoc, Est_StochBs)
  }
  
  
}

RMSE_orig_mat1 = sqrt(colMeans(RMSE_orig_mat))
RMSE_WarpBs_mat1 = sqrt(colMeans(RMSE_WarpBs_mat))
RMSE_StochBs_mat1= sqrt(colMeans(RMSE_StochBs_mat))
RMSE_orig_mat2 = (apply(RMSE_orig_mat, 2, sd))
RMSE_WarpBs_mat2 = (apply(RMSE_WarpBs_mat, 2, sd))
RMSE_StochBs_mat2 = (apply(RMSE_StochBs_mat, 2, sd))
time_com_used = proc.time() -  time_com_begin 

plot(RMSE_orig)
plot(RMSE_WarpBs)
plot(RMSE_StochBs)

#x_seq = c(-1,-0.5,0,0.5)
x_seq = c(-1,-0.5,0,0.5,1)
x_seq = c(0)


# RMSE_orig_data = data.frame(x_seq = x_seq[2:5], RMSE =  RMSE_orig_mat1[2:5], se = RMSE_orig_mat2[2:5]/sqrt(RepeatNum)*(1/2*(colMeans(RMSE_orig_mat))^(-1/2))[2:5])
# RMSE_WarpBs_data = data.frame(x_seq = x_seq[2:5], RMSE =  RMSE_WarpBs_mat1[2:5], se = RMSE_WarpBs_mat2[2:5]/sqrt(RepeatNum)*(1/2*(colMeans(RMSE_WarpBs_mat))^(-1/2))[2:5])
# RMSE_StochBs_data = data.frame(x_seq = x_seq[2:5], RMSE =  RMSE_StochBs_mat1[2:5], se = RMSE_StochBs_mat2[2:5]/sqrt(RepeatNum)*(1/2*(colMeans(RMSE_StochBs_mat))^(-1/2))[2:5])
RMSE_orig_data = data.frame(x_seq = x_seq, RMSE =  RMSE_orig_mat1, se = RMSE_orig_mat2/sqrt(RepeatNum)*(1/2*(colMeans(RMSE_orig_mat))^(-1/2)))
RMSE_WarpBs_data = data.frame(x_seq = x_seq, RMSE =  RMSE_WarpBs_mat1, se = RMSE_WarpBs_mat2/sqrt(RepeatNum)*(1/2*(colMeans(RMSE_WarpBs_mat))^(-1/2)))
RMSE_StochBs_data = data.frame(x_seq = x_seq, RMSE =  RMSE_StochBs_mat1, se = RMSE_StochBs_mat2/sqrt(RepeatNum)*(1/2*(colMeans(RMSE_StochBs_mat))^(-1/2)))

if(method_used=='gwl'){
  RMSE_gwl_data = data.frame(x_seq = x_seq, RMSE =  GWL_estimator_RMSE)
}


p = ggplot(RMSE_orig_data, aes(x_seq, RMSE)) +
  geom_line(aes(color = 'Original Bridge'))+ geom_point(aes(color = 'Original Bridge'),shape = 2) +
  geom_errorbar(aes(x = x_seq, y = RMSE,ymin= RMSE-2*se,ymax=RMSE+2*se,color = 'Original Bridge'),alpha=I(3/7),width=.1,
                position=position_dodge(0.05)) +
  geom_line(data = RMSE_WarpBs_data,aes(x_seq, RMSE,color = 'WarpUBridge'))+geom_point(data = RMSE_WarpBs_data,aes(x_seq, RMSE,color = 'WarpUBridge'),shape = 23)+
  geom_errorbar(data = RMSE_WarpBs_data,aes(x = x_seq, y = RMSE,ymin= RMSE-2*se,ymax=RMSE+2*se,color = 'WarpUBridge'),alpha=I(3/7),width=.1,
                position=position_dodge(0.05)) +
  geom_line(data = RMSE_StochBs_data,aes(x_seq, RMSE,color = 'StochasticsBridge'))+geom_point(data = RMSE_StochBs_data,aes(x_seq, RMSE,color = 'StochasticsBridge'),shape = 21) +
  geom_errorbar(data = RMSE_StochBs_data,aes(x = x_seq, y = RMSE,ymin= RMSE-2*se,ymax=RMSE+2*se,color = 'StochasticsBridge'),alpha=I(3/7),width=.1,
                position=position_dodge(0.05)) 
if(method_used=='gwl'){
  p = p+ geom_line(data = RMSE_gwl_data,aes(x_seq, RMSE,color = 'gwl'))
}
p = p+ theme_bw() + ggtitle('Comparation between Bridge Sampling Methods') +  ylab("RMSE") + xlab('The Number of Evaluation') +
  labs(colour="Method")+
  theme(plot.title = element_text(hjust = 0.5,size = 15))+
  theme(legend.position = c(0.82, 0.81),legend.justification = c(0, 1))+
  scale_color_manual(values = c("green","blue","red","orange"))
p
if(method_used  == 'PL'){
  rm(PL_list)
}

#save(list = ls(), file = paste0('./Results_stored/May6_standsample','com_',com_num,'Repeat_num',RepeatNum,'num_level',num_level,'method_used_',method_used,'time_num_used',time_num_used,'density_fun_used',density_fun_used))


save(list = ls(), file = paste0('./Results_stored/May22_cost_large_Comparison','com_',com_num,'Repeat_num',RepeatNum,'num_level',num_level,'method_used_',method_used,'time_num_used',time_num_used,'density_fun_used',density_fun_used))



##
# library(coda)
# effectiveSize(PL_list[[1]][ , 1,])
# library(smfsb)
# variable_num = 2
# sample_parallel_one_dim = matrix( 0,nrow = 1000000, ncol = length(temps))
# for (i in 1:1000000) {
#   sample_parallel_one_dim[ i, ] = PL_list[[1]][ ,variable_num ,i]
# }
# colnames(sample_parallel_one_dim)=paste("temps=",1/temps,sep="")
# sample_parallel_one_dim[,1:5]
# png(file = 'abc.pdf')
# (mcmcSummary(sample_parallel_one_dim[,1:5],rows=length(temps[1:5])))
# dev.off()

# ## bad data
# plot(density(True_Data[,2],bw = 0.8017), col = 'red', main = '', ylim = c(0, 0.15))
# lines(density(results_sample$w_all[,2],bw = 0.8017 ), col = 'blue', lty = 2)
# 
# legend( x = 10,  y = 0.10, legend=c("True Density", "Bad Samples"), col=c("red", "blue"), lty=c(1:2), cex=0.8, box.lty=0)
# #