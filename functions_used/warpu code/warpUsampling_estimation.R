rm(list=ls())

# Get posterior (including data)
#load(paste("dataset",dataset_num,"_one_planet_posterior.RData",sep=""))
#load("../set posterior function/dataset1_no_planet_posterior.RData")
load("../set posterior function/dataset1_one_planet_posterior.RData")
load('../warpUdraws/sample_stored/List_Final_small_num_sample_animation_com_meanLatest_New_Adaptive_large_sigma50_ourdistance_with_three_distance_N2_com_num_10stage_num_10inser_23density_2N=40000N_Reat1Pre_TRUElarger_than_3')
source('../warpUdraws/functionused/Functionsexo_new.R')
source('../warpUdraws/functionused/simple_case_density.R')

K_prop = com_num
Dim = 2
p1 = q_density_est
com_approx = 5

# Libraries
library(mvtnorm)
library(mnormt)
library(numDeriv)
library(stats4)
library(MASS)
library(mvtnorm)

# Get Warp-U code

source("EMhd_diagonal.R")
source("Bridge_Sampling_funs_eff.R")

# Evaluation costs:
# 1) Posterior draws - this we are taken as "given" i.e. the setting we consider is where posterior draws 
#                      have alrady been collected
# 2) EM - Normals with diagonal covariance matrices 
# 3) Warp-U: relatively few proposal evaluations and some target evaluations

# Split draws
fit_mat = results_sample$w_all
simulation_sample_size = dim(fit_mat)[1] # Loaded stan HMC samples
n = simulation_sample_size
first_half = sample(1:n,size=n/2,replace=FALSE)# 1:(n/2)
second_half = setdiff(1:n,first_half) #(n/2+1):n
data = fit_mat

data_p1_first_half = data[first_half,]
data_p1_second_half = data[second_half,]

# First run EM
sample_size = min(K_prop*50,simulation_sample_size/2)
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

if (1==2){
  for (i in 1:(dim(data)[2]-1)){
    for (j in (i+1):(dim(data)[2])){
      plot(data[,c(i,j)])
      points(mu_prop[,c(i,j)],col=2,pch=16,cex=2)
      readline("pause")
    }
  }
}


# transform the data from the target density
data_tp1 = rtp1(data_p1_second_half)  
proposal_evaluations <- proposal_evaluations+dim(data_p1_second_half)[1]*K_prop*(K_prop+1)



##########
##########
##########
# specify the number of samples from normal
N2 = 10000
data_normal = rmvnorm(N2/2,mean=rep(0,Dim),sigma=diag(1,Dim))
# apply the optimal bridge sampling on data_tp1 and data_normal
# (n/2+N2)K target evaluations
# (n/2+N2)K^2 Normal evaluations
# 2(N2 + (n/2)) standard Normal evaluations
proposal_evaluations <- proposal_evaluations+(dim(data_p1_second_half)[1]+N2/2)*K_prop^2 
target_evaluations <- (dim(data_p1_second_half)[1]+N2/2)*K_prop
standard_normal_evaluations <- (dim(data_p1_second_half)[1]+N2/2)*2
log_bs_Z_hat <- Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1_variant,q2=p2_warp0,r0=10)
log_bs_Z_hat <- Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1_variant_1,q2=p2_warp0,r0=10)
log_bs_Z_hat <- Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1_variant_2,q2=p2_warp0,r0=10)

log_bs_Z_hat <- Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1,q2=p2_warp0,r0=10)





# calculate the MSE -- using boostra
boostrape_num = 10
MSE_1_matrix = matrix(0,nrow = K_prop,boostrape_num)
MSE_2_matrix = matrix(0,nrow = K_prop,boostrape_num)
for (i in 1:boostrape_num) {
  MSE_1_com = rep(0,K_prop)
  MSE_2_com = rep(0,K_prop)
  for(com_approx in 1:K_prop){
    MSE_variant_1 = MSE_tp1_variant_1(data_normal,log = T)
    MSE_variant_2 = MSE_tp1_variant_2(data_normal,log =T)
    MSE_1_com[com_approx]=MSE_variant_1
    MSE_2_com[com_approx]=MSE_variant_2
  }   
  MSE_1_matrix[,i] = MSE_1_com
  MSE_2_matrix[,i] = MSE_2_com
}
mean_MSE_1 = rowMeans(MSE_1_matrix)
mean_MSE_2 = rowMeans(MSE_2_matrix)
#mean_MSE_ori = MSE_tp1(data_normal,log = T)

plot(mean_MSE_1)
#lines(x = 1:K_prop, y = rep(min(),K_prop), col = 'red')
plot(mean_MSE_2)
lines(x = 1:K_prop, y = rep(min(mean_MSE_1),K_prop), col = 'red')



com_approx = 5
## variant, variant1, variant2
est_variant = rep(0,boostrape_num)
est_variant1 = rep(0,boostrape_num)
est_variant2= rep(0,boostrape_num)
for (i in 1:boostrape_num) {
  est_variant[i] <- Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1_variant,q2=p2_warp0,r0=10)
  est_variant1[i] <- Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1_variant_1,q2=p2_warp0,r0=10)
  est_variant2[i] <- Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1_variant_2,q2=p2_warp0,r0=10)
  
}
mean(exp(est_variant))
mean(exp(est_variant1))
mean(exp(est_variant2))


## time compared
library(microbenchmark)

microbenchmark(
  Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1_variant,q2=p2_warp0,r0=10),
  Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1,q2=p2_warp0,r0=10),
  
  times = 3
)
log10_bs_Z_hat <- log_bs_Z_hat/log(10)
Z_hat_store <- exp(log_bs_Z_hat)
log10_bs_Z_hat



## R = r1*pi + r2*p2+...
## split data
split_data = split_data_k(x_in = data_tp1)
aa = Opt_BS_hd_Sto(split_data)

microbenchmark(
  Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1_variant,q2=p2_warp0,r0=10),
  Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1,q2=p2_warp0,r0=10),
  Opt_BS_hd_Sto(split_data),
  times = 3
)


##
## start with w to tilde{w}_k
tildeW_k  = rtp1_split(data_p1_second_half)
cc = Opt_BS_hd_Sto(tildeW_k)


## time compare 
method0 = function(w_in){
  data_trans = rtp1(w_in) 
  return(exp(Opt_BS_hd(x1=data_trans,x2=data_normal, q1=tp1,q2=p2_warp0,r0=10)))
}
method1 = function(w_in){
  data_tp1 = rtp1(w_in)
  split_data = split_data_k(x_in = data_tp1)
  return(Opt_BS_hd_Sto(split_data))
}

# using the same normal distribution
method2 = function(w_in){
  tildeW_k  = rtp1_split(w_in)
  return(Opt_BS_hd_Sto(tildeW_k))
}


# method original 
method_original = function(w_in, normal_in){
  data_trans = rtp1(w_in) 
  return(exp(Opt_BS_hd(x1=data_trans,x2= normal_in, q1=tp1,q2=p2_warp0,r0=10)))
}
# using different normal distribution
method_w_Diff_normal_in = function(w_in, normal_in){
  tildeW_k  = rtp1_split(w_in)
  return(Opt_BS_hd_Sto_splitNorm(split_data = tildeW_k,split_norm = normal_in))
}

# split normal distribution
method3 = function(w_in){
  tildeW_k  = rtp1_split(w_in)
  split_norm = list()
  length(split_norm) = K_prop
  for (kk in 1:K_prop) {
    split_norm[[kk]] = rmvnorm(length(tildeW_k[[kk]]),mean=rep(0,Dim),sigma=diag(1,Dim))
    
  }
  
  return(Opt_BS_hd_Sto_splitNorm(split_data = tildeW_k,split_norm = split_norm))
}



microbenchmark(
  method1(data_p1_second_half),
  method2(data_p1_second_half),
  times = 3
)
library(profvis)
profvis({
  method0(data_p1_second_half)
})
profvis({
  method2(data_p1_second_half)
})

method0(data_p1_second_half)
method1(data_p1_second_half)
method2(data_p1_second_half)


## time comparation
microbenchmark(
  method0(data_p1_second_half),
  method1(data_p1_second_half),
  method2(data_p1_second_half),
  times = 3
)




### compare the variance
boostrape_num =5
N2_vec = seq(from = 5000, to = 29000, by = 2000)
N2_len = length(N2_vec)
r0_mean_vec = rep(0,N2_len)
r2_mean_vec = rep(0,N2_len)
r3_mean_vec = rep(0,N2_len)
r4_mean_vec = rep(0,N2_len)
r0_sd_vec = rep(0,N2_len)
r2_sd_vec = rep(0,N2_len)
r3_sd_vec = rep(0,N2_len)
r4_sd_vec = rep(0,N2_len)
r0_boos = rep(0,boostrape_num)
r2_boos = rep(0,boostrape_num)
r3_boos = rep(0,boostrape_num)
r4_boos = rep(0,boostrape_num)
for ( i_in in 1:N2_len) {
  N2_in = N2_vec[i_in]
  data_normal_all = rmvnorm(N2_in,mean=rep(0,Dim),sigma=diag(1,Dim))
  for (j_in in 1:boostrape_num) {
    data_normal = data_normal_all
    r0_boos[j_in] = method0(data_p1_second_half)
    r2_boos[j_in] = method2(data_p1_second_half)
    r3_boos[j_in] = method3(data_p1_second_half)
    data_normal = data_normal_all[1:floor(N2_in/K_prop), ]
    r4_boos[j_in] = method2(data_p1_second_half)
  }
  r0_mean_vec[i_in] =  mean(r0_boos)
  r2_mean_vec[i_in] =  mean(r2_boos)
  r3_mean_vec[i_in] =  mean(r3_boos)
  r4_mean_vec[i_in] =  mean(r4_boos)
  r0_sd_vec[i_in] = sd(r0_boos)
  r2_sd_vec[i_in] = sd(r2_boos)
  r3_sd_vec[i_in] = sd(r3_boos)
  r4_sd_vec[i_in] = sd(r4_boos)
  
  
}


###
###
error_bar_num = 10
boostrap_num <- 10
#N2_vec <- c(2000,3000)
N2_vec = seq(from = 5000, to = 31000, by = 2000)
Bias2_mat1 = matrix(0, nrow = length(N2_vec), ncol = error_bar_num)
var_mat1 = matrix(0, nrow = length(N2_vec), ncol = error_bar_num)
MSE_mat1 = matrix(0, nrow = length(N2_vec), ncol = error_bar_num)
Bias2_mat2 = matrix(0, nrow = length(N2_vec), ncol = error_bar_num)
var_mat2 = matrix(0, nrow = length(N2_vec), ncol = error_bar_num)
MSE_mat2 = matrix(0, nrow = length(N2_vec), ncol = error_bar_num)
for (err_times in 1:error_bar_num) {
  Results_comppare <- compare_method(boostrap_num = boostrap_num, N2_vec = N2_vec , w_in = data_p1_second_half)
  Bias2_mat1[ , err_times] <- Results_comppare$r1_Bias2
  var_mat1[ , err_times] <- Results_comppare$r1_Var
  MSE_mat1[ , err_times] <- Results_comppare$r1_MSE
  Bias2_mat2[ , err_times] <- Results_comppare$r2_Bias2
  var_mat2[ , err_times] <- Results_comppare$r2_Var
  MSE_mat2[ , err_times] <- Results_comppare$r2_MSE
}

list_result = list()
list_result[[1]] = rowMeans(var_mat1)
list_result[[2]] = rowMeans(var_mat2)
list_result[[3]] = rowMeans(Bias2_mat1)
list_result[[4]] = rowMeans(Bias2_mat2)
list_result[[5]] = rowMeans(MSE_mat1)
list_result[[6]] = rowMeans(MSE_mat2)

list_sd = list()
list_sd[[1]] = apply(var_mat1, 1, sd)
list_sd[[2]] = apply(var_mat2, 1, sd)
list_sd[[3]] = apply(Bias2_mat1, 1, sd)
list_sd[[4]] = apply(Bias2_mat2, 1, sd)
list_sd[[5]] = apply(MSE_mat1, 1, sd)
list_sd[[6]] = apply(MSE_mat2, 1, sd)

name_str = c('var', 'bias2', 'MSE')
for (i in 1:(length(list_result)/2)) {
  miny <- min(list_result[[2*i-1]] - list_sd[[2*i-1]], list_result[[2*i]] -list_sd[[2*i]] )
  maxy <- max(list_result[[2*i-1]] + list_sd[[2*i-1]], list_result[[2*i]] + list_sd[[2*i-1]])
  plot(x = N2_vec, y = list_result[[2*i-1]], col = 'blue', pch = 20,ylim = c(miny, maxy), ylab = name_str[i])
  arrows(x0 = N2_vec, y0 = list_result[[2*i-1]]-list_sd[[2*i-1]], x1 = N2_vec, y1 = list_result[[2*i-1]]+list_sd[[2*i-1]], length=0.05, angle=90, code=3, col = 'blue')
  points(x = N2_vec, y = list_result[[2*i]], col = 'red', pch = 15)
  arrows(x0 = N2_vec, y0 = list_result[[2*i]]-list_sd[[2*i]], x1 = N2_vec, y1 = list_result[[2*i]]+list_sd[[2*i]], length=0.05, angle=90, code=3, col = 'red')
  
  
}



#N2_vec <- c(2000,3000)
# N2_vec = seq(from = 5000, to = 31000, by = 2000)
# Results_comppare <- compare_method(boostrap_num = boostrap_num, N2_vec = N2_vec , w_in = data_p1_second_half)
name_str = c('mean', 'var', 'bias2', 'MSE')
for (i in 1:(length(Results_comppare)/2)) {
  miny <- min(Results_comppare[[2*i-1]], Results_comppare[[2*i]])
  maxy <- max(Results_comppare[[2*i-1]], Results_comppare[[2*i]])
  plot(x = N2_vec, y = Results_comppare[[2*i-1]], col = 'blue', pch = 20,ylim = c(miny, maxy), ylab = name_str[i])
  points(x = N2_vec, y = Results_comppare[[2*i]], col = 'red', pch = 15)
  
}
miny <- min(Results_comppare$r1_MSE, Results_comppare$r2_MSE)
maxy <- max(Results_comppare$r1_MSE, Results_comppare$r2_MSE)
plot(x = N2_vec, y = Results_comppare$r1_MSE, col = 'blue', pch = 20,ylim = c(miny, maxy))

points(x = N2_vec, y = Results_comppare$r2_MSE, col = 'red', pch = 15)

###
###
########
########
###plot##
## mean
plot(x = N2_vec, y = abs((r0_mean_vec-2*pi)), col = 'blue', pch = 20)
plot(x = N2_vec, y = abs((r2_mean_vec-2*pi)), col = 'red', pch = 15)
plot(x = N2_vec, y = abs((r3_mean_vec-2*pi)), col = 'purple', pch = 12)
plot(x = N2_vec, y = abs((r4_mean_vec-2*pi)), col = 'orange', pch = 12)

miny = min(min(abs((r2_mean_vec-2*pi))), min(abs((r0_mean_vec-2*pi))))
miny = min(miny,min(abs((r3_mean_vec-2*pi))))
miny = min(miny,min(abs((r4_mean_vec-2*pi))))
maxy = max(max(abs((r2_mean_vec-2*pi))),max(abs((r0_mean_vec-2*pi))))
maxy = max(maxy,max(abs((r3_mean_vec-2*pi))))
maxy = max(maxy,max(abs((r4_mean_vec-2*pi))))
plot(x = N2_vec, y = abs((r0_mean_vec-2*pi)), col = 'blue', pch = 20,ylim = c(miny, maxy))

points(x = N2_vec, y = abs((r2_mean_vec-2*pi)), col = 'red', pch = 15)
points(x = N2_vec, y = abs((r3_mean_vec-2*pi)), col = 'purple', pch = 10)
points(x = N2_vec, y = abs((r4_mean_vec-2*pi)), col = 'orange', pch = 50)


# variance
plot(x = N2_vec, y = r0_sd_vec, col = 'blue', pch = 20)
plot(x = N2_vec, y = r2_sd_vec, col = 'red', pch = 15)
plot(x = N2_vec, y = r3_sd_vec, col = 'purple', pch = 10)
plot(x = N2_vec, y = r4_sd_vec, col = 'orange', pch = 50)

miny = min(min(r2_sd_vec), min(r0_sd_vec))
miny = min(miny,min(r3_sd_vec))
miny = min(miny,min(r4_sd_vec))
maxy = max(min(r2_sd_vec), max(r0_sd_vec))
maxy = max(maxy, max(r3_sd_vec))
maxy = max(maxy, max(r4_sd_vec))
plot(x = N2_vec, y = r0_sd_vec, col = 'blue', pch = 20,ylim = c(miny, maxy))
points(x = N2_vec, y = r2_sd_vec, col = 'red', pch = 15)
points(x = N2_vec, y = r3_sd_vec, col = 'purple', pch = 10)
points(x = N2_vec, y = r4_sd_vec, col = 'orange', pch = 50)

arrows(x0 = N2_vec, y0 = (r0_mean_vec-r0_sd_vec) , x1= N2_vec, y1 = (r0_mean_vec-r0_sd_vec), length=0.05, angle=90, code=3)


#######################
import_data <- rphi_mix(simulation_sample_size/2)
weights <- p1(import_data)/phi_mix(import_data)
import_sampling_half_log10_Z_hat <- log10(mean(weights))
import_sampling_half_log10_Z_hat

import_data <- rphi_mix(simulation_sample_size)
weights <- p1(import_data)/phi_mix(import_data)
import_sampling_all_log10_Z_hat <- log10(mean(weights))
import_sampling_all_log10_Z_hat

save(list=c("import_sampling_half_log10_Z_hat","import_sampling_all_log10_Z_hat","log10_bs_Z_hat","Z_hat_store","target_evaluations","proposal_evaluations","standard_normal_evaluations","K_prop","simulation_sample_size","sample_size","N2"),file=paste(home,save_name,".RData",sep=""))
