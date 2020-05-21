rm(list=ls())
jid=commandArgs(trailingOnly=T)[1]
jid=as.numeric(jid)

dataset_num <- 1
N_vec <- c(1,2,5,10,20)
iters_vec <- c(1000,3000,5000,10000)
grid_vec <- c(1001,2001,3001)

jid1 <- 1 + floor((jid-1)/(length(iters_vec)*length(grid_vec)))
jid2 <- 1 + floor((jid-1)/length(iters_vec)) - (jid1-1)*length(grid_vec)
jid3 <- jid - ((jid1-1)*length(grid_vec)+(jid2-1))*length(iters_vec)
  
K_prop = N_vec[jid1]# number of normal components in the mixture distribution
K_prop = 10
grid_num <- grid_vec[jid2]-1
grid_num  = 1000
iters <- iters_vec[jid3]
iters = 10000

num_planets <- "one"
save_name <- paste(num_planets,"_planet_full_prior_dataset_num",dataset_num,"_",grid_num,"grid_",iters,"iters_Kprop",K_prop,sep="")

home <- "E:/Dropbox/Exoplanet examples/EPRV code/"
home <- "E:/Dropbox/Warp Transformation/dave code/examples/one_planet/"
home <- "/n/home06/david_jones/exoplanets/"
data_home <- paste(home,"data/",sep="")
run_home <- paste(home,"runs/",sep="")
code_home <- paste(home,"code/",sep="")
setwd(run_home)
#load(paste(num_planets,"_planets_full_prior_fit_only_3chains_",grid_num,"grid_",iters,"iters.RData",sep=""))
load('../Stan HMC output/one_planets_full_prior_fit_only_3chains_1000grid_10000iters.RData')

fit_mat <- as.matrix(fit)

# Get posterior (including data)
#load(paste("dataset",dataset_num,"_one_planet_posterior.RData",sep=""))
#load("../set posterior function/dataset1_no_planet_posterior.RData")
load("../set posterior function/dataset1_one_planet_posterior.RData")


Pmin <- rep(1.25,3)
Pmax <- rep(10^4,3)
paras <- matrix(NA,2,7)
paras[1,] <- c(C=-1,tau=42.4,K=0.1,e=0.1,w=2.0,M0=3.0,sigmaJ=2)
paras[2,] <- c(C=-1,tau=42.4,K=0.1,e=0.1,w=2.0,M0=3.0,sigmaJ=3)
p1(paras)

# Libraries
library(mvtnorm)
library(mnormt)
library(numDeriv)
library(stats4)
library(MASS)
library(mvtnorm)

# Get Warp-U code
setwd(code_home)
source("EMhd_diagonal.R")
source("Bridge_Sampling_funs_eff.R")

# Evaluation costs:
# 1) Posterior draws - this we are taken as "given" i.e. the setting we consider is where posterior draws 
#                      have alrady been collected
# 2) EM - Normals with diagonal covariance matrices 
# 3) Warp-U: relatively few proposal evaluations and some target evaluations

# Split draws
simulation_sample_size = dim(fit_mat)[1] # Loaded stan HMC samples
n = simulation_sample_size
first_half = sample(1:n,size=n/2,replace=FALSE)# 1:(n/2)
second_half = setdiff(1:n,first_half) #(n/2+1):n
data = fit_mat[,1:dim(fit_mat)[2]-1]
for (i in c(3,7)){
  data[,i] <- exp(data[,i])-1
}
i <- 2
data[,i] <- exp(data[,i])
colnames(data) <- c("C","P","K","e","w","M0","sigma")
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
log_bs_Z_hat <- Opt_BS_hd(x1=data_tp1,x2=data_normal, q1=tp1,q2=p2_warp0,r0=10)
log10_bs_Z_hat <- log_bs_Z_hat/log(10)
Z_hat_store <- exp(log_bs_Z_hat)
log10_bs_Z_hat

import_data <- rphi_mix(simulation_sample_size/2)
weights <- p1(import_data)/phi_mix(import_data)
import_sampling_half_log10_Z_hat <- log10(mean(weights))
import_sampling_half_log10_Z_hat

import_data <- rphi_mix(simulation_sample_size)
weights <- p1(import_data)/phi_mix(import_data)
import_sampling_all_log10_Z_hat <- log10(mean(weights))
import_sampling_all_log10_Z_hat

save(list=c("import_sampling_half_log10_Z_hat","import_sampling_all_log10_Z_hat","log10_bs_Z_hat","Z_hat_store","target_evaluations","proposal_evaluations","standard_normal_evaluations","K_prop","simulation_sample_size","sample_size","N2"),file=paste(home,save_name,".RData",sep=""))
