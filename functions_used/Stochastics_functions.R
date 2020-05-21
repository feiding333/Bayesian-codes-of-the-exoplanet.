## get sample from GWL last stage, see details from WarpU Addition in WarpU_2019_version
GetSamples_from_LastStage = function(sample_size, estimate_result_GWL, david_data = 0){
  if(david_data == 0){
    sample_get = c()
    for (i in 1:sample_size){
      print(i)
      ## sample a subregion index k 
      prob_exclude_0 = estimate_result_GWL$g[2:length(estimate_result_GWL$g)]
      prob_exclude_0 = prob_exclude_0/sum(prob_exclude_0)
      one_subregion_ind = sample(x = (2:length(estimate_result_GWL$g)),size = 1,prob = prob_exclude_0)
      G_k_sample = estimate_result_GWL$x_final_stage[estimate_result_GWL$x_final_index == one_subregion_ind,]
      select_ind = sample(x = (1:(dim(G_k_sample)[1])),size = 1)
      sample_get = rbind(sample_get, G_k_sample[select_ind, ])
    }
    return(sample_get)
    
  }else{
    sample_get = c()
    s = 0
    while(s < sample_size){
      print(i)
      ## sample a subregion index k 
      g_estimate = exp(log_g_estimates)
      prob_exclude_0 =g_estimate[2:length(g_estimate)]
      prob_exclude_0 = prob_exclude_0/sum(prob_exclude_0)
      one_subregion_ind = sample(x = (2:length(g_estimate)),size = 1,prob = prob_exclude_0)
      G_k_sample = x_subregion_store[[one_subregion_ind]]
      print('one_subregion_ind')
      print(one_subregion_ind)
      print('dim(G_k_sample)[1]')
      print(dim(G_k_sample)[1])
      if(!is.null(dim(G_k_sample))){
        select_ind = sample(x = (2:(dim(G_k_sample)[1])),size = 1)
        sample_get = rbind(sample_get, G_k_sample[select_ind, ])
        s = s+1
      }
    }
    return(sample_get)
  }

}

## set the q_density for GWL example to apply stochastics Bridge estimation.
q_density_est = function(x_vec, log = F){
  result1 = rep(0,dim(x_vec)[1])
  for (i in (1:dim(x_vec)[1])) {
    x_vec1 = x_vec[i, ]
    tmpphix = (1/3)*exp((-1/2)*sum((x_vec1+5)^2)) + (2/3)*exp((-1/2)*sum((x_vec1-5)^2))
    result1[i] = tmpphix
  }
  if(log){
    return(log(result1))
  }
  if(!log){
    return(result1)
  }
}
