

#install.packages("/home/user/s/scoube/Bidart_1.0.tar.gz", lib = "/home/user/s/scoube/R_packages/")
#install.packages("abind",  "https://cloud.r-project.org", lib = "/home/user/s/scoube/R_packages/")
#install.packages("expm",  "https://cloud.r-project.org", lib = "/home/user/s/scoube/R_packages/")
#install.packages("GpGp",  "https://cloud.r-project.org", lib = "/home/user/s/scoube/R_packages/")
#install.packages("irlba",  "https://cloud.r-project.org", lib = "/home/user/s/scoube/R_packages/")




library(GpGp, lib.loc = "/home/user/s/scoube/R_packages/")
library(FNN, lib.loc = "/home/user/s/scoube/R_packages/")
library(Bidart, lib.loc = "/home/user/s/scoube/R_packages/")
library(abind, lib.loc = "/home/user/s/scoube/R_packages/")
library(irlba, lib.loc = "/home/user/s/scoube/R_packages/")
library(parallel)
library(Matrix)


model_range = c("stat", "nonstat")
model_scale = c("stat", "nonstat")
model_noise = c("stat", "nonstat")
data_range  = c("stat", "nonstat")
data_scale  = c("stat", "nonstat")
data_noise  = c("stat", "nonstat")

seed = seq(1, 30)

inputs  = expand.grid(
  "model_range" = model_range,
  "model_scale" = model_scale,
  "model_noise" = model_noise,
  "data_range" = data_range ,
  "data_scale" = data_scale ,
  "data_noise" = data_noise ,
  "seed" = seed
)


for(tatato in seq(nrow(inputs)))
{
  i = seq(nrow(inputs)) [sapply(sapply(seq(nrow(inputs)), function(i)paste("res", i, "started.RDS", sep = "")), function(x)!x%in%list.files())][1]
  saveRDS("tatato", paste("res", i, "started.RDS", sep = ""))
  print(inputs[i,])
  
  set.seed(inputs$seed[i])
  # locations and latent fields
  locs = 5 * matrix(runif(24000), ncol = 2)
  locs = locs[GpGp::order_maxmin(locs),]
  
  n_KL = 25
  KL_covfun = "matern_isotropic"
  KL_covparms = c(1, .5, 1.5, 0.001)
  KL = get_KL_basis(locs = locs, covparms = KL_covparms, covfun_name = KL_covfun, n_KL = n_KL)
  #plot(KL$KL_decomposition$d^2, ylim= c(0,KL$KL_decomposition$d[1]^2));abline(h=0)
  #Bidart::plot_pointillist_painting(locs, KL$basis[,1])
  
  
  range_beta = matrix(c(log(.1), (inputs$data_range[i] == "nonstat")        * 1 *  rnorm(n_KL)))
  #Bidart::plot_pointillist_painting(locs, cbind(1, KL$basis)%*%range_beta)
  scale_beta = matrix(c(log(10 ),(inputs$data_scale[i] == "nonstat")        * 1 *  rnorm(n_KL)))
  #Bidart::plot_pointillist_painting(locs, cbind(1, KL$basis)%*%scale_beta)
  noise_beta = matrix(c(log(10), (inputs$data_noise[i] == "nonstat")        * 1 *  rnorm(n_KL)))
  #Bidart::plot_pointillist_painting(locs, cbind(1, KL$basis)%*%noise_beta)
  
  # observing observations, with duplicates
  n_obs = 20000
  observation_idx = c(sample(10000, 10000 , replace =F), sample(10000, n_obs - 10000 , replace =T))
  observed_locs = locs[observation_idx,]

  # computing sparse chol
  field_covfun = "nonstationary_exponential_isotropic"
  NNarray = GpGp::find_ordered_nn(locs, 10)
  sparse_chol = compute_sparse_chol(covfun_name = field_covfun, range_beta = range_beta, locs = locs, NNarray = NNarray, range_X = matrix(rep(1, nrow(locs))), KL = KL, use_KL = T, nu = 1.5)
  
  # latent field
  scaled_latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)
  #Bidart::plot_pointillist_painting(locs, scaled_latent_field)
  true_field = (sqrt(variance_field(beta = scale_beta, KL = KL, use_KL = T, X = matrix(rep(1, nrow(locs))))) * scaled_latent_field)
  #Bidart::plot_pointillist_painting(locs, true_field)
  noise = (sqrt(variance_field(beta = noise_beta, KL = KL, use_KL = T, X = matrix(rep(1, nrow(locs))))))[observation_idx]* rnorm(n_obs) 
  #Bidart::plot_pointillist_painting(observed_locs, noise)
  observed_field = true_field[observation_idx]  + noise
  #Bidart::plot_pointillist_painting(observed_locs, observed_field)
  
  # storing inputs
  res = list()
  res$inputs = list()
  res$inputs$inputs = inputs[i,]
  res$inputs$scale_beta  = scale_beta
  res$inputs$noise_beta  = noise_beta
  res$inputs$range_beta  = range_beta
  res$KL = KL
  res$inputs$sparse_chol         = sparse_chol[[1]]
  res$inputs$scaled_latent_field = scaled_latent_field
  res$inputs$observed_field      = observed_field
  res$inputs$true_field          = true_field
  res$inputs$noise               = noise
  res$inputs$locs                = locs
  res$inputs$observation_idx     = observation_idx
  
  # model settings
  noise_KL = F
  scale_KL = F
  range_KL = F
  if(inputs$model_noise[i] == "nonstat") 
  {
    noise_KL =T
  }
  if(inputs$model_scale[i] == "nonstat")
  {
    scale_KL=T
  }
  if(inputs$model_range[i] == "nonstat") 
  {
    range_KL = T
  }
  KL_ = get_KL_basis(locs = observed_locs, covparms = KL_covparms, covfun_name = KL_covfun, n_KL = n_KL)
  
  
  # initializing the model
  res$run = mcmc_nngp_initialize_nonstationary (
    observed_locs = observed_locs, #spatial locations
    observed_field = c(observed_field), # Response variable
    X = NULL, # Covariates per observation
    m = 10, #number of Nearest Neighbors
    reordering = c("maxmin"), #Reordering
    covfun = field_covfun, 
    noise_X = NULL, noise_KL = noise_KL, noise_log_scale_prior = c(-8, 3),
    scale_X = NULL, scale_KL = scale_KL, scale_log_scale_prior = c(-8, 3),
    range_X = NULL, range_KL = range_KL, range_log_scale_prior = c(-8, 3),
    
    KL = KL_, 
    n_chains = 2,  # number of MCMC chains
    seed = 2
  )
  
  
  ##state = res$run$states$chain_2
  ##hierarchical_model = res$run$hierarchical_model
  ##data = res$run$data
  ##vecchia_approx = res$run$vecchia_approx
  
  # running the model
  for(j in seq(20))
  {
    print(i)
    res$run = mcmc_nngp_run_nonstationary(res$run, n_cores = 2, debug_outfile = NULL, big_range = T)
  }
  
  # analysis of performance
  res$performance = list()
  param_field_mse = function(estimated_field, true_field, estimated_intercept, true_intercept) mean(((estimated_field + estimated_intercept) -  (true_field + true_intercept))^2)
  #Bidart::plot_pointillist_painting(observed_locs, observed_field)
  # estimation
  res$estimation = estimate_parameters(mcmc_nngp_list = res$run)
  res$performance$smooth_field_mse = param_field_mse(estimated_field = res$estimation$summaries$field[1,,], true_field = true_field[observation_idx][res$run$vecchia_approx$hctam_scol_1], estimated_intercept = res$estimation$summaries$beta[1,,], true_intercept = 0)
  # prediction
  res$pred = predict_latent_field(mcmc_nngp_list = res$run, predicted_locs = locs[10001:12000,], n_cores = 3)
  res$performance$pred_field_mse = mean(((res$pred$summaries$field[1,] + res$estimation$summaries$beta[1,,]) -  true_field[10001:12000])^2)
  # DIC
  res$performance$DIC = DIC(res$run)
  if(inputs$seed[i]>2)res$run = NULL;res$pred = NULL
  
  saveRDS(res, paste("res", i, "complete.RDS", sep = ""))
}

