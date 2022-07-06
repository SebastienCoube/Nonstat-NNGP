

#install.packages("/home/user/s/scoube/Bidart_1.0.tar.gz", lib = "/home/user/s/scoube/R_packages/")
#install.packages("abind",  "https://cloud.r-project.org", lib = "/home/user/s/scoube/R_packages/")
#install.packages("expm",  "https://cloud.r-project.org", lib = "/home/user/s/scoube/R_packages/")
#install.packages("GpGp",  "https://cloud.r-project.org", lib = "/home/user/s/scoube/R_packages/")




library(GpGp, lib.loc = "/home/user/s/scoube/R_packages/")
library(FNN, lib.loc = "/home/user/s/scoube/R_packages/")
library(Bidart, lib.loc = "/home/user/s/scoube/R_packages/")
library(abind, lib.loc = "/home/user/s/scoube/R_packages/")
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
  
  set.seed(inputs$seed[i])
  # locations and latent fields
  locs = 5 * matrix(runif(24000), ncol = 2)
  locs = locs[GpGp::order_maxmin(locs),]
  latent_field_scale = (inputs$data_scale[i] == "nonstat") * GpGp::fast_Gp_sim(c(.5, .5, 1, 0), locs = locs)
  latent_field_noise = (inputs$data_noise[i] == "nonstat") * GpGp::fast_Gp_sim(c(.5, .5, 1, 0), locs = locs)
  latent_field_range = (inputs$data_range[i] == "nonstat") * GpGp::fast_Gp_sim(c(.5, .5, 1, 0), locs = locs)
  range_beta = log(.1)
  scale_beta = log(1)
  noise_beta = log(.5)
  
  # observing observations, with duplicates
  n_obs = 20000
  observation_idx = c(sample(10000, 10000 , replace =F), sample(10000, n_obs - 10000 , replace =T))
  observed_locs = locs[observation_idx,]
  
  # computing sparse chol
  NNarray = GpGp::find_ordered_nn(locs, 5)
  sparse_chol = compute_sparse_chol(covfun_name = "nonstationary_exponential_isotropic", range_beta = range_beta, locs = locs, NNarray = NNarray, range_field =  latent_field_range, range_X = matrix(rep(1, nrow(locs))))
  
  # latent field
  scaled_latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)
  observed_field = as.vector(exp(.5 * scale_beta + .5*latent_field_scale[observation_idx])* scaled_latent_field[observation_idx] +exp(.5 * noise_beta + .5*latent_field_noise[observation_idx])*rnorm(n_obs))
  
  # storing inputs
  res = list()
  res$inputs = list()
  res$inputs$inputs = inputs[i,]
  res$inputs$latent_field_scale  = latent_field_scale
  res$inputs$latent_field_noise  = latent_field_noise
  res$inputs$latent_field_range  = latent_field_range
  res$inputs$sparse_chol         = sparse_chol[[1]]
  res$inputs$scaled_latent_field = scaled_latent_field
  res$inputs$observed_field      = observed_field
  res$inputs$locs                = locs
  res$inputs$observation_idx     = observation_idx
  
  # model settings
  covfun = "exponential_isotropic"
  noise_range = NULL
  scale_range = NULL
  range_range = NULL
  if(inputs$model_noise[i] == "nonstat") 
  {
    noise_range = .5
  }
  if(inputs$model_scale[i] == "nonstat")
  {
    scale_range = .5
  }
  if(inputs$model_range[i] == "nonstat") 
  {
    range_range = .5
    covfun = "nonstationary_exponential_isotropic"
  }
  
  # initializing the model
  mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
    observed_locs = observed_locs, #spatial locations
    observed_field = c(observed_field), # Response variable
    X = NULL, # Covariates per observation
    m = 5, #number of Nearest Neighbors
    reordering = c("maxmin"), #Reordering
    covfun = covfun, response_model = "Gaussian", # covariance model and response model
    noise_X = NULL, noise_range = noise_range, # range for latent field of parameters, if NULL no latent field
    scale_X = NULL, scale_range = scale_range, # range for latent field of parameters, if NULL no latent field
    range_X = NULL, range_range = range_range, # range for latent field of parameters, if NULL no latent field
    log_NNGP_matern_covfun = "matern_isotropic", # covariance function for the hyperpriors
    log_NNGP_nu = 1, # covariance function for the hyperpriors
    n_chains = 3,  # number of MCMC chains
    seed = 2
  )
  
  # running the model
  print(inputs[i,])
  res$run = mcmc_nngp_list
  for(j in seq(60))
  {
    res$run = mcmc_nngp_run_nonstationary(res$run, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL)
  }
  
  # analysis of performance
  res$performance = list()
  param_field_mse = function(estimated_field, true_field, estimated_intercept, true_intercept)
  {
    if(!is.null(estimated_field) & !is.null(true_field))return(mean(((estimated_field + estimated_intercept) -  (true_field + true_intercept))^2))
    if(is.null(estimated_field) & !is.null(true_field))return(mean(((estimated_intercept) -  (true_field + true_intercept))^2))
    if(!is.null(estimated_field) & is.null(true_field))return(mean(((estimated_field + estimated_intercept) -  (true_intercept))^2))
    if(is.null(estimated_field) & is.null(true_field))return(mean(((estimated_intercept) -  (true_intercept))^2))
  }
  true_field = as.vector(exp(.5*latent_field_scale)* scaled_latent_field)
  # estimation
  res$estimation = estimate_parameters(mcmc_nngp_list = res$run)
  res$performance$smooth_field_mse = param_field_mse(estimated_field = res$estimation$summaries$field[1,,], true_field = true_field[observation_idx][res$run$vecchia_approx$hctam_scol_1], estimated_intercept = res$estimation$summaries$beta[1,,], true_intercept = 0)
  # prediction
  res$pred = predict_latent_field(mcmc_nngp_list = res$run, predicted_locs = locs[10001:12000,], n_cores = 3)
  res$performance$pred_field_mse = mean(((res$pred$summaries$field[1,] + res$estimation$summaries$beta[1,,]) -  true_field[10001:12000])^2)
  # DIC
  res$performance$DIC = DIC(res$run)
  if(inputs$seed[i]>2)res$run = NULL
  
  saveRDS(res, paste("res", i, "complete.RDS", sep = ""))
}

