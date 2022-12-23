install.packages("/home/user/s/scoube/Bidart_1.0.tar.gz", lib = "/home/user/s/scoube/R_packages/")


library(GpGp, lib.loc = "/home/user/s/scoube/R_packages/")
library(FNN, lib.loc = "/home/user/s/scoube/R_packages/")
library(Bidart, lib.loc = "/home/user/s/scoube/R_packages/")
library(abind, lib.loc = "/home/user/s/scoube/R_packages/")
library(irlba, lib.loc = "/home/user/s/scoube/R_packages/")
library(parallel)
library(Matrix)


model_range = c("stat", "nonstat_circular", "nonstat_elliptic")
data_range  = c("stat", "nonstat_circular", "nonstat_elliptic")

seed = seq(1, 30)

inputs  = expand.grid(
  "model_range" = model_range,
  "data_range" = data_range ,
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
  if(inputs$data_range[i]=="stat")             range_beta  = cbind(c(log(.1), rep(0, 25)), c(0, rep(0, 25)), c(0, rep(0,25))) 
  if(inputs$data_range[i]=="nonstat_circular") range_beta  = cbind(c(log(.1), rnorm(25)),  c(0, rep(0, 25)), c(0, rep(0,25)))
  if(inputs$data_range[i]=="nonstat_elliptic") range_beta  = cbind(c(log(.1), rnorm(25)),  c(0, rnorm(25)),  c(0, rnorm(25)))
  
  n_KL = 25
  KL_covfun = "matern_isotropic"
  KL_covparms = c(1, .5, 1.5, 0.001)
  KL = get_KL_basis(locs = locs, covparms = KL_covparms, covfun_name = KL_covfun, n_KL = n_KL)
  
  # observing observations, with duplicates
  n_obs = 20000
  observation_idx = c(sample(10000, 10000 , replace =F), sample(10000, n_obs - 10000 , replace =T))
  observed_locs = locs[observation_idx,]
  
  # computing sparse chol
  NNarray = GpGp::find_ordered_nn(locs, 5)
  sparse_chol = Bidart::compute_sparse_chol(covfun_name = "nonstationary_exponential_anisotropic", range_beta = range_beta, locs = locs, NNarray = NNarray, KL = KL, use_KL = T, range_X = matrix(rep(1, nrow(locs))))
  
  # latent field
  scaled_latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)
  #plot_pointillist_painting(locs, scaled_latent_field, cex = .5)
  observed_field = as.vector(scaled_latent_field[observation_idx] * sqrt(10)) + rnorm(n_obs) * sqrt(10)
  
  
  # storing inputs
  res = list()
  res$inputs = list()
  res$inputs$inputs = inputs[i,]
  res$inputs$sparse_chol         = sparse_chol[[1]]
  res$inputs$observed_field      = observed_field
  res$inputs$locs                = locs
  res$inputs$observation_idx     = observation_idx
  
  # model settings
  covfun = "nonstationary_exponential_isotropic"
  range_KL = F
  if(inputs$model_range[i] == "nonstat_circular") 
  {
    range_KL = T
  }
  if(inputs$model_range[i] == "nonstat_elliptic") 
  {
    range_KL = T
    covfun = "nonstationary_exponential_anisotropic"
  }
  KL_ = get_KL_basis(locs = observed_locs, covparms = KL_covparms, covfun_name = KL_covfun, n_KL = n_KL)
  
  # initializing the model
  mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
    observed_locs = observed_locs, #spatial locations
    observed_field = c(observed_field), # Response variable
    X = NULL, # Covariates per observation
    m = 5, #number of Nearest Neighbors
    reordering = c("maxmin"), #Reordering
    covfun = covfun, 
    range_X = NULL, 
    range_KL = range_KL, KL = KL_, range_log_scale_prior = c(-8, 3),
    n_chains = 2,  # number of MCMC chains
    seed = 1
  )
  
  # running the model
  print(inputs[i,])
  res$run = mcmc_nngp_list
  for(j in seq(60))
  {
    res$run = mcmc_nngp_run_nonstationary(res$run, n_cores = 3, debug_outfile = NULL, big_range = T)
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
  true_field = as.vector(scaled_latent_field)
  # estimation
  res$estimation = estimate_parameters(mcmc_nngp_list = res$run)
  res$performance$smooth_field_mse = param_field_mse(estimated_field = res$estimation$summaries$field[1,,], true_field = true_field[observation_idx][res$run$vecchia_approx$hctam_scol_1], estimated_intercept = res$estimation$summaries$beta[1,,], true_intercept = 0)
  # prediction
  res$pred = predict_latent_field(mcmc_nngp_list = res$run, predicted_locs = locs[10001:12000,], n_cores = 3)
  res$performance$pred_field_mse = mean(((res$pred$summaries$field[1,] + res$estimation$summaries$beta[1,,]) -  true_field[10001:12000])^2)
  # DIC
  res$performance$DIC = DIC(res$run)
  if(inputs$seed[i]>3)res$run = NULL
  
  saveRDS(res, paste("res", i, "complete.RDS", sep = ""))
}

