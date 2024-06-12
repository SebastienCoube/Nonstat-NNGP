Sys.sleep(rpois(1, 20))

libraries_dir = "../R_packages"
library(GpGp, lib.loc = libraries_dir)
library(FNN, lib.loc = libraries_dir)
library(Bidart, lib.loc = libraries_dir)
library(abind, lib.loc = libraries_dir)
library(parallel, lib.loc = libraries_dir)
library(Matrix, lib.loc = libraries_dir)
library(Rcpp, lib.loc = libraries_dir)
library(RcppArmadillo, lib.loc = libraries_dir)
library(magrittr)


model_range = c("stat", "nonstat_circular", "nonstat_elliptic")
data_range  = c("stat", "nonstat_circular", "nonstat_elliptic")

seed = seq(1, 60)

inputs  = expand.grid(
  "model_range" = model_range,
  "data_range" = data_range ,
  "seed" = seed
)


for(tatato in seq(nrow(inputs)))
{
  i = seq(nrow(inputs)) [sapply(sapply(seq(nrow(inputs)), function(i)paste("res", i, "started.RDS", sep = "")), function(x)!x%in%list.files())][1]
  saveRDS("tatato", paste("res", i, "started.RDS", sep = ""))
  
  res = list()
  set.seed(inputs$seed[i])
  # locations and latent fields
  locs = 5 * matrix(runif(24000), ncol = 2)
  locs = locs[GpGp::order_maxmin(locs),]
  
  # getting PP basis
  n_PP = 49
  PP_range = .5
  data_PP = Bidart::get_PP(observed_locs = locs, matern_range = PP_range, lonlat = F, n_PP = n_PP, m = 15)
  res$inputs$data_PP = data_PP
  
  ###  Visualisation of the PP basis
  ###  comparison between PP and NNGP
  seed_vector =  rnorm(data_PP$n_PP + nrow(data_PP$unique_reordered_locs))
  par(mfrow = c(1,2))
  Bidart::plot_pointillist_painting(locs, Bidart::X_PP_mult_right(PP = data_PP, use_PP = T, Y = seed_vector[seq(data_PP$n_PP)]), cex = .3, main ="NNGP into PP")
  points(data_PP$knots, pch = 3, cex = .3)
  Bidart::plot_pointillist_painting(rbind(data_PP$knots, data_PP$unique_reordered_locs), as.vector(Matrix::solve(data_PP$sparse_chol, seed_vector)), cex = .3, main = "NNGP")
  par(mfrow = c(1,1))
  
  
  # model settings
  range_PP = F
  anisotropic  = F
  range_beta = matrix(0, n_PP+1, 3)
  range_beta_1_1 =  log(.1)
  range_beta[1, 1] =  range_beta_1_1
  range_var = 1
  if(inputs$data_range[i] == "nonstat_circular") 
  {
    range_beta[-1,1] = sqrt(range_var)*rnorm(n_PP)
  }
  if(inputs$data_range[i] == "nonstat_elliptic") 
  {
    range_beta[-1,] = sqrt(range_var)*rnorm(3*n_PP)
  }
  
  # computing sparse chol
  NNarray = GpGp::find_ordered_nn(locs, 5)
  sparse_chol = Bidart::compute_sparse_chol(
    range_beta = range_beta, 
    NNarray = NNarray, 
    locs = locs, 
    range_X = matrix(1, nrow(locs)), 
    PP = data_PP, 
    use_PP = T, 
    compute_derivative = T, 
    nu = .5, 
    anisotropic = T, 
    sphere = F, 
    num_threads = parallel::detectCores(), 
    locs_idx = NULL
  )
  
  # latent field
  scaled_latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)
  Bidart::plot_pointillist_painting(locs, scaled_latent_field, cex = .5)
  # doing observations, with duplicates
  n_obs = 20000
  observation_idx = c(seq(10000), sample(10000, n_obs - 10000 , replace =T))
  observed_locs = locs[observation_idx,]
  observed_field = as.vector(scaled_latent_field[observation_idx]) + rnorm(n_obs)
  Bidart::plot_pointillist_painting(observed_locs, observed_field, cex = .5)
  
  # storing inputs
  res$inputs = list()
  res$inputs$inputs = inputs[i,]
  res$inputs$sparse_chol         = sparse_chol[[1]]
  res$inputs$observed_field      = observed_field
  res$inputs$locs                = locs
  res$inputs$observation_idx     = observation_idx
  res$inputs$range_beta     = range_beta
  

  # initializing the model
  PP = Bidart::get_PP(
    observed_locs = observed_locs,
    matern_range = PP_range, 
    lonlat = F, n_PP = n_PP, m = 15)
    
  mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
    observed_locs = observed_locs, #spatial locations
    observed_field = c(observed_field), # Response variable
    X = NULL, # Covariates per observation
    m = 5, #number of Nearest Neighbors
    nu = .5, 
    anisotropic = inputs$model_range[i] == "nonstat_elliptic",
    sphere = F, PP = PP, 
    range_PP = !inputs$model_range[i] == "stat",
    n_chains = 2,  # number of MCMC chains
    seed = 1
  )
  
  # running the model
  print(inputs[i,])
  res$run = mcmc_nngp_list
  for(j in seq(15)){
    res$run = Bidart::mcmc_nngp_run_nonstationary_socket(
      res$run, n_cores = 2,
      debug_outfile = NULL, 
      num_threads_per_chain = parallel::detectCores()/2, 
      thinning = .1, 
      lib.loc = libraries_dir, 
      burn_in = .5
        )
  }
  
  # analysis of performance
  res$performance = list()
  true_field = as.vector(scaled_latent_field)
  
  # estimation
  res$estimation = Bidart::estimate_parameters(mcmc_nngp_list = res$run, get_samples = T, burn_in = .5)
  
  # field MSE at observed locations
  res$performance$smooth_field_mse = 
    mean((res$estimation$summaries$field_at_observed_locs[1,,1] - true_field[observation_idx])^2)
  #plot(res$estimation$summaries$field_at_observed_locs[1,,1] , true_field[observation_idx])
  
  # prediction of the latent field
  res$pred = Bidart::predict_latent_field(mcmc_nngp_list = res$run, predicted_locs = locs[10001:12000,], num_threads_per_chain = parallel::detectCores(), lib.loc = libraries_dir, parallel = T, burn_in = .5)
  res$performance$pred_field_mse = mean((res$pred$summaries$field[1,,1] - true_field[seq(10001, 12000)])^2)
  #plot(res$pred$summaries$field[1,,1], true_field[seq(10001, 12000)])
  
  # DIC
  res$performance$DIC = Bidart::DIC(res$run)
  
  
  # estimation of log-range field
  true_log_range = Bidart::X_PP_mult_right(
    X = matrix(1, nrow(observed_locs)), PP = data_PP, 
    use_PP = inputs$data_range[i]!="stat", 
    Y = range_beta
  )
  estimated_log_range = apply(res$estimation$samples$range_beta, 3, function(x)
    Bidart::X_PP_mult_right(
      X = matrix(1, nrow(observed_locs)), PP = PP, 
      use_PP = inputs$model_range[i]!="stat", 
      Y = x
    ), simplify = F
  ) 
  estimated_log_range = Reduce(function(x, y)abind::abind(x, y, along = 3), estimated_log_range) %>% 
    apply(MARGIN = c(1,2), FUN = mean, simplify = T) 
  ## Bidart::plot_pointillist_painting(locs[observation_idx,], true_log_range[observation_idx,1])
  ## Bidart::plot_pointillist_painting(observed_locs, estimated_log_range[,1])
  ## plot(estimated_log_range[,1],  true_log_range[observation_idx,1])
  ## plot(estimated_log_range[,2],  true_log_range[observation_idx,2])
  res$performance$log_range_det_mse = mean((unlist(estimated_log_range[,1,drop=F]) - true_log_range[observation_idx,1])^2)
  if(inputs$model_range[i]=="nonstat_elliptic"){
    res$performance$log_range_aniso_mse = (
      mean(((estimated_log_range[,2]) - true_log_range[observation_idx,2])^2) + 
        mean(((estimated_log_range[,3]) - true_log_range[observation_idx,3])^2) 
    )/2
  }
  if(inputs$model_range[i]=="nonstat_circular"){
    res$performance$log_range_aniso_mse = (
      mean((true_log_range[observation_idx,2])^2) + 
        mean((true_log_range[observation_idx,3])^2) 
    )/2
  }
  if(inputs$model_range[i]=="stat"){
    res$performance$log_range_aniso_mse = (
      mean((true_log_range[observation_idx,2])^2) + 
        mean((true_log_range[observation_idx,3])^2) 
    )/2
  }
  
  # credible intervals coverage
  res$coverage = list()
  res$coverage$range =  
    (c(range_beta_1_1, 0, 0 )>=c(res$estimation$summaries$range_beta[2,1,], rep(0, 3-length(res$estimation$summaries$range_beta[2,1,]))))&
    (c(range_beta_1_1, 0, 0 )<=c(res$estimation$summaries$range_beta[4,1,], rep(0, 3-length(res$estimation$summaries$range_beta[4,1,]))))
  res$coverage$scale =  
    (0>res$estimation$summaries$scale_beta[2,1,1])&
    (0<res$estimation$summaries$scale_beta[4,1,1])
  res$coverage$noise =  
    (0>res$estimation$summaries$noise_beta[2,1,1])&
    (0<res$estimation$summaries$noise_beta[4,1,1])
  
  true_range_var = rep(0, 6)
  if(inputs$data_range[i] == "nonstat_elliptic")true_range_var[c(1,4,6)] = range_var
  if(inputs$data_range[i] == "nonstat_circular")true_range_var[c(1)] = range_var
  
  if(inputs$model_range[i]=="nonstat_elliptic"){
    range_var_quantiles = 
      apply(res$estimation$samples$range_log_scale, 3, Bidart::expmat, simplify = F) %>% 
      lapply(function(x)x[lower.tri(x, T)]) %>% 
      do.call(what = rbind) %>% 
      apply(2, quantile, probs = c(.025, .975))
  }
  if(inputs$model_range[i]=="nonstat_circular"){
    range_var_quantiles = 
      apply(res$estimation$samples$range_log_scale, 3, Bidart::expmat, simplify = F) %>% 
      lapply(function(x)x[lower.tri(x, T)]) %>% 
      do.call(what = rbind) %>% 
      apply(2, quantile, probs = c(.025, .975))
    range_var_quantiles = cbind(range_var_quantiles, matrix(0, 2,5))
  }
  if(inputs$model_range[i]=="stat")range_var_quantiles = matrix(0, 2,6)
  
  res$coverage$range_var = (true_range_var>= range_var_quantiles[1,])&(true_range_var<= range_var_quantiles[2,])
  
  # saving the run (or not)
  if(inputs$seed[i]>3)res$run = NULL
  
  saveRDS(res, paste("res", i, "complete.RDS", sep = ""))
}

