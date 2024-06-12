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


model_range = c("stat", "nonstat")
model_scale = c("stat", "nonstat")
model_noise = c("stat", "nonstat")
data_range  = c("stat", "nonstat")
data_scale  = c("stat", "nonstat")
data_noise  = c("stat", "nonstat")

seed = seq(1, 60)

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
  started_complete = readRDS("started_complete.RDS")
  i = seq(nrow(inputs)) [-started_complete$started][1]
  started_complete$started = c(started_complete$started, i)
  saveRDS(started_complete, "started_complete.RDS")
  
  res = list()
  print(inputs[i,])
  # locations and latent fields
  locs = 5 * matrix(runif(24000), ncol = 2)
  locs = locs[GpGp::order_maxmin(locs),]
  
  # getting PP basis
  n_PP = 49
  PP_range = .5
  data_PP = Bidart::get_PP(observed_locs = locs, matern_range = PP_range, lonlat = F, n_PP = n_PP, m = 15)
  Bidart::compare_PP_NNGP(data_PP)
  res$inputs$data_PP = data_PP
  
  params_var = 1
  range_beta_1_1 = log(.1)
  range_beta = matrix(c(range_beta_1_1, (inputs$data_range[i] == "nonstat") * sqrt(params_var) *  rnorm(n_PP)))
  scale_beta = matrix(c(0,(inputs$data_scale[i] == "nonstat")   * sqrt(params_var) *  rnorm(n_PP)))
  noise_beta = matrix(c(0, (inputs$data_noise[i] == "nonstat")  * sqrt(params_var) *  rnorm(n_PP)))
  
  # observing observations, with duplicates
  n_obs = 20000
  observation_idx = c(sample(10000, 10000 , replace =F), sample(10000, n_obs - 10000 , replace =T))
  observed_locs = locs[observation_idx,]

  # computing sparse chol
  NNarray = GpGp::find_ordered_nn(locs, 10)
  sparse_chol = Bidart::compute_sparse_chol(
    anisotropic = F, 
    range_beta = range_beta, 
    locs = locs, 
    NNarray = NNarray, 
    range_X = matrix(rep(1, nrow(locs))), PP = data_PP, use_PP = T, nu = .5)
  
  true_field_log_range = log(Bidart::variance_field(beta = range_beta, PP = data_PP, use_PP = T, X = matrix(rep(1, nrow(locs)))))
    Bidart::plot_pointillist_painting(locs, true_field_log_range)
  true_field_log_var =   log(Bidart::variance_field(beta = scale_beta, PP = data_PP, use_PP = T, X = matrix(rep(1, nrow(locs)))))
    Bidart::plot_pointillist_painting(locs, true_field_log_var)
  true_noise_log_var =   log(Bidart::variance_field(beta = noise_beta, PP = data_PP, use_PP = T, X = matrix(rep(1, nrow(locs)))))
    Bidart::plot_pointillist_painting(locs, true_noise_log_var)
  
  # latent field
  scaled_latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)
    Bidart::plot_pointillist_painting(locs, scaled_latent_field, cex = .5)
  true_field = (exp(.5*true_field_log_var) * scaled_latent_field)
    Bidart::plot_pointillist_painting(locs, true_field, cex = .5)
  noise = (exp(.5*true_noise_log_var)[observation_idx]* rnorm(n_obs)) 
    Bidart::plot_pointillist_painting(observed_locs, noise, cex = .5)
    Bidart::plot_pointillist_painting(observed_locs, exp(.5*true_noise_log_var)[observation_idx], cex = .5)
  observed_field = true_field[observation_idx]  + noise
    Bidart::plot_pointillist_painting(observed_locs, observed_field)
  
  # storing inputs
  res$inputs = list()
  res$inputs$inputs = inputs[i,]
  res$inputs$scale_beta  = scale_beta
  res$inputs$noise_beta  = noise_beta
  res$inputs$range_beta  = range_beta
  res$inputs$sparse_chol         = sparse_chol[[1]]
  res$inputs$scaled_latent_field = scaled_latent_field
  res$inputs$observed_field      = observed_field
  res$inputs$true_field          = true_field
  res$inputs$noise               = noise
  res$inputs$locs                = locs
  res$inputs$observation_idx     = observation_idx
  
  # model settings
  noise_PP = F
  scale_PP = F
  range_PP = F
  if(inputs$model_noise[i] == "nonstat"){
    noise_PP =T
  }
  if(inputs$model_scale[i] == "nonstat"){
    scale_PP=T
  }
  if(inputs$model_range[i] == "nonstat") {
    range_PP = T
  }
  PP = Bidart::get_PP(observed_locs = observed_locs, matern_range = PP_range, lonlat = F, n_PP = n_PP, m = 15)
  
  # initializing the model
  res$run = Bidart::mcmc_nngp_initialize_nonstationary (
    observed_locs = observed_locs, #spatial locations
    observed_field = c(observed_field), # Response variable
    X = NULL, # Covariates per observation
    m = 5, #number of Nearest Neighbors
    noise_X = NULL, noise_PP = noise_PP,
    scale_X = NULL, scale_PP = scale_PP,
    range_X = NULL, range_PP = range_PP,
    PP = PP, 
    n_chains = 2,  # number of MCMC chains
    seed = 2, nu = .5
  )
  
  
  ##state = res$run$states$chain_2
  ##hierarchical_model = res$run$hierarchical_model
  ##data = res$run$data
  ##vecchia_approx = res$run$vecchia_approx
  
  # running the model
  for(j in seq(15)){
    print(i)
    res$run = Bidart::mcmc_nngp_run_nonstationary_socket(
      res$run, 
      n_cores = 2, 
      debug_outfile = NULL, 
      num_threads_per_chain = ceiling(parallel::detectCores()/2), 
      lib.loc = libraries_dir, 
      burn_in = .5
      )
  }
  
  # analysis of performance
  res$performance = list()
  
  # estimation
  res$estimation = Bidart::estimate_parameters(mcmc_nngp_list = res$run, get_samples = T, burn_in = .5)
  # field MSE at observed
  res$performance$smooth_field_mse = 
    mean((res$estimation$summaries$field_at_observed_locs[1,,1] - true_field[observation_idx])^2)
  
  
  # Estimating fields of parameters
  estimated_log_range = apply(
    res$estimation$samples$range_beta, 3, 
    function(x)
      Bidart::X_PP_mult_right(
        X = matrix(1, nrow(observed_locs)), PP = PP, 
        use_PP = inputs$model_range[i]!="stat", 
        Y = x
      )
  ) %>% apply(MARGIN = 1, FUN = mean)
  #Bidart::plot_pointillist_painting(observed_locs, estimated_log_range)
  #Bidart::plot_pointillist_painting(observed_locs, true_field_log_range[observation_idx])
  res$performance$log_range_mse = mean((estimated_log_range - true_field_log_range[observation_idx])^2)
  
  estimated_log_scale = apply(res$estimation$samples$scale_beta, 3, function(x)
    Bidart::X_PP_mult_right(
      X = matrix(1, nrow(observed_locs)), PP = PP, 
      use_PP = inputs$model_scale[i]!="stat", 
      Y = x
    )
  ) %>% apply(MARGIN = 1, FUN = mean)
  #Bidart::plot_pointillist_painting(observed_locs, estimated_log_scale)
  #Bidart::plot_pointillist_painting(observed_locs, true_field_log_var[observation_idx])
  res$performance$log_scale_mse = mean((estimated_log_scale - true_field_log_var[observation_idx])^2)
  
  estimated_log_noise = apply(res$estimation$samples$noise_beta, 3, function(x)
    Bidart::X_PP_mult_right(
      X = matrix(1, nrow(observed_locs)), PP = PP, 
      use_PP = inputs$model_noise[i]!="stat", 
      Y = x
    )
  ) %>% apply(MARGIN = 1, FUN = mean)
  res$performance$log_noise_mse = mean((estimated_log_noise - true_noise_log_var[observation_idx])^2)
  #Bidart::plot_pointillist_painting(observed_locs, estimated_log_noise)
  #Bidart::plot_pointillist_painting(observed_locs, true_noise_log_var[observation_idx])
  # prediction
  res$pred = Bidart::predict_latent_field(mcmc_nngp_list = res$run, predicted_locs = locs[10001:12000,], num_threads_per_chain = parallel::detectCores(), lib.loc = libraries_dir, parallel = T, burn_in = .5)
  res$performance$pred_field_mse = mean((res$pred$summaries$field[1,,1] - true_field[seq(10001, 12000)])^2)
  #plot(res$pred$summaries$field[1,,1], true_field[seq(10001, 12000)])
  # DIC
  res$performance$DIC = Bidart::DIC(res$run)
  
  
  # coverage
  res$coverage = list()
  res$coverage$range_beta = 
    (range_beta_1_1 >= res$estimation$summaries$range_beta[2,1,])&
    (range_beta_1_1 <= res$estimation$summaries$range_beta[4,1,])
  res$coverage$scale_beta = 
    (0 >= res$estimation$summaries$scale_beta[2,1,])&
    (0 <= res$estimation$summaries$scale_beta[4,1,])
  res$coverage$noise_beta = 
    (0 >= res$estimation$summaries$noise_beta[2,1,])&
    (0 <= res$estimation$summaries$noise_beta[4,1,])
  range_var_bounds = c(0,0)
  if(res$inputs$inputs$model_range== "nonstat")range_var_bounds = exp(res$estimation$summaries$range_log_scale[c(2,4),,])
  true_range_var =  params_var * (res$inputs$inputs$data_range == "nonstat")
  res$coverage$range_var = 
    (true_range_var >= range_var_bounds[1])&
    (true_range_var <= range_var_bounds[2])
  scale_var_bounds = c(0,0)
  if(res$inputs$inputs$model_scale== "nonstat")scale_var_bounds = exp(res$estimation$summaries$scale_log_scale[c(2,4),,])
  true_scale_var =  params_var * (res$inputs$inputs$data_scale == "nonstat")
  res$coverage$scale_var = 
    (true_scale_var >= scale_var_bounds[1])&
    (true_scale_var <= scale_var_bounds[2])
  noise_var_bounds = c(0,0)
  if(res$inputs$inputs$model_noise== "nonstat")noise_var_bounds = exp(res$estimation$summaries$noise_log_scale[c(2,4),,])
  true_noise_var =  params_var * (res$inputs$inputs$data_noise == "nonstat")
  res$coverage$noise_var = 
    (true_noise_var >= noise_var_bounds[1])&
    (true_noise_var <= noise_var_bounds[2])

  
  # saving the run (or not)
  if(inputs$seed[i]>3)res$run = NULL
  
  saveRDS(res, paste("res", i, "complete.RDS", sep = ""))
  started_complete = readRDS("started_complete.RDS")
  started_complete$complete = c(started_complete$complete, i)
  saveRDS(started_complete, "started_complete.RDS")
}

