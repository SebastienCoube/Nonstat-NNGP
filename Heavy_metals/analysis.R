# noise, scale, range

X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL


mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_isotropic_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = X, # range for latent field of parameters, if NULL no latent field
  scale_X = X, # range for latent field of parameters, if NULL no latent field
  range_X = X, # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 1
)

remove(X)
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, burn_in = .2)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_nsr")


# noise, scale
remove(mcmc_nngp_list) ; gc()

X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
 
mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "exponential_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = X, # range for latent field of parameters, if NULL no latent field
  scale_X = X, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X)
for(i in seq(60))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_ns")


# stat
remove(mcmc_nngp_list) ; gc()

X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL

mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "exponential_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = NULL, # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X)
for(i in seq(50))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_stat")


 