
#####################
# Full data set run #
#####################

  ####################
  # KL decomposition #
  ####################

locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs
KL_decomposition = Bidart::get_KL_basis(locs, m = 10, n_PP = 500, n_KL = 25, covparms = c(1, 3, 1, 0))
saveRDS(KL_decomposition, "Heavy_metals/KL_decomposition.RDS")

  #######
  # run #
  #######

# stat
remove(mcmc_nngp_list) ; gc()

X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL

mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 10, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "exponential_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = NULL, # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X)
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, swap_range_scale = T)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_stat")

# noise + range nonstat
gc()

X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
KL_decomposition = readRDS("Heavy_metals/KL_decomposition.RDS")

mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 10, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_isotropic_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(KL_decomposition$basis[,seq(25)], X)), 
  scale_X = NULL,
  range_X = as.data.frame(cbind(KL_decomposition$basis[,seq(25)], X[, c("gcarb", "globedem", "twi")])), 
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X)
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_nr")





##############
# Validation #
##############

####################
# train-test split #
####################
if(length(grep("validation_train.RDS", list.files("Heavy_metals/")))==0)
{
  # loading data
  X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
  X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
  locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs
  field = readRDS("Heavy_metals/processed_data.RDS")$observed_field
  
  # shuffle locs and select 50000 train locs
  set.seed(1234)
  locs_ = unique(locs);locs_ = locs_[GpGp::order_maxmin(locs_),];locs_ = locs_[seq(50000),]
  train_locs = locs_  [na.omit(match(split(locs, row(locs)), split(locs_, row(locs_)))),]
  train_field = field [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_))))  ]
  train_X     = X     [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))), ]
  # get KL basis
  train_KL_decomposition      = Bidart::get_KL_basis(train_locs, lonlat = T, covparms = c(1, 3, 1, 0), m = 10)
  # save
  saveRDS(list(train_locs = train_locs, train_field = train_field, train_X = train_X, train_KL_decomposition = train_KL_decomposition), "Heavy_metals/validation_train.RDS")
  
  # get test data set
  test_locs = locs  [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),]  # get locs
  test_KL = Bidart::predict_KL_basis(test_locs, train_KL_decomposition) # predict KL basis
  # save
  saveRDS(list(
    test_locs = test_locs, 
    test_KL = test_KL, 
    test_X    = X     [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),],
    test_obs  = field [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))) ]
  ), "Heavy_metals/test_data_set.RDS")
  
}

##########################################################
# stationary with spatial basis and removed observations #
##########################################################
train_data_set = readRDS("Heavy_metals/validation_train.RDS")
mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs  = train_data_set$train_locs, #spatial locations
  observed_field = train_data_set$train_field, # Response variable
  X = train_data_set$train_X, # Covariates per observation
  m = 10, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "exponential_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = NULL, # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 1
)
remove(train_data_set)
for(i in seq(30))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, swap_range_scale = T)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_stat_validation.RDS")
 
#################################
# nonstationary noise and range #
#################################
train_data_set = readRDS("Heavy_metals/validation_train.RDS")
mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs  = train_data_set$train_locs, #spatial locations
  observed_field = train_data_set$train_field, # Response variable
  X = train_data_set$train_X, # Covariates per observation
  m = 10, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_isotropic_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(train_data_set$train_KL_decomposition$basis[,seq(25)], train_data_set$train_X)), # range for latent field of parameters, if NULL no latent field
  range_X = as.data.frame(cbind(train_data_set$train_KL_decomposition$basis[,seq(25)], train_data_set$train_X[, c("gcarb", "globedem", "twi")])), # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)
remove(train_data_set)
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_nr_basis_validation.RDS")


#######################
# nonstationary noise #
#######################
train_data_set = readRDS("Heavy_metals/validation_train.RDS")
mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs  = train_data_set$train_locs, #spatial locations
  observed_field = train_data_set$train_field, # Response variable
  X = train_data_set$train_X, # Covariates per observation
  m = 10, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "exponential_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(train_data_set$train_KL_decomposition$basis[,seq(25)], train_data_set$train_X)), # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)
remove(train_data_set)
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, swap_range_scale = T)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_n_validation.RDS")

