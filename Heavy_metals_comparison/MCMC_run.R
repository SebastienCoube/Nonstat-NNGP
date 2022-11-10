
#####################
# Full data set run #
#####################

  ####################
  # KL decomposition #
  ####################

#locs = readRDS("Heavy_metals_comparison/processed_data.RDS")$observed_locs
#KL_decomposition = Bidart::get_KL_basis(locs, m = 5, n_PP = 500, n_KL = 100, covparms = c(1, 5, 2.5, .0005))
#plot(KL_decomposition$KL_decomposition$d^2)
#plot(c(0, cumsum(KL_decomposition$KL_decomposition$d^2)))
#KL_decomposition = Bidart::get_KL_basis(locs, m = 5, n_PP = 500, n_KL = 15, covparms =  c(1, 5, 2.5, .0005))
#saveRDS(KL_decomposition, "Heavy_metals_comparison/KL_decomposition.RDS")

  #######
  # run #
  #######

# stat
remove(mcmc_nngp_list) ; gc()

X = readRDS("Heavy_metals_comparison/processed_data.RDS")$X_locs


mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals_comparison/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals_comparison/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 10, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "matern_sphere", nu = 1.5,  # covariance model
  noise_X = NULL, # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X)
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL)
}
saveRDS(mcmc_nngp_list, "Heavy_metals_comparison/run_stat")

# noise + range + scale nonstat
gc()

X = readRDS("Heavy_metals_comparison/processed_data.RDS")$X_locs
KL_decomposition = readRDS("Heavy_metals_comparison/KL_decomposition.RDS")
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL


mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals_comparison/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals_comparison/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 10, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_matern_isotropic_sphere", nu = 1.5, # covariance model
  noise_X = X,                                   noise_KL = T,  
  scale_X = X[, c("gcarb", "globedem", "twi")],  scale_KL = T, 
  range_X = X[, c("gcarb", "globedem", "twi")],  range_KL = T,  
  KL = KL_decomposition,
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X)
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL)
}
saveRDS(mcmc_nngp_list, "Heavy_metals_comparison/run_nsr")





##############
# Validation #
##############

####################
# train-test split #
####################
if(length(grep("validation_train.RDS", list.files("Heavy_metals_comparison/")))==0)
{
  # loading data
  X = readRDS("Heavy_metals_comparison/processed_data.RDS")$X_locs
  X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
  locs = readRDS("Heavy_metals_comparison/processed_data.RDS")$observed_locs
  field = readRDS("Heavy_metals_comparison/processed_data.RDS")$observed_field
  
  # shuffle locs and select 50000 train locs
  set.seed(1234)
  locs_ = unique(locs);locs_ = locs_[GpGp::order_maxmin(locs_),];locs_ = locs_[seq(50000),]
  train_locs = locs_  [na.omit(match(split(locs, row(locs)), split(locs_, row(locs_)))),]
  train_field = field [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_))))  ]
  train_X     = X     [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))), ]
  # get KL basis
  train_KL_decomposition = Bidart::get_KL_basis(train_locs, m = 10, n_PP = 500, n_KL = 20, covparms =  c(1, 4, 1.5, .001))
  Bidart::plot_pointillist_painting(train_locs, train_KL_decomposition$basis[,1])
  Bidart::plot_pointillist_painting(train_locs, train_KL_decomposition$basis[,2])
  Bidart::plot_pointillist_painting(train_locs, train_KL_decomposition$basis[,10])
  Bidart::plot_pointillist_painting(train_locs, train_KL_decomposition$basis[,20])
  plot(train_KL_decomposition$KL_decomposition$d^2)
  abline(h=0)
  # save
  saveRDS(list(train_locs = train_locs, train_field = train_field, train_X = train_X, train_KL_decomposition = train_KL_decomposition), "Heavy_metals_comparison/validation_train.RDS")
  # get test data set
  test_locs = locs  [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),]  # get locs
  test_KL = Bidart::predict_KL_basis(test_locs, train_KL_decomposition) # predict KL basis
  # save
  saveRDS(list(
    test_locs = test_locs, 
    test_KL = test_KL, 
    test_X    = X     [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),],
    test_obs  = field [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))) ]
  ), "Heavy_metals_comparison/test_data_set.RDS")
  
}


############################################
# nonstationary noise and range matern 15 #
###########################################
train_data_set = readRDS("Heavy_metals_comparison/validation_train.RDS")
TF = c(T, F)
for(nonstat_noise in TF)
{
for(nonstat_range in TF)
{
for(nonstat_scale in TF)
{
  if(length(grep(paste("run_noise", nonstat_noise, "scale", nonstat_scale, "range", nonstat_range, ".RDS", sep = "_"), list.files("Heavy_metals_comparison/")))==0)
  {
    scale_X = NULL
    noise_X = NULL
    range_X = NULL 
    if(nonstat_noise)noise_X = as.data.frame(train_data_set$train_X)                             
    if(nonstat_scale)scale_X = as.data.frame(train_data_set$train_X[, c("gcarb", "globedem", "twi")])                             
    if(nonstat_range)range_X = as.data.frame(train_data_set$train_X[, c("gcarb", "globedem", "twi")])                             
    covfun_name = "matern_sphere"
    if(nonstat_range)covfun_name = "nonstationary_matern_isotropic_sphere"
    mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
      observed_locs  = train_data_set$train_locs, #spatial locations
      observed_field = train_data_set$train_field, # Response variable
      X = train_data_set$train_X, # Covariates per observation
      m = 10, #number of Nearest Neighbors
      reordering = c("maxmin"), #Reordering
      covfun = covfun_name, 
      nu = 1.5, # Matérn smoothness of the covariance model
      scale_X = scale_X, scale_KL = nonstat_scale,
      noise_X = noise_X, noise_KL = nonstat_noise,  
      range_X = range_X, range_KL = nonstat_range, 
      KL = train_data_set$train_KL_decomposition, 
      n_chains = 2,  # number of MCMC chains
      seed = 1
    )
    for(i in seq(30))
    { 
      mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, debug_outfile = NULL)
    }
    #mcmc_nngp_update_Gaussian(mcmc_nngp_list$data, mcmc_nngp_list$hierarchical_model, mcmc_nngp_list$vecchia_approx, mcmc_nngp_list$states$chain_1, 100)
    saveRDS(mcmc_nngp_list, paste("Heavy_metals_comparison/run_noise", nonstat_noise, "scale", nonstat_scale, "range", nonstat_range, ".RDS", sep = "_"))
  }
}
}
}


