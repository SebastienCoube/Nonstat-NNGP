
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
  PP      = Bidart::get_PP(train_locs, lonlat = T, matern_range = 10)
  # save
  saveRDS(list(train_locs = train_locs, train_field = train_field, train_X = train_X, PP = PP), "Heavy_metals/validation_train.RDS")
  
  # get test data set
  test_locs = locs  [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),]  # get locs
  # save
  saveRDS(list(
    test_locs = test_locs, 
    test_X    = X     [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),],
    test_obs  = field [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))) ]
  ), "Heavy_metals/test_data_set.RDS")
  
}

#################################
# nonstationary noise and range #
#################################
train_data_set= readRDS(file = "Heavy_metals/validation_train.RDS")
for(noise_nonstat in c(TRUE, FALSE)){
for(scale_nonstat in c(TRUE, FALSE)){
  scale_beta = NULL
  if(scale_nonstat ==T)scale_X = train_data_set$train_X[,c("gcarb","globedem","twi")]
  noise_beta = NULL
  if(noise_nonstat ==T)noise_X = train_data_set$train_X
 mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
   observed_locs = train_data_set$train_locs, 
   observed_field = train_data_set$train_field, 
   X = train_data_set$train_X, m = 10, nu = 1.5, anisotropic = F, 
   sphere = T, PP = train_data_set$PP, n_chains = 4, 
   noise_PP = noise_nonstat, noise_X = noise_X, noise_log_scale_prior = c(-10, 10),
   scale_PP = scale_nonstat, scale_X = scale_X, scale_log_scale_prior = c(-10, 10),
 ) 
 for(iter in seq(20))mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(
   mcmc_nngp_list, seed = 2, n_cores = 2, num_threads_per_chain = 10
 )
 saveRDS(mcmc_nngp_list, paste("Heavy_metals/run_noise", noise_nonstat, "scale", scale_nonstat, ".RDS", sep ="_"))
}
}
  
