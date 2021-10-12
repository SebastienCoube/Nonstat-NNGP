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





# noise, scale, range, with spatial basis

X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs

set.seed(1)
locs_ = unique(locs)
locs_ = locs_[GpGp::order_maxmin(locs_),]
#NNarray = GpGp::find_ordered_nn(locs_, 5)
#PP_basis = Matrix::solve(
#  Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], 
#                       j = NNarray[!is.na(NNarray)], 
#                       x = GpGp::vecchia_Linv(c(1, 3, 1, 0), "matern_isotropic", locs_, NNarray)[!is.na(NNarray)], 
#                       triangular = T
#  ), 
#  diag(1, nrow(locs_), 1000)
#)
#KL_decomposition = irlba::irlba(PP_basis, nu = 40, nv = 40)
#saveRDS(KL_decomposition, "Heavy_metals/KL_decomposition.RDS")
KL_decomposition = readRDS("Heavy_metals/KL_decomposition.RDS")


plot(KL_decomposition$d^2)

KL_basis = KL_decomposition$u[match(split(locs, row(locs)), split(locs_, row(locs_))),]

Bidart::plot_pointillist_painting(locs, KL_basis[,1])
#Bidart::plot_pointillist_painting(locs_, (PP_basis%*%KL_decomposition$v[,seq(20)]%*%diag(1/KL_decomposition$d[seq(20)]))[,1])
Bidart::plot_pointillist_painting(locs, KL_basis[,2])
#Bidart::plot_pointillist_painting(locs_, (PP_basis%*%KL_decomposition$v[,seq(20)]%*%diag(1/KL_decomposition$d[seq(20)]))[,2])
Bidart::plot_pointillist_painting(locs, KL_basis[,3])
Bidart::plot_pointillist_painting(locs, KL_basis[,4])
Bidart::plot_pointillist_painting(locs, KL_basis[,5])
Bidart::plot_pointillist_painting(locs, KL_basis[,6])
Bidart::plot_pointillist_painting(locs, KL_basis[,7])
Bidart::plot_pointillist_painting(locs, KL_basis[,8])
Bidart::plot_pointillist_painting(locs, KL_basis[,9])
Bidart::plot_pointillist_painting(locs, KL_basis[,10])
Bidart::plot_pointillist_painting(locs, KL_basis[,20])
Bidart::plot_pointillist_painting(locs, KL_basis[,25])
Bidart::plot_pointillist_painting(locs, KL_basis[,30])

seed = Sys.time()
par(mfrow = c(2, 1))
set.seed(seed)
Bidart::plot_pointillist_painting(locs, KL_basis[,seq(40)]%*% rnorm(40)%*% KL_decomposition$d^2)
set.seed(seed)
Bidart::plot_pointillist_painting(locs, KL_basis[,seq(30)]%*% rnorm(30)%*% KL_decomposition$d[seq(30)]^2)




mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_isotropic_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  scale_X = as.data.frame(cbind(X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  range_X = as.data.frame(cbind(X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X);remove(locs);remove(locs_);remove(NNarray);remove(PP_basis);remove(KL_basis);remove(KL_decomposition) 
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, burn_in = .2)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_nsr_basis")




# noise, scale, range, only spatial basis

X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs

set.seed(1)
locs_ = unique(locs)
locs_ = locs_[GpGp::order_maxmin(locs_),]
#NNarray = GpGp::find_ordered_nn(locs_, 5)
#PP_basis = Matrix::solve(
#  Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], 
#                       j = NNarray[!is.na(NNarray)], 
#                       x = GpGp::vecchia_Linv(c(1, 3, 1, 0), "matern_isotropic", locs_, NNarray)[!is.na(NNarray)], 
#                       triangular = T
#  ), 
#  diag(1, nrow(locs_), 1000)
#)
#KL_decomposition = irlba::irlba(PP_basis, nu = 40, nv = 40)
#saveRDS(KL_decomposition, "Heavy_metals/KL_decomposition.RDS")
KL_decomposition = readRDS("Heavy_metals/KL_decomposition.RDS")


plot(KL_decomposition$d^2)

KL_basis = KL_decomposition$u[match(split(locs, row(locs)), split(locs_, row(locs_))),]

mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_isotropic_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  scale_X = as.data.frame(KL_basis[,seq(20)]), # range for latent field of parameters, if NULL no latent field
  range_X = as.data.frame(KL_basis[,seq(20)]), # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X);remove(locs);remove(locs_);remove(NNarray);remove(PP_basis);remove(KL_basis);remove(KL_decomposition) 
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, burn_in = .2)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_nsr_basis_only")












##############
# Validation #
##############

# noise, scale, range, with spatial basis and removed observations for validation

# loading data
X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs
field = readRDS("Heavy_metals/processed_data.RDS")$observed_field

# shuffle locs and select 50000 train locs
set.seed(1)
locs_ = unique(locs)
locs_ = locs_[GpGp::order_maxmin(locs_),]
locs_ = locs_[seq(50000),]

# getting KL basis
NNarray = GpGp::find_ordered_nn(locs_, 5)
PP_basis = Matrix::solve(
  Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], 
                       j = NNarray[!is.na(NNarray)], 
                       x = GpGp::vecchia_Linv(c(1, 3, 1, 0), "matern_isotropic", locs_, NNarray)[!is.na(NNarray)], 
                       triangular = T
  ), 
  diag(1, nrow(locs_), 1000)
)
KL_decomposition = irlba::irlba(PP_basis, nu = 40, nv = 40)
saveRDS(KL_decomposition, "Heavy_metals/KL_decomposition_validation.RDS")
KL_decomposition = readRDS("Heavy_metals/KL_decomposition_validation.RDS")

# matching KL basis with selected regressors and obs
observed_locs = locs_              [na.omit(match(split(locs, row(locs)), split(locs_, row(locs_)))),]
KL_basis      = KL_decomposition$u [na.omit(match(split(locs, row(locs)), split(locs_, row(locs_)))),]
observed_field = field [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_))))  ]
observed_X     = X     [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))), ]
# check
#Bidart::plot_pointillist_painting(observed_locs, KL_basis[,1])
#Bidart::plot_pointillist_painting(observed_locs, observed_X[,1])

# save test data
saveRDS(list(
  "test_locs" = locs  [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),], 
  "test_X"    = X     [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),], 
  "test_obs"  = field [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))) ]
), "Heavy_metals/test_data_set.RDS")

mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs  = observed_locs, #spatial locations
  observed_field = observed_field, # Response variable
  X = observed_X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_isotropic_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(observed_X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  scale_X = as.data.frame(cbind(observed_X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  range_X = as.data.frame(cbind(observed_X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X);remove(locs);remove(locs_);remove(NNarray);remove(PP_basis);remove(KL_basis);remove(KL_decomposition); remove(observed_X); remove(observed_locs) 
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, burn_in = .2)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_nsr_basis_validation.RDS")

# noise, scale, range, with only spatial basis and removed observations for validation

# loading data
X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs
field = readRDS("Heavy_metals/processed_data.RDS")$observed_field

# shuffle locs and select 50000 train locs
set.seed(1)
locs_ = unique(locs)
locs_ = locs_[GpGp::order_maxmin(locs_),]
locs_ = locs_[seq(50000),]

# getting KL basis
NNarray = GpGp::find_ordered_nn(locs_, 5)
PP_basis = Matrix::solve(
  Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], 
                       j = NNarray[!is.na(NNarray)], 
                       x = GpGp::vecchia_Linv(c(1, 3, 1, 0), "matern_isotropic", locs_, NNarray)[!is.na(NNarray)], 
                       triangular = T
  ), 
  diag(1, nrow(locs_), 1000)
)
KL_decomposition = irlba::irlba(PP_basis, nu = 40, nv = 40)
saveRDS(KL_decomposition, "Heavy_metals/KL_decomposition_validation.RDS")
KL_decomposition = readRDS("Heavy_metals/KL_decomposition_validation.RDS")

# matching KL basis with selected regressors and obs
observed_locs = locs_              [na.omit(match(split(locs, row(locs)), split(locs_, row(locs_)))),]
KL_basis      = KL_decomposition$u [na.omit(match(split(locs, row(locs)), split(locs_, row(locs_)))),]
observed_field = field [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_))))  ]
observed_X     = X     [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))), ]
# check
#Bidart::plot_pointillist_painting(observed_locs, KL_basis[,1])
#Bidart::plot_pointillist_painting(observed_locs, observed_X[,1])

# save test data
saveRDS(list(
  "test_locs" = locs  [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),], 
  "test_X"    = X     [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),], 
  "test_obs"  = field [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))) ]
), "Heavy_metals/test_data_set.RDS")

mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs  = observed_locs, #spatial locations
  observed_field = observed_field, # Response variable
  X = observed_X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_isotropic_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(observed_X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  scale_X = as.data.frame(KL_basis[,seq(20)]), # range for latent field of parameters, if NULL no latent field
  range_X = as.data.frame(KL_basis[,seq(20)]), # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X);remove(locs);remove(locs_);remove(NNarray);remove(PP_basis);remove(KL_basis);remove(KL_decomposition); remove(observed_X); remove(observed_locs) 
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, burn_in = .2)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_nsr_basis_validation_basis_only.RDS")



# noise, scale, range, with spatial basis, chosen regressors and removed observations for validation

# loading data
X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs
field = readRDS("Heavy_metals/processed_data.RDS")$observed_field

# shuffle locs and select 50000 train locs
set.seed(1)
locs_ = unique(locs)
locs_ = locs_[GpGp::order_maxmin(locs_),]
locs_ = locs_[seq(50000),]

# getting KL basis
NNarray = GpGp::find_ordered_nn(locs_, 5)
PP_basis = Matrix::solve(
  Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], 
                       j = NNarray[!is.na(NNarray)], 
                       x = GpGp::vecchia_Linv(c(1, 3, 1, 0), "matern_isotropic", locs_, NNarray)[!is.na(NNarray)], 
                       triangular = T
  ), 
  diag(1, nrow(locs_), 1000)
)
KL_decomposition = irlba::irlba(PP_basis, nu = 40, nv = 40)
saveRDS(KL_decomposition, "Heavy_metals/KL_decomposition_validation.RDS")
KL_decomposition = readRDS("Heavy_metals/KL_decomposition_validation.RDS")

# matching KL basis with selected regressors and obs
observed_locs = locs_              [na.omit(match(split(locs, row(locs)), split(locs_, row(locs_)))),]
KL_basis      = KL_decomposition$u [na.omit(match(split(locs, row(locs)), split(locs_, row(locs_)))),]
observed_field = field [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_))))  ]
observed_X     = X     [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))), ]
# check
#Bidart::plot_pointillist_painting(observed_locs, KL_basis[,1])
#Bidart::plot_pointillist_painting(observed_locs, observed_X[,1])

# save test data
saveRDS(list(
  "test_locs" = locs  [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),], 
  "test_X"    = X     [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))),], 
  "test_obs"  = field [is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))) ]
), "Heavy_metals/test_data_set.RDS")

mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs  = observed_locs, #spatial locations
  observed_field = observed_field, # Response variable
  X = observed_X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_isotropic_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(observed_X, KL_basis[,seq(25)])), # range for latent field of parameters, if NULL no latent field
  scale_X = as.data.frame(cbind(observed_X[, c("gcarb", "globedem", "twi")], KL_basis[,seq(25)])), # range for latent field of parameters, if NULL no latent field
  range_X = as.data.frame(cbind(observed_X[, c("gcarb", "globedem", "twi")], KL_basis[,seq(25)])), # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)

remove(X);remove(locs);remove(locs_);remove(NNarray);remove(PP_basis);remove(KL_basis);remove(KL_decomposition); remove(observed_X); remove(observed_locs) 
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, burn_in = .2)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_nsr_basis_validation_chosen.RDS")



# stationary with spatial basis and removed observations


# loading data
X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs
field = readRDS("Heavy_metals/processed_data.RDS")$observed_field
# shuffle locs and select 50000 train locs
set.seed(1)
locs_ = unique(locs)
locs_ = locs_[GpGp::order_maxmin(locs_),]
locs_ = locs_[seq(50000),]
# matching KL basis with selected regressors and obs
observed_locs = locs_              [na.omit(match(split(locs, row(locs)), split(locs_, row(locs_)))),]
observed_field = field [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_))))  ]
observed_X              = X     [!is.na(match(split(locs, row(locs)), split(locs_, row(locs_)))), ]
# check
#Bidart::plot_pointillist_painting(observed_locs, X[,1])
mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs  = observed_locs, #spatial locations
  observed_field = observed_field, # Response variable
  X = observed_X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "exponential_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = NULL, # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 2
)
remove(X);remove(locs);remove(locs_);remove(NNarray);remove(PP_basis);remove(KL_basis);remove(KL_decomposition) 
for(i in seq(50))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_stat_validation.RDS")




# noise, scale, with spatial basis

X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs


KL_decomposition = readRDS("Heavy_metals/KL_decomposition.RDS")
set.seed(1)
locs_ = unique(locs)
locs_ = locs_[GpGp::order_maxmin(locs_),]
KL_basis = KL_decomposition$u[match(split(locs, row(locs)), split(locs_, row(locs_))),]



mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "exponential_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  scale_X = as.data.frame(cbind(X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 1
)

remove(X);remove(locs);remove(locs_);remove(NNarray);remove(PP_basis);remove(KL_basis);remove(KL_decomposition) 
for(i in seq(60))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, burn_in = .2)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_ns_basis")


