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
saveRDS(KL_decomposition, "Heavy_metals/KL_decomposition.RDS")

plot(KL_decomposition$d^2)

KL_basis = KL_decomposition$u[match(split(locs, row(locs)), split(locs_, row(locs_))),]

Bidart::plot_pointillist_painting(locs, KL_basis[,1])
#Bidart::plot_pointillist_painting(locs_, (PP_basis%*%KL_decomposition$v[,seq(20)]%*%diag(1/KL_decomposition$d[seq(20)]))[,1])
Bidart::plot_pointillist_painting(locs, KL_basis[,2])
#Bidart::plot_pointillist_painting(locs_, (PP_basis%*%KL_decomposition$v[,seq(20)]%*%diag(1/KL_decomposition$d[seq(20)]))[,2])
Bidart::plot_pointillist_painting(locs, KL_basis[,3])
Bidart::plot_pointillist_painting(locs, KL_basis[,4])
Bidart::plot_pointillist_painting(locs, KL_basis[,20])


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
  seed = 1
)

remove(X);remove(locs);remove(locs_);remove(NNarray);remove(PP_basis);remove(KL_basis);remove(KL_decomposition) 
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, burn_in = .2)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_nsr_basis")




# noise, scale, with spatial basis

X = readRDS("Heavy_metals/processed_data.RDS")$X_locs
X$MAJOR1 = NULL ; X$glwd31 = NULL ; X$minotype = NULL
locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs


KL_decomposition = readRDS("Heavy_metals/KL_decomposition.RDS")



mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary (
  observed_locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs, #spatial locations
  observed_field = readRDS("Heavy_metals/processed_data.RDS")$observed_field, # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_isotropic_sphere", response_model = "Gaussian", # covariance model and response model
  noise_X = as.data.frame(cbind(X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  scale_X = as.data.frame(cbind(X, KL_basis[,seq(20)])), # range for latent field of parameters, if NULL no latent field
  n_chains = 3,  # number of MCMC chains
  seed = 1
)

remove(X);remove(locs);remove(locs_);remove(NNarray);remove(PP_basis);remove(KL_basis);remove(KL_decomposition) 
for(i in seq(40))
{ 
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, field_n_mala = 0, burn_in = .2)
}
saveRDS(mcmc_nngp_list, "Heavy_metals/run_ns_basis")


