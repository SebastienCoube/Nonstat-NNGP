
#################
# Range + noise #
#################
remove(list = ls()) ; gc()
set.seed(2)
n_locs = 10000
n_obs =  20000

# locations and latent fields
locs = cbind(10 * runif(n_locs), 10 * runif(n_locs))

# observing observations, with duplicates
observation_idx = sample(seq(n_locs), n_obs, replace =T)
observed_locs = locs[observation_idx,]

# covariates
X = locs[observation_idx,1]
X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5), X))
#X = as.data.frame(scale(X, scale = F, center = T))
colnames(X) = c("binomial", "first_coordinate")
source("Bidart/R/Useful_stuff.R")
KL = get_KL_basis(locs = observed_locs, n_KL = 12, covfun_name = "matern15_isotropic", covparms = c(1, 2, .0001))

# IMPORTANT REMEMBER TO IMPLEMENT EFFICIENT KL MULT
#beta = rnorm(50)
#V = rep(0, nrow(KL_scale$sparse_chol))
#V[seq(500)]=KL_scale$KL_decomposition$v %*% beta
#t1 = Sys.time()
#for(i in seq(100))tatato = as.matrix(Matrix::solve(KL_scale$sparse_chol, V, triangular = T))
#Sys.time()-t1
#
#t1 = Sys.time()
#for(i in seq(100))tatato = KL_scale$KL_decomposition$u %*% diag(KL_scale$KL_decomposition$d) %*% beta
#Sys.time()-t1

KL_noise = as.data.frame(KL$basis %*% diag(KL$KL_decomposition$d))
KL_ = get_KL_basis(locs = locs, n_KL = 12, covfun_name = "matern15_isotropic", covparms = c(1, 2, .0001))
KL_range = as.data.frame(KL_$basis %*% diag(KL_$KL_decomposition$d))


# regression coeffs
beta = c(1, .01,  -.01)
beta_noise = 0*c(1, rnorm(ncol(KL_noise)))
beta_range = 1*cbind(c(-1, rnorm(ncol(KL_noise))), c(0, rnorm(ncol(KL_noise))), c(0, rnorm(ncol(KL_noise))))



# actual fields 
log_noise_variance = as.matrix(cbind(1, KL_noise)) %*% beta_noise
log_range = as.matrix(cbind(1, KL_range)) %*% beta_range 

NNarray = GpGp::find_ordered_nn(locs, 10)
#latent_field = GpGp::fast_Gp_sim_Linv(
#  Bidart::compute_sparse_chol(covfun_name = "nonstationary_matern_isotropic", range_beta =  matrix(beta_range), NNarray = NNarray, locs = locs, range_X = as.matrix(cbind(1, KL_range)), compute_derivative = F, nu = 1.5)[[1]], 
#  NNarray )
Bidart::compute_sparse_chol(covfun_name = "nonstationary_matern_anisotropic", range_beta =  beta_range, NNarray = NNarray, locs = locs, range_X = matrix(1, nrow(locs)), KL = KL_, use_KL = T, compute_derivative = F, nu = 1.5)[[1]]
Bidart::compute_sparse_chol(covfun_name = "nonstationary_matern_anisotropic", range_beta =  beta_range, NNarray = NNarray, locs = locs, range_X = matrix(1, nrow(locs)), KL = KL_, use_KL = T, compute_derivative = T, nu = 1.5)[[2]]

latent_field = GpGp::fast_Gp_sim_Linv(
  Bidart::compute_sparse_chol(covfun_name = "nonstationary_exponential_anisotropic", range_beta =  beta_range, NNarray = NNarray, locs = locs, range_X = matrix(1, nrow(locs)), KL = KL_, use_KL = T, compute_derivative = F, nu = 1.5)[[1]], 
  NNarray )

observed_field = latent_field[observation_idx] +exp(.5*log_noise_variance)*rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta

Bidart::plot_pointillist_painting(locs, latent_field)



source("Bidart/R/Coloring.R")
source("Bidart/R/Useful_stuff.R")
source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
  observed_locs = observed_locs, #spatial locations
  X = X, observed_field = as.vector(observed_field), 
  m = 10, 
  reordering = "maxmin", covfun = "nonstationary_matern_anisotropic", nu = 1.5,
  noise_X = NULL, noise_KL = T,
  scale_X = NULL, scale_KL = F, 
  range_X = NULL, range_KL = T, 
  KL = KL
)



state = mcmc_nngp_list$states$chain_2
hierarchical_model = mcmc_nngp_list$hierarchical_model
data = mcmc_nngp_list$data
vecchia_approx = mcmc_nngp_list$vecchia_approx


mcmc_nngp_list$hierarchical_model$beta_priors$range_beta_mean
mcmc_nngp_list$hierarchical_model$beta_priors$noise_beta_mean
mcmc_nngp_list$hierarchical_model$beta_priors$scale_beta_mean
mcmc_nngp_list$hierarchical_model$beta_priors$range_beta_precision
mcmc_nngp_list$hierarchical_model$beta_priors$scale_beta_precision
mcmc_nngp_list$hierarchical_model$beta_priors$noise_beta_precision

#mcmc_nngp_list
#plot(mcmc_nngp_list$data$locs[,1], 
#     field_from_KL_coeffs(mcmc_nngp_list$states$chain_1$params$scale_field, 
#                          mcmc_nngp_list$hierarchical_model$hyperprior_covariance$scale_NNGP_prior, 
#                          mcmc_nngp_list$states$chain_1$params$scale_log_scale)
#)

state  = mcmc_nngp_list$states$chain_1

mcmc_nngp_list$states$chain_1$transition_kernels$range_beta_ancillary 
state$params$range_log_scale = 0
source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
#test =  mcmc_nngp_update_Gaussian(data = mcmc_nngp_list$data,
#                                      hierarchical_model = mcmc_nngp_list$hierarchical_model, vecchia_approx = mcmc_nngp_list$vecchia_approx, # model architecture
#                                      state = state, # model state
#                                      n_iterations_update = 1000, thinning = .2, iter_start = 0, seed = 2# practical settings
#)

source("Bidart/R/visualisation.R")
source("Bidart/R/mcmc_nngp_run.R")
mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 200, n_cycles = 1, seed = 2)


#source("Bidart/R/mcmc_nngp_predict_nonstationary.R")
#test = estimate_parameters(mcmc_nngp_list)
#