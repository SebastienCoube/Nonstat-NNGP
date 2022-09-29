
#################
# Scale + noise #
#################
remove(list = ls()) ; gc()
set.seed(2)
n_locs = 3000
n_obs = 10000

# locations and latent fields
locs = cbind(10 * runif(n_locs), .000001*runif(n_locs))
latent_field = GpGp::fast_Gp_sim(c(1, .005, 1.5, 0), locs = locs)

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
KL_scale = as.data.frame(KL$basis %*% diag(KL$KL_decomposition$d))

# regression coeffs
beta = c(1, .01,  -.01)
beta_noise = c(4, rnorm(ncol(KL_scale))) +  1*c(0, rnorm(ncol(KL_scale)))
beta_scale = c(4, rnorm(ncol(KL_scale))) +  1*c(0, rnorm(ncol(KL_scale)))



# actual fields 
log_noise_variance = as.matrix(cbind(1, KL_noise)) %*% beta_noise
log_scale = as.matrix(cbind(1, KL_scale)) %*% beta_scale 
observed_field = exp(.5*log_scale)* latent_field[observation_idx] +exp(.5*log_noise_variance)*rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta

par(mfrow = c(2, 2))
plot(observed_locs[,1], log_scale)
plot(observed_locs[,1], log_noise_variance)

par(mfrow = c(2, 2))
plot(observed_locs[,1], observed_field, col = 2)
points(observed_locs[,1], exp(.5*log_scale)* latent_field[observation_idx] + cbind(1, as.matrix(X))%*%beta)
plot(observed_locs[,1], latent_field[observation_idx])
plot(observed_locs[,1], exp(log_noise_variance/2))
plot(observed_locs[,1], exp(log_scale/2))


source("Bidart/R/Coloring.R")
source("Bidart/R/Useful_stuff.R")
source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
  observed_locs = observed_locs, #spatial locations
  X = X, observed_field = as.vector(observed_field), 
  m = 10, 
  reordering = "maxmin", covfun = "matern_isotropic", nu = 1.5,
  noise_X = NULL, noise_KL = T,
  scale_X = NULL, scale_KL = T, 
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

#state  = mcmc_nngp_list$states$chain_1
##state$params$scale_log_scale = 0
#source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
#source("Bidart/R/Useful_stuff.R")
#test =      mcmc_nngp_update_Gaussian(data = mcmc_nngp_list$data,
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