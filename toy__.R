
#################
# Range + noise #
#################
#remove(list = ls()) ; gc()
set.seed(2)
n_locs = 10000
n_obs = 20000

# locations and latent fields
locs = cbind(10 * runif(n_locs), .000001*runif(n_locs))

# observing observations, with duplicates
observation_idx = sample(seq(n_locs), n_obs, replace =T)
observed_locs = locs[observation_idx,]

# covariates
X = locs[observation_idx,1]
X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5), X))
#X = as.data.frame(scale(X, scale = F, center = T))
colnames(X) = c("binomial", "first_coordinate")
source("Bidart/R/Useful_stuff.R")
KL = get_KL_basis(locs = observed_locs, n_KL = 6, covfun_name = "matern15_isotropic", covparms = c(1, 1, .0001))

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
KL_ = get_KL_basis(locs = locs, n_KL = 6, covfun_name = "matern15_isotropic", covparms = c(1, 1, .0001))
KL_range = as.data.frame(KL_$basis %*% diag(KL_$KL_decomposition$d))


# regression coeffs
beta = c(1000, 10,  -10)
beta_noise = 1*c(-1, 1, .4, rnorm(ncol(KL_noise)))
beta_range = 1*c(-5, .6* rnorm(ncol(KL_noise)))



# actual fields 
log_noise_variance = as.matrix(cbind(1, X, KL_noise)) %*% beta_noise
log_range = as.matrix(cbind(1, KL_range)) %*% beta_range 

NNarray = GpGp::find_ordered_nn(locs, 10)
#latent_field = GpGp::fast_Gp_sim_Linv(
#  Bidart::compute_sparse_chol(covfun_name = "nonstationary_matern_isotropic", range_beta =  matrix(beta_range), NNarray = NNarray, locs = locs, range_X = as.matrix(cbind(1, KL_range)), compute_derivative = F, nu = 1.5)[[1]], 
#  NNarray )
latent_field = GpGp::fast_Gp_sim_Linv(
  Bidart::compute_sparse_chol(covfun_name = "nonstationary_matern_isotropic", range_beta =  matrix(beta_range), NNarray = NNarray, locs = locs, range_X = matrix(1, nrow(locs)), KL = KL_, use_KL = T, compute_derivative = F, nu = 1.5)[[1]], 
  NNarray )

observed_field = latent_field[observation_idx] +exp(.5*log_noise_variance)*rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta


par(mfrow = c(2, 2))
plot(locs[,1], log_range)
plot(observed_locs[,1], log_noise_variance)

par(mfrow = c(2, 2))
plot(observed_locs[,1], observed_field, col = 2)
points(observed_locs[,1], latent_field[observation_idx] + cbind(1, as.matrix(X))%*%beta)
plot(observed_locs[,1], latent_field[observation_idx])
plot(observed_locs[,1], exp(log_noise_variance/2))
plot(locs[,1], exp(log_range))

par(mfrow = c(1, 1))
plot(observed_locs[,1], latent_field[observation_idx])

source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")

mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
  observed_locs = observed_locs, #spatial locations
  X = X, observed_field = as.vector(observed_field), 
  m = 10, 
  reordering = "maxmin", covfun = "nonstationary_matern_isotropic", nu = 1.5,
  noise_X = X, noise_KL = T,
  scale_X = NULL, scale_KL = F, 
  range_X = NULL, range_KL = T, 
  KL = KL, n_chains = 2
)


state = mcmc_nngp_list$states$chain_2
hierarchical_model = mcmc_nngp_list$hierarchical_model
data = mcmc_nngp_list$data
vecchia_approx = mcmc_nngp_list$vecchia_approx

#source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
#mcmc_nngp_update_Gaussian(data = data, hierarchical_model = hierarchical_model, vecchia_approx = vecchia_approx, state = state, n_iterations_update = 100)


mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, seed = 2)


predicted_locs = seq(0, 10, .005);predicted_locs = cbind(predicted_locs, rnorm(length(predicted_locs), 0, .00001))
predicted_locs = predicted_locs[sample(seq(nrow(predicted_locs)), 2*nrow(predicted_locs), T),]

X_range_pred = NULL
X_scale_pred = NULL
burn_in = .5
n_cores = 1
predict_range = F
predict_scale = F




estimations = Bidart::estimate_parameters(mcmc_nngp_list)

#remove(KL)

pred = Bidart::predict_latent_field(
  mcmc_nngp_list = mcmc_nngp_list, predicted_locs = predicted_locs, 
  X_range_pred = NULL, X_scale_pred = NULL, 
  burn_in = .5, n_cores = 3
    )



plot(predicted_locs[,1], pred$summaries$log_range[1,], col =2)
points(locs[,1], log_range)

plot(predicted_locs[,1], pred$summaries$field[1,], col =2, cex = .5)
points(locs[,1], latent_field, cex = .3, pch = 15)


plot(predicted_locs[,1], pred$summaries$field[5,], col =2)



pred_fixed = Bidart::predict_fixed_effects(mcmc_nngp_list, X_pred = X)
plot(
as.matrix(cbind(1, X)) %*% beta,
pred_fixed$summaries[1,,])
abline(a = 1, b = 1, col = 2)


X_noise_pred = cbind(rbinom(nrow(predicted_locs), 1, .5), predicted_locs[,1])
pred_noise = Bidart::predict_noise(mcmc_nngp_list = mcmc_nngp_list, X_noise_pred = X_noise_pred, predicted_locs = predicted_locs)

plot(pred_noise$predicted_locs[,1], pred_noise$summaries[1,,], col =2)
points(observed_locs[,1], log_noise_variance)



