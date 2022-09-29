

#################
# Scale + noise #
#################
remove(list = ls()) ; gc()
set.seed(3)
n_locs = 4000
n_obs = 6000

# locations and latent fields
locs = cbind(10 * runif(n_locs), 1)
latent_field = GpGp::fast_Gp_sim(c(1, .01, 1.5, 0), locs = locs)
latent_field_scale = GpGp::fast_Gp_sim(c(1, 1, 1.5, 0), locs = locs)
latent_field_noise = GpGp::fast_Gp_sim(c(1, 1, 1.5, 0), locs = locs)

# observing observations, with duplicates
observation_idx = sample(seq(n_locs), n_obs, replace =T)
observed_locs = locs[observation_idx,]

# covariates
X = locs[observation_idx,1]
X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5), X))
#X = as.data.frame(scale(X, scale = F, center = T))
colnames(X) = c("binomial", "first_coordinate")
X_noise = X
X_scale = X
X_scale[,1] = NULL

# regression coeffs
beta = c(1, .01,  -.01)
beta_noise = c(1, -1, 0)
beta_scale = c(0,0)

# actual fields 
log_noise_variance = as.matrix(cbind(1, X_noise)) %*% beta_noise + latent_field_noise[observation_idx]
log_scale = as.matrix(cbind(1,X_scale)) %*% beta_scale + latent_field_scale[observation_idx]
observed_field = exp(.5*log_scale)* latent_field[observation_idx] +exp(.5*log_noise_variance)*rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta

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
  noise_X = X_noise, noise_KL = T,
  scale_X = X_scale, scale_KL = T,
  KL = KL
)


mcmc_nngp_list
plot(mcmc_nngp_list$data$locs[,1], 
     field_from_KL_coeffs(mcmc_nngp_list$states$chain_1$params$scale_field, 
                          mcmc_nngp_list$hierarchical_model$hyperprior_covariance$scale_NNGP_prior, 
                          mcmc_nngp_list$states$chain_1$params$scale_log_scale)
     )

source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
#mcmc_nngp_list$states$chain_1$params$scale_log_scale =  -.7
#mcmc_nngp_list$states$chain_1$params$noise_log_scale =  -.7
test =      mcmc_nngp_update_Gaussian(data = mcmc_nngp_list$data,
                                      hierarchical_model = mcmc_nngp_list$hierarchical_model, vecchia_approx = mcmc_nngp_list$vecchia_approx, # model architecture
                                      state = mcmc_nngp_list$states$chain_1, # model state
                                      n_iterations_update = 1000, thinning = .2, iter_start = 0, seed = 2,# practical settings
                                      field_n_chromatic = 1,  # number of chromatic steps for the latent field at each iteration
                                      field_n_mala = 1,  # number of mala steps for the latent field at each iteration
                                      )

source("Bidart/R/visualisation.R")
mcmc_nngp_list = mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1, seed = 2)


source("Bidart/R/mcmc_nngp_predict_nonstationary.R")
test = estimate_parameters(mcmc_nngp_list)


predicted_locs =  cbind(seq(1000)/100, 1) 
test = predict_latent_field(
  mcmc_nngp_list = mcmc_nngp_list, 
  predicted_locs = predicted_locs, 
  X_scale_pred = cbind( predicted_locs[,-2]),  
  burn_in = .5, n_cores = 3)



dev.off()
plot(observed_locs[,1], exp(.5*log_scale)* latent_field[observation_idx])
lines(predicted_locs[,1], test$summaries$field[1,,], col = 2)

plot(locs[,1], latent_field_scale)
lines(predicted_locs[,1], test$summaries$scale_field[1,,], col = 2)

plot(observed_locs[,1], log_scale)
lines(predicted_locs[,1], test$summaries$scale_first_coordinate[1,,], col = 2)
lines(predicted_locs[,1], test$summaries$`scale_(Intercept)`[1,,], col = 2)
lines(predicted_locs[,1], test$summaries$log_scale[1,,], col = 2)

source("Bidart/R/mcmc_nngp_predict_nonstationary.R")
test = mcmc_nngp_estimate(mcmc_nngp_list = mcmc_nngp_list)
test$range_beta
test$scale_beta


source("Bidart/R/mcmc_nngp_predict_nonstationary.R")
test = predict_fixed_effects(mcmc_nngp_list = mcmc_nngp_list, predicted_locs = predicted_locs, X_pred = cbind(rbinom(nrow(predicted_locs), 1, .5), predicted_locs[,-2]))
test$summaries$summed_fised_effects
#########
# rang + noise 1D
#########

source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
source("Bidart/R/Coloring.R")


set.seed(1)


# locations and latent fields
#locs = as.matrix(expand.grid(seq(0, 5, .05), seq(0, 5, .05)))

locs = cbind(10 * runif(n_locs), 0.00001* runif(n_locs))
predicted_locs =  as.matrix(expand.grid(seq(0, 5, .1), seq(0, 5, .1)) )

latent_field_scale =  0* GpGp::fast_Gp_sim(c(.5, 1, 1.5, 0), locs = locs)
latent_field_noise =     GpGp::fast_Gp_sim(c(.5, 1, 1.5, 0), locs = locs)
latent_field_range =     GpGp::fast_Gp_sim(c(.5, 1, 1.5, 0), locs = locs)

n_obs = 2 * nrow(locs)
# observing observations, with duplicates
observation_idx = c(sample(nrow(locs), nrow(locs) , replace =F), sample(nrow(locs), n_obs - nrow(locs) , replace =T))
observed_locs = locs[observation_idx,]

# covariates
X = locs[observation_idx,]
X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5), X))
X = as.data.frame(scale(X, center = T, scale = F))
colnames(X) = c("binomial", "first_coordinate", "second_coordinate")
X_noise = X
X_scale = X
X_scale[,1] = NULL
X_range = X_scale


# regression coeffs
beta = c(100, 10,  -10, 5)
beta_noise = c(.5, 1, .1, .1)
beta_scale = c(1,.00, .00)
beta_range = c(-4,.0, .0)

# actual fields 
log_noise_variance = as.matrix(cbind(1, X_noise)) %*% beta_noise + latent_field_noise[observation_idx]
log_scale = as.matrix(cbind(1,X_scale)) %*% beta_scale + latent_field_scale[observation_idx]
log_range = as.matrix(cbind(1, X_range)) %*% beta_range + latent_field_range[observation_idx]

NNarray = GpGp::find_ordered_nn(locs, 10)
sparse_chol = Bidart::compute_sparse_chol(covfun_name = "nonstationary_matern_isotropic", range_beta = beta_range, locs = locs, NNarray = NNarray, range_field =  latent_field_range, range_X = as.matrix(cbind(1, locs)), nu = 1.5)
latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)


observed_field = as.vector(exp(.5*log_scale)* latent_field[observation_idx] +exp(.5*log_noise_variance)*rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta)

par(mfrow = c(2, 2))
plot(observed_locs[,1], observed_field, col = 2)
points(observed_locs[,1], exp(.5*log_range)* latent_field[observation_idx] + cbind(1, as.matrix(X))%*%beta)
plot(observed_locs[,1], latent_field[observation_idx])
plot(observed_locs[,1], exp(log_noise_variance/2))
plot(observed_locs[,1], exp(log_range/2))

source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
source("Bidart/R/Coloring.R")
source("Bidart/R/Useful_stuff.R")
mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
  observed_locs = observed_locs, #spatial locations
  observed_field = c(observed_field), # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_matern_isotropic", response_model = "Gaussian", nu = 1.5,# covariance model and response model
  noise_X = X_noise, noise_range = 1, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, range_range = 1, # range for latent field of parameters, if NULL no latent field
  log_NNGP_matern_covfun = "matern_isotropic", # covariance function for the hyperpriors
  log_NNGP_nu = 1.5, # covariance function for the hyperpriors
  n_chains = 3,  # number of MCMC chains
  seed = 4, n_KL = 100, n_PP = 1000
)


source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
source("Bidart/R/mcmc_nngp_run.R")
source("Bidart/R/visualisation.R")
source("Bidart/R/Useful_stuff.R")
source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
#mcmc_nngp_list$states$chain_1$params$range_log_scale = -.7
test =      mcmc_nngp_update_Gaussian(data = mcmc_nngp_list$data,
                                      hierarchical_model = mcmc_nngp_list$hierarchical_model, vecchia_approx = mcmc_nngp_list$vecchia_approx, # model architecture
                                      state = mcmc_nngp_list$states$chain_1, # model state
                                      n_iterations_update = 1000, thinning = .2, iter_start = 0, seed = 2,# practical settings
                                      field_n_chromatic = 1,  # number of chromatic steps for the latent field at each iteration
                                      field_n_mala = 1,  # number of mala steps for the latent field at each iteration
)
#########
# On 2D #
#########

source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
source("Bidart/R/Coloring.R")


set.seed(1)


# locations and latent fields
#locs = as.matrix(expand.grid(seq(0, 5, .05), seq(0, 5, .05)))
locs = 5*matrix(runif(20000), 10000)
predicted_locs =  as.matrix(expand.grid(seq(0, 5, .1), seq(0, 5, .1)) )

latent_field_scale =  0* GpGp::fast_Gp_sim(c(.5, 1, 1.5, 0), locs = locs)
latent_field_noise =     GpGp::fast_Gp_sim(c(.5, 1, 1.5, 0), locs = locs)
latent_field_range =     GpGp::fast_Gp_sim(c(.5, 1, 1.5, 0), locs = locs)

n_obs = 2 * nrow(locs)
# observing observations, with duplicates
observation_idx = c(sample(nrow(locs), nrow(locs) , replace =F), sample(nrow(locs), n_obs - nrow(locs) , replace =T))
observed_locs = locs[observation_idx,]

# covariates
X = locs[observation_idx,]
X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5), X))
X = as.data.frame(scale(X, center = T, scale = F))
colnames(X) = c("binomial", "first_coordinate", "second_coordinate")
X_noise = X
X_scale = X
X_scale[,1] = NULL
X_range = X_scale


# regression coeffs
beta = c(100, 10,  -10, 5)
beta_noise = c(.5, 1, .1, .1)
beta_scale = c(1,.00, .00)
beta_range = c(-4,.0, .0)

# actual fields 
log_noise_variance = as.matrix(cbind(1, X_noise)) %*% beta_noise + latent_field_noise[observation_idx]
log_scale = as.matrix(cbind(1, X_scale)) %*% beta_scale + latent_field_scale[observation_idx]
log_range = as.matrix(cbind(1, X_range)) %*% beta_range + latent_field_range[observation_idx]

NNarray = GpGp::find_ordered_nn(locs, 10)
sparse_chol = Bidart::compute_sparse_chol(covfun_name = "nonstationary_matern_isotropic", range_beta = beta_range, locs = locs, NNarray = NNarray, range_field =  latent_field_range, range_X = as.matrix(cbind(1, locs)), nu = 1.5)
latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)


observed_field = as.vector(exp(.5*log_scale)* latent_field[observation_idx] +exp(.5*log_noise_variance)*rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta)

get_colors = function(x)heat.colors(100)[round((x - min(x))/(max(x)-min(x))*100)+1]
par(mfrow = c(2, 1))

plot(locs[,], col = get_colors(latent_field_range), main = "range latent field", pch = 15)
plot(observed_locs[,], col = get_colors((log_range)), main = "log range field", pch = 15)

plot(locs[,], col = get_colors(latent_field), main = "latent field", pch = 15)
plot(observed_locs[,], col = get_colors(observed_field), main = "observed field", pch = 15)

plot(locs[,], col = get_colors(latent_field_noise), main = "noise latent field", pch = 15)
plot(observed_locs[,], col = get_colors(log_noise_variance), main = "log noise variance", pch = 15)

source("Bidart/R/visualisation.R")
par(mfrow = c(1, 1))
plot_pointillist_painting(locs, latent_field)
plot_pointillist_painting(observed_locs, log_range)
plot_pointillist_painting(locs, latent_field_noise)

source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
source("Bidart/R/Coloring.R")
source("Bidart/R/Useful_stuff.R")
mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
  observed_locs = observed_locs, #spatial locations
  observed_field = c(observed_field), # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_matern_isotropic", response_model = "Gaussian", nu = 1.5,# covariance model and response model
  noise_X = X_noise, noise_range = 1, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, range_range = 1, # range for latent field of parameters, if NULL no latent field
  log_NNGP_matern_covfun = "matern_isotropic", # covariance function for the hyperpriors
  log_NNGP_nu = 1.5, # covariance function for the hyperpriors
  n_chains = 3,  # number of MCMC chains
  seed = 4, n_KL = 100, n_PP = 1000
)


source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
source("Bidart/R/mcmc_nngp_run.R")
source("Bidart/R/visualisation.R")
source("Bidart/R/Useful_stuff.R")
source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
#mcmc_nngp_list$states$chain_1$params$range_log_scale = -.7
##test =      mcmc_nngp_update_Gaussian(data = mcmc_nngp_list$data,
##                                      hierarchical_model = mcmc_nngp_list$hierarchical_model, vecchia_approx = mcmc_nngp_list$vecchia_approx, # model architecture
##                                      state = mcmc_nngp_list$states$chain_1, # model state
##                                      n_iterations_update = 1000, thinning = .2, iter_start = 0, seed = 2,# practical settings
##                                      field_n_chromatic = 1,  # number of chromatic steps for the latent field at each iteration
##                                      field_n_mala = 1,  # number of mala steps for the latent field at each iteration
##)


mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_iterations_update = 100, n_cycles = 1)

#########
# On 2D anisotropic range random #
#########


remove(list = ls())
gc()
#Rcpp::sourceCpp("Bidart/R/vecchia_nonstat.cpp")
source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
source("Bidart/R/Coloring.R")


set.seed(2)


# locations and latent fields

locs = as.matrix(expand.grid(seq(0, 5, .05), seq(0, 5, .05)))
locs = locs[GpGp::order_maxmin(locs),]

latent_field_scale = 0 * GpGp::fast_Gp_sim(c(.5, 1, .5, 0), locs = locs)
latent_field_noise = 0 * GpGp::fast_Gp_sim(c(.5, 1, .5, 0), locs = locs)
latent_field_range =  
  cbind(GpGp::fast_Gp_sim(c(.5, 1, 1, 0), locs = locs), 
        GpGp::fast_Gp_sim(c(.5, 1, 1, 0), locs = locs), 
        GpGp::fast_Gp_sim(c(.5, 1, 1, 0), locs = locs))


n_obs = nrow(locs) + 10000
# observing observations, with duplicates
observation_idx = c(sample(nrow(locs), nrow(locs) , replace =F), sample(nrow(locs), n_obs - nrow(locs) , replace =T))
observed_locs = locs[observation_idx,]

# covariates
X = locs[observation_idx,]
X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5), X))
X = as.data.frame(scale(X, center = T, scale = F))
colnames(X) = c("binomial", "first_coordinate", "second_coordinate")
X_noise = X
X_scale = X
X_scale[,1] = NULL
X_range = X_scale


# regression coeffs
beta = c(100, 10,  -10, 5)
beta_noise = c(-5, 0, 0, 0)
beta_scale = c(0,.00, .00)
beta_range = matrix(c(-1,.0, .0, -1,.0, .0, 0 ,.0, .0), ncol = 3)

# actual fields 
log_noise_variance = as.matrix(cbind(1, X_noise)) %*% beta_noise + latent_field_noise[observation_idx]
log_scale = as.matrix(cbind(1,X_scale)) %*% beta_scale + latent_field_scale[observation_idx]
log_range = as.matrix(cbind(1, X_range)) %*% beta_range + latent_field_range[observation_idx,]

NNarray = GpGp::find_ordered_nn(locs, 10)
sparse_chol = Bidart::compute_sparse_chol(covfun_name = "nonstationary_exponential_anisotropic", range_beta = beta_range, locs = locs, NNarray = NNarray, range_field =  latent_field_range, range_X = as.matrix(cbind(1, locs)))
latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)




observed_field = 
  as.vector(exp(.5*log_scale)* latent_field[observation_idx]
            +exp(.5*log_noise_variance)*rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta
  )




source("Bidart/R/visualisation.R")
par(mfrow = c(1, 1))
#plot_ellipses(locs, latent_field_range, shrink = .005)

#plot_ellipses(locs, latent_field_range + cbind(1, locs) %*% beta_range, shrink = .05)

plot_pointillist_painting(locs = locs, field = latent_field_range[,1]+ latent_field_range[,2])
par(mfrow = c(2,2))
plot_pointillist_painting(locs = locs, field = latent_field, cex = .8)
plot_pointillist_painting(locs = locs, field = latent_field_range[,1])
plot_pointillist_painting(locs = locs, field = latent_field_range[,2])
plot_pointillist_painting(locs = locs, field = latent_field_range[,3])



source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
source("Bidart/R/Useful_stuff.R")
mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
  observed_locs = observed_locs, #spatial locations
  observed_field = c(observed_field), # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_anisotropic", response_model = "Gaussian", # covariance model and response model
  noise_X = NULL, noise_range = NULL, # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, scale_range = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, range_range = .5, # range for latent field of parameters, if NULL no latent field
  log_NNGP_matern_covfun = "matern_isotropic", # covariance function for the hyperpriors
  log_NNGP_nu = 1, # covariance function for the hyperpriors
  n_chains = 3,  # number of MCMC chains
  seed = 4, n_PP = 1000, n_KL = 100
)

gc()

source("Bidart/R/mcmc_nngp_run.R")
source("Bidart/R/visualisation.R")
source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")

test =      mcmc_nngp_update_Gaussian(data = mcmc_nngp_list$data,
                                      hierarchical_model = mcmc_nngp_list$hierarchical_model, vecchia_approx = mcmc_nngp_list$vecchia_approx, # model architecture
                                      state = mcmc_nngp_list$states$chain_1, # model state
                                      n_iterations_update = 1000, thinning = .2, iter_start = 0, seed = 2,# practical settings
                                      field_n_chromatic = 1,  # number of chromatic steps for the latent field at each iteration
                                      field_n_mala = 1,  # number of mala steps for the latent field at each iteration
)


for(i in seq(10))mcmc_nngp_list = mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1)

plot(mcmc_nngp_list$records$chain_1$range_log_scale[1,1,])
points(mcmc_nngp_list$records$chain_2$range_log_scale[1,1,], col = 2)
points(mcmc_nngp_list$records$chain_3$range_log_scale[1,1,], col = 3)


plot(mcmc_nngp_list$records$chain_1$range_log_scale[6,,])
points(mcmc_nngp_list$records$chain_2$range_log_scale[6,,], col = 2)
points(mcmc_nngp_list$records$chain_3$range_log_scale[6,,], col = 3)

plot(apply(mcmc_nngp_list$records$chain_1$range_log_scale, 3, Bidart::expmat)[1,])
points(apply(mcmc_nngp_list$records$chain_2$range_log_scale, 3, Bidart::expmat)[1,], col = 2)
points(apply(mcmc_nngp_list$records$chain_3$range_log_scale, 3, Bidart::expmat)[1,], col = 3)
plot(apply(mcmc_nngp_list$records$chain_1$range_log_scale, 3, Bidart::expmat)[2,])
points(apply(mcmc_nngp_list$records$chain_2$range_log_scale, 3, Bidart::expmat)[2,], col = 2)
points(apply(mcmc_nngp_list$records$chain_3$range_log_scale, 3, Bidart::expmat)[2,], col = 3)



plot(apply(mcmc_nngp_list$records$chain_1$range_log_scale, 3, Bidart::expmat)[3,])
points(apply(mcmc_nngp_list$records$chain_2$range_log_scale, 3, Bidart::expmat)[3,], col = 2)
points(apply(mcmc_nngp_list$records$chain_3$range_log_scale, 3, Bidart::expmat)[3,], col = 3)



plot(apply(mcmc_nngp_list$records$chain_1$range_log_scale, 3, Bidart::expmat)[5,])
points(apply(mcmc_nngp_list$records$chain_2$range_log_scale, 3, Bidart::expmat)[5,], col = 2)
points(apply(mcmc_nngp_list$records$chain_3$range_log_scale, 3, Bidart::expmat)[5,], col = 3)
plot(apply(mcmc_nngp_list$records$chain_1$range_log_scale, 3, Bidart::expmat)[9,])
points(apply(mcmc_nngp_list$records$chain_2$range_log_scale, 3, Bidart::expmat)[9,], col = 2)
points(apply(mcmc_nngp_list$records$chain_3$range_log_scale, 3, Bidart::expmat)[9,], col = 3)




plot_pointillist_painting(locs = locs, field = latent_field, cex = 1)

plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_3$field[,,200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_2$field[,,200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_1$field[,,200], cex = 1)


par(mfrow = c(1, 1))
plot_pointillist_painting(locs = locs, field = latent_field_range[,1])
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_1$range_field[,1,200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_2$range_field[,1,200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_3$range_field[,1,200], cex = 1)

plot_pointillist_painting(locs = locs, field = latent_field_range[,2])
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_1$range_field[,2, 200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_2$range_field[,2, 200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_3$range_field[,2, 200], cex = 1)


plot_pointillist_painting(locs = locs, field = latent_field_range[,1]-latent_field_range[,2] )
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_1$range_field[,1,200] -  mcmc_nngp_list$records$chain_1$range_field[,2,200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_2$range_field[,1,200] -  mcmc_nngp_list$records$chain_2$range_field[,2,200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_3$range_field[,1,200] -  mcmc_nngp_list$records$chain_3$range_field[,2,200], cex = 1)


plot_pointillist_painting(locs = locs, field = latent_field_range[,3])
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_1$range_field[,3,200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_2$range_field[,3,200], cex = 1)
plot_pointillist_painting(locs = mcmc_nngp_list$data$locs, field = mcmc_nngp_list$records$chain_3$range_field[,3,200], cex = 1)



source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
test =      mcmc_nngp_update_Gaussian(data = mcmc_nngp_list$data,
                                      hierarchical_model = mcmc_nngp_list$hierarchical_model, vecchia_approx = mcmc_nngp_list$vecchia_approx, # model architecture
                                      state = mcmc_nngp_list$states$chain_1, # model state
                                      n_iterations_update = 2500, thinning = .1, iter_start = 1, seed = 2,# practical settings
                                      field_n_chromatic = 2,  # number of chromatic steps for the latent field at each iteration
)
source("Bidart/R/mcmc_nngp_run.R")
test = mcmc_nngp_run_nonstationary("mcmc_nngp_list")

par(mfrow = c(2, 1))
plot(c(test$params_records$range_log_scale)[-seq(50)])
abline(h = -log(2))
plot(cumsum (c(test$params_records$range_log_scale)[-seq(50)] - mean(c(test$params_records$range_log_scale)[-seq(50)])))

plot(c(test$params_records$noise_log_scale)[-seq(50)] )
abline(h = -log(2))
plot(cumsum (c(test$params_records$noise_log_scale)[-seq(50)] - mean(c(test$params_records$noise_log_scale)[-seq(50)])))

plot(test$params_records$range_beta[1,, -seq(50)] )
abline(h = beta_range[1])
plot(cumsum (test$params_records$range_beta[1,, -seq(50)] - mean(test$params_records$range_beta[1,, -seq(50)])))
plot(test$params_records$range_beta[2,, -seq(50)] )
abline(h = beta_range[2])
plot(cumsum (test$params_records$range_beta[2,, -seq(50)] - mean(test$params_records$range_beta[2,, -seq(50)])))
plot(test$params_records$range_beta[3,, -seq(50)] )
abline(h = beta_range[3])
plot(cumsum (test$params_records$range_beta[3,, -seq(50)] - mean(test$params_records$range_beta[3,, -seq(50)])))


par(mfrow = c(2, 2))
plot(
  mcmc_nngp_list$data$locs, 
  col = get_colors(apply(test$params_records$noise_field[,,- seq(50)],1 , function(x)quantile(x, .025))),
  ylab = ""
)
plot(
  mcmc_nngp_list$data$locs, 
  col = get_colors(apply(test$params_records$noise_field[,,- seq(50)],1 , function(x)quantile(x, .975))),
  ylab = ""
)
plot(locs[,], col = get_colors(latent_field_noise), main = "noise latent field")
plot(
  mcmc_nngp_list$data$locs, 
  col = 
    1 + 
    (apply(test$params_records$noise_field[,,- seq(50)],1 , function(x)quantile(x, .025))>latent_field_noise) +
    2 * (apply(test$params_records$noise_field[,,- seq(50)],1 , function(x)quantile(x, .975))<latent_field_noise)
  ,
  ylab = ""
)




#########
# On 2D , iso data, aniso model #
#########
remove(list = ls())
gc()
Rcpp::sourceCpp("Bidart/R/vecchia_nonstat.cpp")
source("Bidart/R/Coloring.R")
source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")


set.seed(2)


# locations and latent fields
locs = as.matrix(expand.grid(seq(0, 5, .05), seq(0, 5, .05)))

locs = locs[GpGp::order_maxmin(locs),]

latent_field_scale = 0 * GpGp::fast_Gp_sim(c(.1, .5, 1, 0), locs = locs)
latent_field_noise =  GpGp::fast_Gp_sim(c(.5, .5, 1, 0), locs = locs)
latent_field_range =  GpGp::fast_Gp_sim(c(.5, .5, 1, 0), locs = locs)


n_obs = nrow(locs) + 10000
# observing observations, with duplicates
observation_idx = c(sample(nrow(locs), nrow(locs) , replace =F), sample(nrow(locs), n_obs - nrow(locs) , replace =F))
observed_locs = locs[observation_idx,]

# covariates
X = locs[observation_idx,]
X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5), X))
X = as.data.frame(scale(X, center = T, scale = F))
colnames(X) = c("binomial", "first_coordinate", "second_coordinate")
X_noise = X
X_scale = X
X_scale[,1] = NULL
X_range = X_scale


# regression coeffs
beta = c(100, 10,  -10, 5)
beta_noise = c(-1, 1, -.1, .1)
beta_scale = c(1,.00, .00)
beta_range = c(-1.5,.0, .0)

# actual fields 
log_noise_variance = as.matrix(cbind(1, X_noise)) %*% beta_noise + latent_field_noise[observation_idx]
log_scale = as.matrix(cbind(1,X_scale)) %*% beta_scale + latent_field_scale[observation_idx]
log_range = as.matrix(cbind(1, X_range)) %*% beta_range + latent_field_range[observation_idx]

NNarray = GpGp::find_ordered_nn(locs, 10)
sparse_chol = Bidart::compute_sparse_chol(covfun_name = "nonstationary_exponential_isotropic", range_beta = beta_range, locs = locs, NNarray = NNarray, range_field =  latent_field_range, range_X = as.matrix(cbind(1, locs)))
latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)


observed_field = as.vector(exp(.5*log_scale)* latent_field[observation_idx] +exp(.5*log_noise_variance)*rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta)

get_colors = function(x)heat.colors(100)[round((x - min(x))/(max(x)-min(x))*100)+1]
par(mfrow = c(1, 1))

plot(locs[,], col = get_colors(latent_field_range), main = "range latent field", pch = 15)
plot(observed_locs[,], col = get_colors((log_range)), main = "log range field", pch = 15)

plot(locs[,], col = get_colors(latent_field), main = "latent field", pch = 15)
plot(observed_locs[,], col = get_colors(observed_field), main = "observed field", pch = 15)

plot(locs[,], col = get_colors(latent_field_noise), main = "noise latent field", pch = 15)
plot(observed_locs[,], col = get_colors(log_noise_variance), main = "log noise variance", pch = 15)

source("Bidart/R/visualisation.R")
par(mfrow = c(1, 1))
plot_ellipses(locs, log_range, shrink = .01)
plot_ellipses(locs, as.matrix(latent_field_range), shrink = .005)

source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
  observed_locs = observed_locs, #spatial locations
  observed_field = c(observed_field), # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_anisotropic", response_model = "Gaussian", # covariance model and response model
  noise_X = X_noise, noise_range = .5, # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, scale_range = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = NULL, range_range = .5, # range for latent field of parameters, if NULL no latent field
  log_NNGP_matern_covfun = "matern_isotropic", # covariance function for the hyperpriors
  log_NNGP_matern_smoothness = 1, # covariance function for the hyperpriors
  n_chains = 3,  # number of MCMC chains
  seed = 4
)

mcmc_nngp_list$states$chain_1$params$range_beta

source("Bidart/R/mcmc_nngp_run.R")
source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
mcmc_nngp_list = mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1)




source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
log_range = cbind(log_range, log_range, 0)
test =      mcmc_nngp_update_Gaussian(data = mcmc_nngp_list$data,
                                      hierarchical_model = mcmc_nngp_list$hierarchical_model, vecchia_approx = mcmc_nngp_list$vecchia_approx, # model architecture
                                      state = mcmc_nngp_list$states$chain_1, # model state
                                      n_iterations_update = 2500, thinning = .1, iter_start = 1, seed = 2,# practical settings
                                      field_n_chromatic = 2,  # number of chromatic steps for the latent field at each iteration
                                      field_n_mala = 2,  # number of mala steps for the latent field at each iteration
                                      scale_n_chromatic_sufficient = 0, scale_n_mala_sufficient =1, scale_n_chromatic_ancillary = 0, scale_n_mala_ancillary = 1,  
                                      noise_n_chromatic = 0,            noise_n_mala =2, 
                                      range_sufficient = T, range_ancillary =F)


M = matrix(rnorm(100), 10)
M[,11] = seq(10)

par(mfrow = c(2, 1))
plot(c(test$params_records$range_log_scale)[-seq(50)])
abline(h = -log(2))
plot(cumsum (c(test$params_records$range_log_scale)[-seq(50)] - mean(c(test$params_records$range_log_scale)[-seq(50)])))

plot(c(test$params_records$noise_log_scale)[-seq(50)] )
abline(h = -log(2))
plot(cumsum (c(test$params_records$noise_log_scale)[-seq(50)] - mean(c(test$params_records$noise_log_scale)[-seq(50)])))

plot(test$params_records$range_beta[1,, -seq(50)] )
abline(h = beta_range[1])
plot(cumsum (test$params_records$range_beta[1,, -seq(50)] - mean(test$params_records$range_beta[1,, -seq(50)])))
plot(test$params_records$range_beta[2,, -seq(50)] )
abline(h = beta_range[2])
plot(cumsum (test$params_records$range_beta[2,, -seq(50)] - mean(test$params_records$range_beta[2,, -seq(50)])))
plot(test$params_records$range_beta[3,, -seq(50)] )
abline(h = beta_range[3])
plot(cumsum (test$params_records$range_beta[3,, -seq(50)] - mean(test$params_records$range_beta[3,, -seq(50)])))


par(mfrow = c(2, 2))
plot(
  mcmc_nngp_list$data$locs, 
  col = get_colors(apply(test$params_records$noise_field[,,- seq(50)],1 , function(x)quantile(x, .025))),
  ylab = ""
)
plot(
  mcmc_nngp_list$data$locs, 
  col = get_colors(apply(test$params_records$noise_field[,,- seq(50)],1 , function(x)quantile(x, .975))),
  ylab = ""
)
plot(locs[,], col = get_colors(latent_field_noise), main = "noise latent field")
plot(
  mcmc_nngp_list$data$locs, 
  col = 
    1 + 
    (apply(test$params_records$noise_field[,,- seq(50)],1 , function(x)quantile(x, .025))>latent_field_noise) +
    2 * (apply(test$params_records$noise_field[,,- seq(50)],1 , function(x)quantile(x, .975))<latent_field_noise)
  ,
  ylab = ""
)






#####################################
# On 2D anisotropic range no random #
#####################################

set.seed(1)


# locations and latent fields
locs = as.matrix(expand.grid(seq(0, 5, .05), seq(0, 5, .05)))


n_obs = nrow(locs) + 1000
# observing observations, with duplicates
observation_idx = c(sample(nrow(locs), nrow(locs) , replace =F), sample(nrow(locs), n_obs - nrow(locs) , replace =T))
observed_locs = locs[observation_idx,]

# covariates
X = locs[observation_idx,]
X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5), X))
X = as.data.frame(scale(X, center = T, scale = F))
colnames(X) = c("binomial", "first_coordinate", "second_coordinate")
X_noise = X
X_scale = X
X_scale[,1] = NULL
X_range = X_scale


# regression coeffs
beta = c(100, 10,  -10, 5)
beta_noise = c(.5, 1, .1, .1)
beta_scale = c(1,.00, .00)
beta_range = matrix(c(-2,.2, -.2, -2, -.2, .2, 0 ,.2, -.2), ncol = 3)

# actual fields 
log_range = as.matrix(cbind(1, X_range)) %*% beta_range + latent_field_range[observation_idx]

NNarray = GpGp::find_ordered_nn(locs, 10)
sparse_chol = Bidart::compute_sparse_chol(covfun_name = "nonstationary_exponential_anisotropic", range_beta = beta_range, locs = locs, NNarray = NNarray, range_field =  NULL, range_X = as.matrix(cbind(1, locs)))
latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)
log_noise_variance = as.matrix(cbind(1, X_noise)) %*% beta_noise


observed_field = 
  as.vector(latent_field[observation_idx]
            +exp(.5*log_noise_variance)* rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta
  )




source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
source("Bidart/R/visualisation.R")
par(mfrow = c(1, 1))

plot_ellipses(locs, cbind(1, locs) %*% beta_range, shrink = .02)

plot_pointillist_painting(locs = locs, field = latent_field, cex = 2)

plot_pointillist_painting(locs = observed_locs, field = observed_field)



source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
  observed_locs = observed_locs, #spatial locations
  observed_field = c(observed_field), # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_anisotropic", response_model = "Gaussian", # covariance model and response model
  noise_X = X_noise, noise_range = NULL, # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, scale_range = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = X_range, range_range = NULL, # range for latent field of parameters, if NULL no latent field
  log_NNGP_matern_covfun = "matern_isotropic", # covariance function for the hyperpriors
  log_NNGP_matern_smoothness = 1, # covariance function for the hyperpriors
  n_chains = 3,  # number of MCMC chains
  seed = 4
)

source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
source("Bidart/R/mcmc_nngp_run.R")
mcmc_nngp_list = mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1)


#####################################
# On 2D anisotropic range no random, with spatial basis#
#####################################
remove(list = ls()); gc()
set.seed(1)


# locations and latent fields
locs = as.matrix(expand.grid(seq(0, 5, .05), seq(0, 5, .05)))+0.001
locs = locs[GpGp::order_maxmin(locs),]


n_obs = nrow(locs) + 10000
# observing observations, with duplicates
observation_idx = c(sample(nrow(locs), nrow(locs) , replace =F), sample(nrow(locs), n_obs - nrow(locs) , replace =T))
observed_locs = locs[observation_idx,]

# covariates
X = locs[observation_idx,]
X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5), X))
X = as.data.frame(scale(X, center = T, scale = F))
colnames(X) = c("binomial", "first_coordinate", "second_coordinate")
X_noise = X



knots = expand.grid(seq(0, 5, 1), seq(0, 5, 1))
NNarray = GpGp::find_ordered_nn(rbind(as.matrix(knots), locs), 10)
X_range = Matrix::solve(Matrix::sparseMatrix(
  i = row(NNarray)[!is.na(NNarray)], 
  j = NNarray[!is.na(NNarray)], 
  x = GpGp::vecchia_Linv(c(1, 1, 0), "exponential_isotropic", rbind(as.matrix(knots), locs), NNarray)[!is.na(NNarray)], 
  triangular = T
), 
Matrix::sparseMatrix(i = seq(nrow(knots)), j = seq(nrow(knots)), x = 1, dims =  c(nrow(rbind(as.matrix(knots), locs)), nrow(knots)))
)[-seq(nrow(knots)),]
X_range = as.matrix(scale(X_range, center = T, scale = F))
plot(GpGp::exponential_isotropic(c(1, 1, 0), cbind(seq(200)/100, 1))[,1])
# regression coeffs
beta = c(5, 10,  -10, 5)
beta_noise = c(-1, .1, -.1, -.1)
beta_range = cbind(c(-2.4, .7 * rnorm(nrow(knots))), c(-2.4, .7 * rnorm(nrow(knots))), c(0, .7 * rnorm(nrow(knots))))


NNarray = GpGp::find_ordered_nn(locs, 10)
sparse_chol = Bidart::compute_sparse_chol(covfun_name = "nonstationary_exponential_anisotropic", range_beta = beta_range, locs = locs, NNarray = NNarray, range_field =  NULL, range_X = as.matrix(cbind(1, X_range)))
latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)
log_noise_variance = as.matrix(cbind(1, X_noise)) %*% beta_noise


observed_field = 
  as.vector(latent_field[observation_idx]
            +exp(.5*log_noise_variance)* rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta
  )




source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
source("Bidart/R/visualisation.R")
par(mfrow = c(1, 1))

plot_ellipses(locs, cbind(1, X_range) %*% beta_range, shrink = .02)

plot_pointillist_painting(locs = locs, field = as.vector(cbind(1, X_range)%*% beta_range[,1]), cex = 2)
plot_pointillist_painting(locs = locs, field = as.vector(cbind(1, X_range)%*% beta_range[,2]), cex = 2)

plot_pointillist_painting(locs = locs, field = latent_field, cex = 1)

plot_pointillist_painting(locs = observed_locs, field = observed_field)



source("Bidart/R/mcmc_nngp_initialize_nonstationary.R")
source("Bidart/R/Coloring.R")
mcmc_nngp_list = mcmc_nngp_initialize_nonstationary (
  observed_locs = observed_locs, #spatial locations
  observed_field = c(observed_field), # Response variable
  X = X, # Covariates per observation
  m = 5, #number of Nearest Neighbors
  reordering = c("maxmin"), #Reordering
  covfun = "nonstationary_exponential_anisotropic", response_model = "Gaussian", # covariance model and response model
  noise_X = X_noise, noise_range = NULL, # range for latent field of parameters, if NULL no latent field
  scale_X = NULL, scale_range = NULL, # range for latent field of parameters, if NULL no latent field
  range_X = as.data.frame(as.matrix(X_range)[observation_idx,]), range_range = NULL, # range for latent field of parameters, if NULL no latent field
  log_NNGP_matern_covfun = "matern_isotropic", # covariance function for the hyperpriors
  log_NNGP_matern_smoothness = 1, # covariance function for the hyperpriors
  n_chains = 3,  # number of MCMC chains
  seed = 4
)

source("Bidart/R/mcmc_nngp_update_nonstationary_Gaussian.R")
source("Bidart/R/mcmc_nngp_run.R")
mcmc_nngp_list = mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_cores = 3, n_iterations_update = 100, n_cycles = 1)



mcmc_nngp_list$records$chain_1$range_beta

