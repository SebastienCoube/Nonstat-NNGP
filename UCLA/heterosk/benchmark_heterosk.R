orig_par <- par()
benchmark = list()
seed = 1

###################
# Simulating data #
###################

set.seed(seed)
# spatial locations
locs = cbind(runif(11000), runif(11000))
# noiseless GP 
w_range = .025
w = GpGp::fast_Gp_sim(c(1,w_range, 1.5, 0), locs = rbind(locs), m=10)
# heteroscedastic noise variance
log_tau_var = 1
log_tau_range = .2
tau = exp(GpGp::fast_Gp_sim(c(log_tau_var,.2, 1.5, 0), locs = rbind(locs), m=10) - 1)
# adding noise to get observations
z = (w[seq(nrow(locs))] + sqrt(tau) * rnorm(nrow(locs)))
# splitting train test
predicted_locs = locs[-seq(10000),]
predicted_z = z[-seq(10000)]
locs = locs[seq(10000),]
z = z[seq(10000)]
tau = tau[seq(10000)]

# visualizing data
plot_pointillist = function(locs, x, main)plot(locs, xlab = "", ylab = "", main = main, pch  =15, col = heat.colors(n = 100)[as.numeric(cut(x, breaks = seq(min(x), max(x), length.out = 98)))], axes = F)
par(mfrow = c(1, 3))
par(mar = c(1,1,1,1))
plot_pointillist(locs, w, "noiseless field")
plot_pointillist(locs, log(tau), "log variance of the noise")
plot_pointillist(locs, z, "observed noisy field")

############
# Local GP #
############
## library(hetGP); library(lhs)
## library(liGP); library(laGP)
## t1 = Sys.time()
## lhs_design <- randomLHS(10,2)
## n <- 80
## Xmt <- scale_ipTemplate(locs, n, space_fill_design=lhs_design, method='qnorm')$Xm.t
## g_prior <- garg(list(mle=TRUE), z)
## liGP_n_smooth <- liGP(XX=locs, X=locs, Y=z, Xm=Xmt, N=n, theta = w_range,
##                       g = g_prior, epsK=1e-5, num_thread = 10)
## liGP_n_pred <- liGP(XX=predicted_locs, X=locs, Y=z, Xm=Xmt, N=n, theta = w_range,
##                     g = g_prior, epsK=1e-5, num_thread = 10)
## log(liGP_n_smooth$nu * liGP_n_smooth$var) - log(tau)
## # saving result
## benchmark$local_GP = c(
##   w_smooth_MSE = mean((liGP_n_smooth$mean - w[seq(nrow(locs))])^2),
##   w_pred_MSE =  mean((liGP_n_pred$mean - w[-seq(nrow(locs))])^2), 
##   log_tau_MSE =  mean((log(liGP_n_smooth$nu * liGP_n_smooth$var) -  log(tau))^2), 
##   time = Sys.time()-t1, 
##   seed = seed
## )
## # mean(dnorm(z, liGP_n_smooth$mean, sqrt(liGP_n_smooth$var), log=T))
## # mean(dnorm(predicted_z, liGP_n_pred$mean, sqrt(liGP_n_pred$var), log=T))
## ## Nice retrieval of noise variance 
## # par(mfrow=  c(2, 1))
## # par(mar = c(2,2,1,1))
## # plot_pointillist(locs, log(liGP_n_smooth$nu * liGP_n_smooth$var), "estimated log var liGP")
## # plot_pointillist(locs, log(tau), "true log var")
## # par(mfrow=  c(1, 1))
## # par(mar = c(4,4,1,1))
## # plot(log(liGP_n_smooth$nu * liGP_n_smooth$var), log(tau), 
## #      xlab = "estimated log nugget", ylab = "true log nugget")
## # abline(a = 0, b=1)




################
# Nonstat NNGP #
################
par(orig_par)
t1 = Sys.time()
PP = Bidart::get_PP(observed_locs = locs, matern_range = .2, n_PP = 30)
Bidart::compare_PP_NNGP(PP)
mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
  observed_locs = locs, observed_field = z, 
  noise_PP = T, PP = PP
)
for(i in seq(15))mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(mcmc_nngp_list, num_threads_per_chain = 10)
estimation = Bidart::estimate_parameters(mcmc_nngp_list, get_samples = T)
fixed_effects = Bidart::predict_fixed_effects(mcmc_nngp_list)
fixed_effects = outer(rep(1, nrow(locs) + nrow(predicted_locs)), t(matrix(fixed_effects$predicted_samples)))
noise_NNGP = Bidart::predict_noise(
  mcmc_nngp_list = mcmc_nngp_list, 
  X_noise_pred = NULL, 
  predicted_locs = rbind(locs, predicted_locs))
NNGP_pred = Bidart::predict_latent_field(
  predicted_locs = predicted_locs, mcmc_nngp_list = mcmc_nngp_list, 
  num_threads_per_chain = 10)

benchmark$nonstat_NNGP = c(
  w_smooth_MSE = mean((estimation$summaries$beta[1,1,1] + estimation$summaries$field_at_observed_locs[1,,] - w[seq(nrow(locs))])^2),
  w_pred_MSE =   mean((estimation$summaries$beta[1,1,1] + NNGP_pred$summaries$field[1,,] - w[-seq(nrow(locs))])^2), 
  log_tau_MSE =  mean((noise_NNGP$summaries[1,seq(nrow(locs)),1] -  log(tau))^2), 
  time = Sys.time()-t1, 
  seed = seed
)


### log-scores
#log_score_smooth = Bidart::log_score_Gaussian(
#  observed_field = z, 
#  latent_field_samples = estimation$samples$field_at_observed_locs, 
#  log_noise_samples = noise_NNGP$predicted_samples[seq(nrow(locs)),,,drop=F], 
#  fixed_effects_samples = fixed_effects[seq(nrow(locs)),,,drop=F])
#log_score_pred = Bidart::log_score_Gaussian(
#  observed_field = predicted_z, 
#  latent_field_samples = NNGP_pred$predicted_samples$field, 
#  log_noise_samples = noise_NNGP$predicted_samples[-seq(nrow(locs)),,,drop=F], 
#  fixed_effects_samples = fixed_effects[-seq(nrow(locs)),,,drop=F])
#(log_score_smooth$total)/nrow(locs)
#(log_score_pred$total)/nrow(predicted_locs)

########
# INLA #
########

t1 = Sys.time()
# making mesh
prmesh = INLA::inla.mesh.2d(locs, cutoff = .015)
prmesh$n
plot(prmesh)
# making noise
PP = Bidart::get_PP(observed_locs = locs, matern_range = .2, n_PP = 30)
Bidart::compare_PP_NNGP(PP)
X_noise = cbind(1, as.matrix(Matrix::solve(PP$sparse_chol, diag(1, nrow(PP$sparse_chol), PP$n_PP))[-seq(PP$n_PP),][PP$idx,]))
plot_pointillist(locs, X_noise[,2], "PP 1")

# X_noise corresponds to observations of the covariates
# ***at the coordinates of the observations***
A_noise = Matrix::Diagonal(n=nrow(X_noise), x=1)
# A noise is an identity matrix with
# as many rows and columns as the number of observations
spde_noise <- INLA::inla.spde2.generic(
  M0 = Matrix::Diagonal(n=nrow(X_noise), x=1.0),
  M1 = Matrix::Diagonal(n=nrow(X_noise), x=0.0),
  M2 = Matrix::Diagonal(n=nrow(X_noise), x=0.0),
  B0 = cbind(0, X_noise), B1 = matrix(0, 1, 1+ ncol(X_noise)),
  B2 = matrix(0, 1, 1+ ncol(X_noise)),
  theta.mu = rep(0,ncol(X_noise)),
  theta.Q = diag(rep(1/log_tau_var,ncol(X_noise))), # true prior
  transform = "identity")  


# latent effects
A_field = 
  INLA::inla.spde.make.A(
    prmesh, # SPDE mesh
    loc=locs # coordinates of the observations of the interest variable
  )
spde_spatially_coherent = 
  INLA::inla.spde2.matern(mesh = prmesh, alpha=2)

# stacking

stk <- INLA::inla.stack(
  data = list(y = z),
  A = list(A_noise,
           A_field
  ),
  effects=list(
    INLA::inla.spde.make.index(
      name = "heteroscedastic.noise",
      spde_noise$n.spde), #the noise index
    # the length of the noise index is equal to the number of observations
    INLA::inla.spde.make.index(
      name = "spatial.field",
      spde_spatially_coherent$n.spde)  #the spatial field index
    # the length of the spatial index is equal to the number of mesh nodes
  ))

formula = z ~
  f(heteroscedastic.noise, model = spde_noise)+ # heteroskedastic Gaussian noise
  f(spatial.field, model = spde_spatially_coherent) # spatially coherent random effect

clik <- list(hyper = list(prec = list(initial = 20,
                                      fixed = TRUE)))
res <- INLA::inla(
  formula, data = INLA::inla.stack.data(stk), control.family = clik,
  control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(stk)), 
  verbose = T, control.compute = list(cpo = TRUE, dic = TRUE)
)


#Bidart::plot_pointillist_painting(prmesh$loc, res$summary.random$spatial.field$mean)
Bidart::plot_pointillist_painting(locs, w)
Bidart::plot_pointillist_painting(
  locs, 
  as.vector(A_field %*% res$summary.random$spatial.field$mean) + res$summary.fixed$mean
) 
Bidart::plot_pointillist_painting(locs, -2*spde_noise$param.inla$B0[,-1] %*% res$summary.hyperpar$mode[seq(PP$n_PP + 1)])
Bidart::plot_pointillist_painting(locs, log(tau))
plot(-2*spde_noise$param.inla$B0[,-1] %*% res$summary.hyperpar$mode[seq(PP$n_PP + 1)], log(tau))
abline(a=0, b=1)

benchmark$INLA = c(
  w_smooth_MSE = mean((as.vector(A_field %*% res$summary.random$spatial.field$mean) + res$summary.fixed$mean - w[seq(nrow(locs))])^2),
  w_pred_MSE = mean((as.vector(INLA::inla.spde.make.A(prmesh, loc=predicted_locs) %*% res$summary.random$spatial.field$mean) + res$summary.fixed$mean - w[-seq(nrow(locs))])^2),
  log_tau_MSE =  mean((log(tau) + 2*spde_noise$param.inla$B0[,-1] %*% res$summary.hyperpar$mode[seq(PP$n_PP + 1)])^2), 
  time = Sys.time()-t1, 
  seed = seed
)



###########
# Deep GP #
###########

t1 = Sys.time()

library(deepgp)
deepGP_fit = fit_two_layer(locs, z, nmcmc = 10000, verb = T, vecchia = T)
deepGP_fit = trim(deepGP_fit, deepGP_fit$nmcmc/5, 10)
deepgp_smooth = predict(object = deepGP_fit, x_new = locs, cores = 10)
par(mfrow = c(1,2))
plot_pointillist(locs, deepgp_smooth$mean, "deepGP smooth mean")
plot_pointillist(locs, w[seq(10000)], "true w")
deepgp_pred =  predict(object = deepGP_fit, x_new = predicted_locs, cores = 10)
plot_pointillist(predicted_locs, deepgp_pred$mean, "deepGP smooth mean")
plot_pointillist(predicted_locs, w[-seq(10000)], "true w pred")
benchmark$deepGP = c(
  w_smooth_MSE = mean((deepgp_smooth$mean - w[seq(nrow(locs))])^2),
  w_pred_MSE =   mean((deepgp_pred$mean - w[-seq(nrow(locs))])^2), 
  log_tau_MSE = NA,
  time = Sys.time()-t1, 
  seed = seed
)
