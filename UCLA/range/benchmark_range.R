orig_par <- par()
benchmark = list()
seed = 40

###################
# Simulating data #
###################

# simulatig data
set.seed(seed)
# spatial locations
locs = cbind(runif(11000), runif(11000))
# nonstat log range
log_range = (GpGp::fast_Gp_sim(c(1,.2, 1.5, 0), locs = rbind(locs), m=10) -  8)*.5
# latent field
NNarray = GpGp::find_ordered_nn(locs, 10) 
w = 
  GpGp::fast_Gp_sim_Linv(
  Bidart::nonstat_vecchia_Linv(
    log_range = 2*as.matrix(log_range), 
    covfun_name = "nonstationary_matern_isotropic", 
    sphere = F, locs = locs, NNarray = NNarray, 
    num_threads = 10, compute_derivative = F
  )[[1]], NNarray = NNarray, z = rnorm(nrow(locs))
)
# noise variance
tau = exp(0)
# adding noise to get observations
z = (w[seq(nrow(locs))] + sqrt(tau) * rnorm(nrow(locs)))
# splitting train test
predicted_locs = locs[-seq(10000),]
predicted_z = z[-seq(10000)]
locs = locs[seq(10000),]
z = z[seq(10000)]
log_range = log_range[seq(10000)]

# visualizing data
plot_pointillist = function(locs, x, main)plot(locs, xlab = "", ylab = "", main = main, pch  =15, col = heat.colors(n = 100)[as.numeric(cut(x, breaks = seq(min(x), max(x), length.out = 98)))], axes = F)
par(mfrow = c(1, 3))
par(mar = c(1,1,1,1))
plot_pointillist(locs, w, "noiseless field")
plot_pointillist(locs, log_range, "log range")
plot_pointillist(locs, z, "observed noisy field")


############
# Local GP #
############
library(hetGP); library(lhs)
library(liGP); library(laGP)

t1 = Sys.time()
# local GP design
lhs_design <- randomLHS(10,2)
n <- 80
Xmt <- scale_ipTemplate(locs, n, space_fill_design=lhs_design, method='qnorm')$Xm.t
theta_prior <- darg(NULL, locs)
liGP_r_smooth <- liGP(XX=locs, X=locs, Y=z, Xm=Xmt, N=n, theta = theta_prior,
             g = tau, epsK=1e-5, num_thread = 10)
liGP_r_pred <- liGP(XX=predicted_locs, X=locs, Y=z, Xm=Xmt, N=n, theta = theta_prior,
             g = tau, epsK=1e-5, num_thread = 10)

# plot
par(mar = c(2,2,1,1))
par(mfrow=  c(2, 1))
plot_pointillist(locs, log_range, "true log range")
plot_pointillist(locs, log(liGP_r_smooth$mle$theta), "estimated lengthscale parameter")
cor(log(liGP_r_smooth$mle$theta), log_range)

benchmark$local_GP = c(
  w_smooth_MSE = mean((liGP_r_smooth$mean - w[seq(nrow(locs))])^2),
  w_pred_MSE =  mean((liGP_r_pred$mean - w[-seq(nrow(locs))])^2), 
  log_range_MSE =  mean((log(liGP_r_smooth$mle$theta) -log(2) - 0.5 * (log(1.5))  -  log_range)^2), 
  time = Sys.time()-t1, 
  seed = seed
)
mean((log(liGP_r_smooth$mle$theta) -log(2) - 0.5 * (log(1.5))  -  log_range))



################
# Nonstat NNGP #
################
par(orig_par)
t1 = Sys.time()
PP = Bidart::get_PP(observed_locs = locs, matern_range = .2, n_PP = 30)
Bidart::compare_PP_NNGP(PP)
mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
  observed_locs = locs, observed_field = z, 
  range_PP = T, PP = PP
)
for(i in seq(12))mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(mcmc_nngp_list, num_threads_per_chain = 10)
estimation = Bidart::estimate_parameters(
  mcmc_nngp_list, 
  get_samples = T, burn_in = .2)
fixed_effects = 
  Bidart::predict_fixed_effects(
    mcmc_nngp_list, 
    burn_in = .2)
fixed_effects = outer(
  rep(1, nrow(locs) + nrow(predicted_locs)), 
  t(matrix(fixed_effects$predicted_samples))
)
NNGP_pred = Bidart::predict_latent_field(
  predicted_locs = predicted_locs, 
  mcmc_nngp_list = mcmc_nngp_list, 
  num_threads_per_chain = 10, burn_in = .2)

estimated_log_range_NNGP = 
  Bidart::X_PP_mult_right(
    X = matrix(1, nrow(locs)), 
    PP = PP, use_PP = T, 
    Y = estimation$summaries$range_beta[1,,]) 
par(mfrow = c(1,2))
plot_pointillist(locs, estimated_log_range_NNGP, main = "estimated log range")
plot_pointillist(locs, log_range, main = "true log range")

benchmark$nonstat_NNGP = c(
  w_smooth_MSE = mean((estimation$summaries$beta[1,1,1] + estimation$summaries$field_at_observed_locs[1,,] - w[seq(nrow(locs))])^2),
  w_pred_MSE =   mean((estimation$summaries$beta[1,1,1] + NNGP_pred$summaries$field[1,,] - w[-seq(nrow(locs))])^2), 
  log_range_MSE = mean((log_range- estimated_log_range_NNGP)^2),
  time = Sys.time()-t1, 
  seed = seed
)

########
# INLA #
########

t1 = Sys.time()
library(INLA)
# mesh
prmesh = INLA::inla.mesh.2d(locs, cutoff = .015)
plot(prmesh)
# spatial basis functions
PP_INLA = Bidart::get_PP(observed_locs = prmesh$loc[,-3], matern_range = .2, n_PP = 30)
Bidart::compare_PP_NNGP(PP_INLA)
PP_basis = Matrix::solve(
  PP_INLA$sparse_chol, 
  Matrix::sparseMatrix(
    i = seq(PP_INLA$n_PP),
    j = seq(PP_INLA$n_PP),
    x = rep(1, PP_INLA$n_PP),
    dims = c(nrow(PP_INLA$sparse_chol), PP_INLA$n_PP)
  )
)[-seq(PP_INLA$n_PP),][(PP_INLA$idx),]
plot_pointillist(prmesh$loc[,-3], PP_basis[,1], "PP1")
plot_pointillist(prmesh$loc[,-3], PP_basis[,2], "PP2")

# Matern smoothness
nu <- 1 
alpha = nu + 1
# SPDE model
B.kappa = as.matrix(cbind(0, 0, 1, PP_basis))
B.tau =   as.matrix(cbind(0, 1, 0, -nu*PP_basis) )
spde <- inla.spde2.matern(
  prmesh, 
  B.tau = B.tau, 
  B.kappa = B.kappa, 
  theta.prior.mean = rep(0, ncol(B.kappa) - 1), 
  theta.prior.prec = rep(1, ncol(B.kappa) - 1), 
  alpha = nu+1
)

# stack
projloc = inla.mesh.projector(prmesh, locs)
stk <- inla.stack(
  data = list(y = z),
  A = list(projloc$proj$A),
  effects = list(data.frame(i = 1:prmesh$n))
)
# formula without intercept
formula <- z ~ 0 + f(i, model = spde)
# running
inla_nonstat <- INLA::inla(formula, data = inla.stack.data(stk), 
                     family = "gaussian", 
                     control.predictor = 
                       list(
                         compute = TRUE, 
                         A = inla.stack.A(stk))
)

# visualizing field interpolation.
# interpolation seems all good !
origin_params = par()
par(mfrow = c(1, 2))
par(mar = c(1,1,1,1))
plot_pointillist(locs, w, "noiseless field")
plot_pointillist(locs, as.vector(projloc$proj$A %*% inla_nonstat$summary.random$i$mean), "mean SPDE effects")
par(origin_params)
plot(w[seq(nrow(locs))],xlab =  "noiseless field", as.vector(projloc$proj$A %*% inla_nonstat$summary.random$i$mean), ylab =  "mean SPDE effects")
abline(a =0, b=1)

# estimating nonstat Matérn range and marginal variance 
# starting with log tau and log kappa
log_tau_estimated_nonstat   = B.tau[,-1]   %*% inla_nonstat$summary.hyperpar$mean[-1]
log_kappa_estimated_nonstat = B.kappa[,-1] %*% inla_nonstat$summary.hyperpar$mean[-1]

# plotting 
par(mfrow = c(1,2))
plot_pointillist(locs, -as.vector(projloc$proj$A %*% log_kappa_estimated_nonstat), main =  "- log kappa")
plot_pointillist(locs, log_range, "log range")

## Marginal variance is σ2 =Γ(ν)/(Γ(α) (4π)^(d/2) κ^(2ν) τ^2 )
## check that it is constant
# estimated_matern_marginal_var = 
#   exp(
#     + lgamma(nu) - lgamma(nu + 1) - log(4*pi)
#     - 2 * nu * B.kappa[,-1] %*% inla_nonstat$summary.hyperpar$mean[-1]
#     - 2 *      B.tau[,-1]   %*% inla_nonstat$summary.hyperpar$mean[-1]
#   )

# converting log-range because of change in smoothness
# using rho NNGP = rho INLA 

log_range_INLA = .5 * (log(nu) - 1.5) - log_kappa_estimated_nonstat
plot_pointillist(locs, as.vector(projloc$proj$A %*% log_range_INLA), main =  "INLA range")
plot_pointillist(locs, log_range, "log range")
par(mfrow = c(1,1))
plot(projloc$proj$A %*%log_range_INLA, log_range) 
abline(a=0, b=1)

benchmark$INLA = c(
  w_smooth_MSE = mean((as.vector(projloc$proj$A %*% inla_nonstat$summary.random$i$mean) - w[seq(nrow(locs))])^2),
  w_pred_MSE = mean((as.vector(inla.mesh.projector(prmesh, predicted_locs)$proj$A %*% inla_nonstat$summary.random$i$mean) - w[-seq(nrow(locs))])^2),
  log_range_MSE = mean((log_range- -as.vector(projloc$proj$A %*% log_kappa_estimated_nonstat))^2),
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
  log_range_MSE = NA,
  time = Sys.time()-t1, 
  seed = seed
)
