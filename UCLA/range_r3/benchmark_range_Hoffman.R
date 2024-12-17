orig_par <- par()
benchmark = list()

.libPaths("../R_packages/")


library(spam)
library(FNN)
library(fields)
library(GpGp)
library(Bidart)
library(laGP)
library(liGP)
library(hetGP)
library(ellipse)
library(coda)
library(lhs)
library(snow)
library(doParallel)
library(doSNOW)
library(sp)

seed = 0
seed = seq(51)[which(! (paste("range", seq(51), ".RDS", sep = "") %in% list.files()))][1]

saveRDS("tatato", paste("range", seed, ".RDS", sep = ""))
###################
# Simulating data #
###################

# simulatig data
set.seed(seed)
# spatial locations
locs = cbind(runif(11000), runif(11000))
# nonstat log range
log_range = (GpGp::fast_Gp_sim(c(1,.2, 1.5, 0), locs = rbind(locs), m=10) -  3.7)
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
tau = exp(-1)
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


########
# INLA #
########

library(INLA)
t1 = Sys.time()
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
projloc_pred = inla.mesh.projector(prmesh, predicted_locs)
#obs
stk = 
  inla.stack(
    data = list(z = z),
    A = list(INLA::inla.spde.make.A(
      prmesh, # SPDE mesh
      loc=locs # coordinates of the observations of the interest variable
    )),
    effects = list(data.frame(i = 1:prmesh$n))
  )

#pred
stk_pred = 
  inla.stack(
    data = list(z = NA),
    A = list(INLA::inla.spde.make.A(
      prmesh, # SPDE mesh
      loc=predicted_locs # coordinates of the observations of the interest variable
    )),
    effects = list(data.frame(i = 1:prmesh$n))
  )

stk = inla.stack(stk, stk_pred)
# formula without intercept
formula <- z ~ 0 + f(i, model = spde)
# running
inla_nonstat <- inla(formula, data = inla.stack.data(stk), 
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
plot_pointillist(locs, as.vector(projloc$proj$A %*% inla_nonstat$summary.random$i$mean), "mean SPDE effects using A")
plot_pointillist(locs, as.vector(inla_nonstat$summary.fitted.values$mean[seq(10000)]), "mean SPDE effects using summary")
plot(as.vector(projloc$proj$A %*% inla_nonstat$summary.random$i$mean), xlab = "mean SPDE effects using A",
     as.vector(inla_nonstat$summary.fitted.values$mean[seq(10000)]), ylab = "mean SPDE effects using summary"
)

#checking $fitted_values and $summary.random are the same (at mesh nodes)
match_mesh_locs = match(
  split(prmesh$loc[,-3], row(prmesh$loc[,-3])),
  split(locs, row(locs))
)
plot(
  inla_nonstat$summary.random$i$sd[!is.na(match_mesh_locs)],
inla_nonstat$summary.fitted.values$sd[
  match_mesh_locs[!is.na(match_mesh_locs)]]
)

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

log_range_INLA = -.5 * log(1.5) - log_kappa_estimated_nonstat
plot_pointillist(locs, as.vector(projloc$proj$A %*% log_range_INLA), main =  "INLA range")
plot_pointillist(locs, log_range, "log range")
par(mfrow = c(1,1))
plot(projloc$proj$A %*% log_range_INLA, log_range) 
abline(a=0, b=1)


# elpd and coverage
INLA_noise_var = 1/inla_nonstat$summary.hyperpar[1,"mode"]
# observed
plot(inla_nonstat$summary.fitted.values$mean[seq(10000)], w[seq(10000)])
INLA_var_obs = inla_nonstat$summary.fitted.values$sd[seq(10000)]^2 + INLA_noise_var
INLA_mean_obs = inla_nonstat$summary.fitted.values$mean[seq(10000)]
obs_elpd = mean(dnorm(z, INLA_mean_obs, sqrt(INLA_var_obs), log=T))
obs_cover = mean((z > INLA_mean_obs - 1.96*sqrt(INLA_var_obs)) &(z < INLA_mean_obs + 1.96*sqrt(INLA_var_obs)))

plot(inla_nonstat$summary.fitted.values$mean[seq(10001, 11000)], w[-seq(10000)])
INLA_var_pred =  inla_nonstat$summary.fitted.values$sd[seq(10001, 11000)]^2 + INLA_noise_var
INLA_mean_pred = inla_nonstat$summary.fitted.values$mean[seq(10001, 11000)]
pred_elpd = mean(dnorm(predicted_z, INLA_mean_pred, sqrt(INLA_var_pred), log=T))
pred_cover = mean((predicted_z > INLA_mean_pred - 1.96*sqrt(INLA_var_pred)) &(predicted_z < INLA_mean_pred + 1.96*sqrt(INLA_var_pred)))


benchmark$INLA = c(
  w_smooth_MSE = mean((as.vector(projloc$proj$A %*% inla_nonstat$summary.random$i$mean) - w[seq(nrow(locs))])^2),
  w_pred_MSE = mean((as.vector(inla.mesh.projector(prmesh, predicted_locs)$proj$A %*% inla_nonstat$summary.random$i$mean) - w[-seq(nrow(locs))])^2),
  log_range_MSE = mean((log_range- as.vector(projloc$proj$A %*%log_range_INLA))^2),
  obs_coverage = obs_cover, 
  pred_coverage = pred_cover,
  obs_elpd = obs_elpd, 
  pred_elpd = pred_elpd,
  time = as.numeric(Sys.time()-t1, unit = "mins"), 
  seed = seed
)


print("INLA done !")


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
for(i in seq(15))mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(mcmc_nngp_list, num_threads_per_chain = 10, lib.loc = "../R_packages/")
estimation = Bidart::estimate_parameters(
  mcmc_nngp_list, 
  get_samples = T, lib.loc = "../R_packages/",
  burn_in = .2)
fixed_effects = 
  Bidart::predict_fixed_effects(
    mcmc_nngp_list, lib.loc = "../R_packages/",
    burn_in = .2)
fixed_effects = outer(
  rep(1, nrow(locs) + nrow(predicted_locs)), 
  t(matrix(fixed_effects$predicted_samples))
)
NNGP_pred = Bidart::predict_latent_field(
  predicted_locs = predicted_locs, 
  mcmc_nngp_list = mcmc_nngp_list, 
  num_threads_per_chain = 10, lib.loc = "../R_packages/",
  burn_in = .2)
noise_NNGP = Bidart::predict_noise(
  mcmc_nngp_list = mcmc_nngp_list, 
  X_noise_pred = NULL, 
  predicted_locs = rbind(locs, predicted_locs),
  burn_in = .2
)


estimated_log_range_NNGP = 
  Bidart::X_PP_mult_right(
    X = matrix(1, nrow(locs)), 
    PP = PP, use_PP = T, 
    Y = estimation$summaries$range_beta[1,,]) 
par(mfrow = c(1,2))
plot_pointillist(locs, estimated_log_range_NNGP, main = "estimated log range")
plot_pointillist(locs, log_range, main = "true log range")


# coverage
preds_mean_samples = fixed_effects[-seq(10000),,,drop=F] + NNGP_pred$predicted_samples$field
preds_noisevar_samples = exp(noise_NNGP$predicted_samples[-seq(10000),,,drop = F]) 
response_preds_samples = (preds_mean_samples + sqrt(preds_noisevar_samples) * rnorm(length(preds_noisevar_samples)))[,1,]
quantiles_preds = apply(response_preds_samples, 1, function(x)quantile(x, probs = c(.025, .975)))
pred_coverage = mean((predicted_z > quantiles_preds[1,])&(predicted_z <quantiles_preds[2,]))

obs_mean_samples = fixed_effects[seq(10000),,,drop=F] + estimation$samples$field_at_observed_locs
obs_noisevar_samples = exp(noise_NNGP$predicted_samples[seq(10000),,,drop = F]) 
response_obs_samples = (obs_mean_samples + sqrt(obs_noisevar_samples) * rnorm(length(obs_noisevar_samples)))[,1,]
quantiles_obs = apply(response_obs_samples, 1, function(x)quantile(x, probs = c(.025, .975)))
coverage = mean((z > quantiles_obs[1,])&(z <quantiles_obs[2,]))


# elpd
NNGP_obs_elpd = Bidart::log_score_Gaussian(
  observed_field = z,
  latent_field_samples = estimation$samples$field_at_observed_locs, 
  log_noise_samples = noise_NNGP$predicted_samples[seq(10000),,,drop = F],
  fixed_effects_samples =  fixed_effects[seq(10000),,,drop=F]
)$total / 10000

NNGP_pred_elpd = Bidart::log_score_Gaussian(
  observed_field = predicted_z,
  latent_field_samples = NNGP_pred$predicted_samples$field, 
  log_noise_samples = noise_NNGP$predicted_samples[-seq(10000),,,drop = F],
  fixed_effects_samples =  fixed_effects[-seq(10000),,,drop=F]
)$total / 1000

benchmark$nonstat_NNGP = c(
  w_smooth_MSE = mean((estimation$summaries$beta[1,1,1] + estimation$summaries$field_at_observed_locs[1,,] - w[seq(nrow(locs))])^2),
  w_pred_MSE =   mean((estimation$summaries$beta[1,1,1] + NNGP_pred$summaries$field[1,,] - w[-seq(nrow(locs))])^2), 
  log_range_MSE = mean((log_range- estimated_log_range_NNGP)^2),
  obs_coverage = coverage, 
  pred_coverage = pred_coverage,
  obs_elpd = NNGP_obs_elpd, 
  pred_elpd = NNGP_pred_elpd,
  time = as.numeric(Sys.time()-t1, unit = "mins"), 
  seed = seed
)

pdf(paste("nonstat_NNGP", seed, ".pdf", sep = ""), width = 28, height = 7)
Bidart::diagnostic_plots(mcmc_nngp_list, burn_in = .2, starting_proportion = .2) 
dev.off()
print("Nonstat NNGP done !")

############
# Local GP #
############

t1 = Sys.time()
# local GP design
lhs_design <- randomLHS(10,2)
n <- 80
Xmt <- scale_ipTemplate(locs, n, space_fill_design=lhs_design, method='qnorm')$Xm.t
theta_prior <- darg(NULL, locs)
liGP_r_smooth <- liGP(XX=locs, X=locs, Y=z, Xm=Xmt, N=n, theta = theta_prior,
                      g = tau, epsK=1e-5, num_thread = 12)
liGP_r_pred <- liGP(XX=predicted_locs, X=locs, Y=z, Xm=Xmt, N=n, theta = theta_prior,
                    g = tau, epsK=1e-5, num_thread = 12)

# plot
par(mar = c(2,2,1,1))
par(mfrow=  c(2, 1))
plot_pointillist(locs, log_range, "true log range")
plot_pointillist(locs, log(liGP_r_smooth$mle$theta), "estimated lengthscale parameter")
cor(log(liGP_r_smooth$mle$theta), log_range)
plot(log(liGP_r_smooth$mle$theta), log_range)

# coverage and elpd
locGP_obs_coverage = mean((z<liGP_r_smooth$mean + qnorm(.975)*sqrt(liGP_r_smooth$nu * liGP_r_smooth$var))&(z>liGP_r_smooth$mean + qnorm(.025)*sqrt(liGP_r_smooth$nu * liGP_r_smooth$var)))
locGP_pred_coverage = mean((predicted_z<liGP_r_pred$mean + qnorm(.975)*sqrt(liGP_r_pred$nu * liGP_r_pred$var))&(predicted_z>liGP_r_pred$mean + qnorm(.025)*sqrt(liGP_r_pred$nu * liGP_r_pred$var)))
locGP_obs_elpd = mean(dnorm(z, liGP_r_smooth$mean, sqrt(liGP_r_smooth$nu * liGP_r_smooth$var), log= T))
locGP_pred_elpd =  mean(dnorm(predicted_z, liGP_r_pred$mean, sqrt(liGP_r_pred$nu * liGP_r_pred$var), log= T))




benchmark$local_GP = c(
  w_smooth_MSE = mean((liGP_r_smooth$mean - w[seq(nrow(locs))])^2),
  w_pred_MSE =  mean((liGP_r_pred$mean - w[-seq(nrow(locs))])^2), 
  log_range_MSE =  mean((log(liGP_r_smooth$mle$theta) -log(2) - 0.5 * (log(1.5))  -  log_range)^2), 
  obs_coverage = locGP_obs_coverage, 
  pred_coverage = locGP_pred_coverage,
  obs_elpd = locGP_obs_elpd, 
  pred_elpd = locGP_pred_elpd, 
  time = as.numeric(Sys.time()-t1, unit = "secs"), 
  seed = seed
)
###########
# Deep GP #
###########

t1 = Sys.time()

library(deepgp)
deepGP_fit = fit_two_layer(locs, z, nmcmc = 25000, verb = T, vecchia = T)
deepGP_fit = trim(deepGP_fit, deepGP_fit$nmcmc/5, 10)
deepgp_smooth = predict(object = deepGP_fit, x_new = locs, cores = 1)
par(mfrow = c(1,2))
plot_pointillist(locs, deepgp_smooth$mean, "deepGP smooth mean")
plot_pointillist(locs, w[seq(10000)], "true w")
deepgp_pred =  predict(object = deepGP_fit, x_new = predicted_locs, cores = 1)
plot_pointillist(predicted_locs, deepgp_pred$mean, "deepGP smooth mean")
plot_pointillist(predicted_locs, w[-seq(10000)], "true w pred")

deepGP_cover_smooth = mean(z > (deepgp_smooth$mean - 1.96 * sqrt(deepgp_smooth$s2)) & z < (deepgp_smooth$mean + 1.96 * sqrt(deepgp_smooth$s2)))
deepGP_cover_pred = mean(predicted_z > (deepgp_pred$mean - 1.96 * sqrt(deepgp_pred$s2)) & predicted_z < (deepgp_pred$mean + 1.96 * sqrt(deepgp_pred$s2)))

deepGP_elpd_smooth = mean(dnorm(z, deepgp_smooth$mean, sqrt(deepgp_smooth$s2), log = T))
deepGP_elpd_pred = mean(dnorm(predicted_z, deepgp_pred$mean, sqrt(deepgp_pred$s2), log = T))



benchmark$deepGP = c(
  w_smooth_MSE = mean((deepgp_smooth$mean - w[seq(nrow(locs))])^2),
  w_pred_MSE =   mean((deepgp_pred$mean - w[-seq(nrow(locs))])^2), 
  log_range_MSE = NA,
  obs_coverage = deepGP_cover_smooth, 
  pred_coverage = deepGP_cover_pred,
  obs_elpd = deepGP_elpd_smooth, 
  pred_elpd = deepGP_elpd_pred,
  time = as.numeric(Sys.time()-t1, unit = "hours"), 
  seed = seed
)




pdf(paste("deepGP", seed, ".pdf", sep = ""), width = 28, height = 7)
plot(deepGP_fit) 
dev.off()
print("deepGP done!")

saveRDS(benchmark, paste("range", seed, ".RDS", sep = ""))
seed = seq(51)[which(! (paste("range", seq(51), ".RDS", sep = "") %in% list.files()))][1]

