orig_par <- par()
benchmark = list()

.libPaths("../R_packages/")


library(spam)
library(FNN)
library(fields)
library(hetGP)
library(ellipse)
library(coda)
library(GpGp)
library(Bidart)
library(laGP)
library(liGP)
library(lhs)
library(snow)
library(doParallel)
library(sp)

seed = 0
seed = seq(51)[which(! (paste("het", seq(51), ".RDS", sep = "") %in% list.files()))][1]

saveRDS("tatato", paste("het", seed, ".RDS", sep = ""))
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

########
# INLA #
########

library(INLA)
#inla.binary.install()
#library(INLA)
t1 = Sys.time()
# making mesh
prmesh = INLA::inla.mesh.2d(rbind(locs, predicted_locs)#, cutoff = .015
)
loc_idx_in_mesh = match(split(locs, row(locs)), split(prmesh$loc[,c(1,2)], row(prmesh$loc[,c(1,2)])))
predicted_loc_idx_in_mesh = match(split(predicted_locs, row(predicted_locs)), split(prmesh$loc[,c(1,2)], row(prmesh$loc[,c(1,2)])))
prmesh$n
plot(prmesh)
# making noise
PP = Bidart::get_PP(observed_locs = rbind(locs, predicted_locs), matern_range = .2, n_PP = 30)
Bidart::compare_PP_NNGP(PP)
X_noise = cbind(1, as.matrix(Matrix::solve(PP$sparse_chol, diag(1, nrow(PP$sparse_chol), PP$n_PP))[-seq(PP$n_PP),][PP$idx,]))
pred_X_noise = X_noise[-seq(10000),]
X_noise = X_noise[seq(10000),]

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

formula = z ~ -1 + 
  f(heteroscedastic.noise, model = spde_noise)+ # heteroskedastic Gaussian noise
  f(spatial.field, model = spde_spatially_coherent) # spatially coherent random effect

clik <- list(hyper = list(prec = list(initial = 20,
                                      fixed = TRUE)))
res <- INLA::inla(
  formula, data = INLA::inla.stack.data(stk), control.family = clik,
  control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(stk)), 
  verbose = T, control.compute = list(cpo = TRUE, dic = TRUE)
)


# plotting stuff just to be sure
#Bidart::plot_pointillist_painting(prmesh$loc, res$summary.random$spatial.field$mean)
# Bidart::plot_pointillist_painting(locs, w)
# Bidart::plot_pointillist_painting(
#   locs, 
#   as.vector(A_field %*% res$summary.random$spatial.field$mean)
# ) 
# Bidart::plot_pointillist_painting(locs, -2*spde_noise$param.inla$B0[,-1] %*% res$summary.hyperpar$mode[seq(PP$n_PP + 1)])
# Bidart::plot_pointillist_painting(locs, log(tau))
# plot(-2*spde_noise$param.inla$B0[,-1] %*% res$summary.hyperpar$mode[seq(PP$n_PP + 1)], log(tau))
# abline(a=0, b=1)

# coverage
GP_mean = res$summary.random$spatial.field$mean[loc_idx_in_mesh]
Bidart::plot_pointillist_painting(locs, GP_mean)
GP_var = res$summary.random$spatial.field$sd[loc_idx_in_mesh]^2
Bidart::plot_pointillist_painting(locs, GP_var)
noise_var = exp(-2*spde_noise$param.inla$B0[,-1] %*% res$summary.hyperpar$mode[seq(PP$n_PP + 1)])
Bidart::plot_pointillist_painting(locs, noise_var)
coverage = 
  mean(
    (z > GP_mean + sqrt(GP_var + noise_var)*qnorm(.025))&
      (z < GP_mean - sqrt(GP_var + noise_var)*qnorm(.025))
  )

pred_GP_mean = res$summary.random$spatial.field$mean[predicted_loc_idx_in_mesh]
pred_GP_var = res$summary.random$spatial.field$sd[predicted_loc_idx_in_mesh]^2
pred_noise_var = exp(- 2*pred_X_noise %*% res$summary.hyperpar$mode[seq(PP$n_PP + 1)])
pred_coverage = 
  mean(
    (predicted_z > pred_GP_mean + sqrt(pred_GP_var + pred_noise_var)*qnorm(.025))&
      (predicted_z < pred_GP_mean - sqrt(pred_GP_var + pred_noise_var)*qnorm(.025))
  )

# elpd using approximation
INLA_pred_elpd =  mean(dnorm(predicted_z, pred_GP_mean, sqrt(pred_GP_var + pred_noise_var), log = T))
INLA_train_elpd = mean(dnorm(z,           GP_mean,      sqrt(GP_var      + noise_var), log = T))



benchmark$INLA = c(
  w_smooth_MSE = mean((GP_mean  - w[seq(nrow(locs))])^2),
  w_pred_MSE = mean((pred_GP_mean - w[-seq(nrow(locs))])^2),
  log_tau_MSE =  mean((log(tau) + 2*spde_noise$param.inla$B0[,-1] %*% res$summary.hyperpar$mode[seq(PP$n_PP + 1)])^2), 
  obs_coverage = coverage, 
  pred_coverage = pred_coverage,
  obs_elpd = INLA_train_elpd, 
  pred_elpd = INLA_pred_elpd,
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
  noise_PP = T, PP = PP
)


#mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_nonpar(mcmc_nngp_list, num_threads_per_chain = 10)
for(i in seq(25))mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(mcmc_nngp_list, num_threads_per_chain = 6, lib.loc = "../R_packages/")
estimation = Bidart::estimate_parameters(mcmc_nngp_list, get_samples = T, burn_in = .2)
fixed_effects = Bidart::predict_fixed_effects(mcmc_nngp_list, lib.loc = "../R_packages/", burn_in = .2)
fixed_effects = outer(
  rep(1, nrow(locs) + nrow(predicted_locs)), 
  t(matrix(fixed_effects$predicted_samples))
)
noise_NNGP = Bidart::predict_noise(
  mcmc_nngp_list = mcmc_nngp_list, 
  X_noise_pred = NULL, 
  predicted_locs = rbind(locs, predicted_locs),
  burn_in = .2
)
NNGP_pred = Bidart::predict_latent_field(
  predicted_locs = predicted_locs, mcmc_nngp_list = mcmc_nngp_list, 
  num_threads_per_chain = 5, lib.loc = "../R_packages/",
  burn_in = .2)


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
  log_tau_MSE =  mean((noise_NNGP$summaries[1,seq(nrow(locs)),1] -  log(tau))^2), 
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
print("Nonstat NNGP done !")


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
  obs_coverage = deepGP_cover_smooth, 
  pred_coverage = deepGP_cover_pred,
  obs_elpd = deepGP_elpd_smooth, 
  pred_elpd = deepGP_elpd_pred,
  log_tau_MSE = NA,
  time = as.numeric(Sys.time()-t1, unit = "hours") , 
  seed = seed
)
pdf(paste("deepGP", seed, ".pdf", sep = ""), width = 28, height = 7)
plot(deepGP_fit) 
dev.off()
print("deepGP done!")

############
# Local GP #
############

t1 = Sys.time()
lhs_design <- randomLHS(10,2)
n <- 80
Xmt <- scale_ipTemplate(locs, n, space_fill_design=lhs_design, method='qnorm')$Xm.t
g_prior <- garg(list(mle=TRUE), z)
liGP_n_smooth <- liGP(XX=locs, X=locs, Y=z, Xm=Xmt, N=n, theta = exp(log(w_range) + log(2) + 0.5*log(1.5)),
                      g = g_prior, epsK=1e-5, num_thread = 12)
liGP_n_pred <- liGP(XX=predicted_locs, X=locs, Y=z, Xm=Xmt, N=n, theta = exp(log(w_range) + log(2) + 0.5*log(1.5)),
                    g = g_prior, epsK=1e-5, num_thread = 12)
#log(liGP_n_smooth$nu * liGP_n_smooth$var) - log(tau)

# coverage and elpd
locGP_obs_coverage = mean((z<liGP_n_smooth$mean + qnorm(.975)*sqrt(liGP_n_smooth$nu * liGP_n_smooth$var))&(z>liGP_n_smooth$mean + qnorm(.025)*sqrt(liGP_n_smooth$nu * liGP_n_smooth$var)))
locGP_pred_coverage = mean((predicted_z<liGP_n_pred$mean + qnorm(.975)*sqrt(liGP_n_pred$nu * liGP_n_pred$var))&(predicted_z>liGP_n_pred$mean + qnorm(.025)*sqrt(liGP_n_pred$nu * liGP_n_pred$var)))
locGP_obs_elpd = mean(dnorm(z, liGP_n_smooth$mean, sqrt(liGP_n_smooth$nu * liGP_n_smooth$var), log= T))
locGP_pred_elpd =  mean(dnorm(predicted_z, liGP_n_pred$mean, sqrt(liGP_n_pred$nu * liGP_n_pred$var), log= T))


# saving result
benchmark$local_GP = c(
  w_smooth_MSE = mean((liGP_n_smooth$mean - w[seq(nrow(locs))])^2),
  w_pred_MSE =  mean((liGP_n_pred$mean - w[-seq(nrow(locs))])^2), 
  log_tau_MSE =  mean((log(liGP_n_smooth$nu * liGP_n_smooth$mle$g) -  log(tau))^2), 
  obs_coverage = locGP_obs_coverage, 
  pred_coverage = locGP_pred_coverage,
  obs_elpd = locGP_obs_elpd, 
  pred_elpd = locGP_pred_elpd,
  time = as.numeric(Sys.time()-t1, unit = "secs"), 
  seed = seed
)
# mean(dnorm(z, liGP_n_smooth$mean, sqrt(liGP_n_smooth$var), log=T))
# mean(dnorm(predicted_z, liGP_n_pred$mean, sqrt(liGP_n_pred$var), log=T))
## Nice retrieval of noise variance 
# par(mfrow=  c(2, 1))
# par(mar = c(2,2,1,1))
# plot_pointillist(locs, log(liGP_n_smooth$nu * liGP_n_smooth$var), "estimated log var liGP")
# plot_pointillist(locs, log(tau), "true log var")
# par(mfrow=  c(1, 1))
# par(mar = c(4,4,1,1))
# plot(log(liGP_n_smooth$nu * liGP_n_smooth$var), log(tau), 
#      xlab = "estimated log nugget", ylab = "true log nugget")
# abline(a = 0, b=1)


print("Local GP done !")

########
# Done # 
########


saveRDS(benchmark, paste("het", seed, ".RDS", sep = ""))
