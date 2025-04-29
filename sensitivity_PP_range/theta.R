setwd("sensitivity_PP_range/")

lib.loc = NULL


library(FNN    , lib.loc = lib.loc)
library(fields , lib.loc = lib.loc)
library(ellipse, lib.loc = lib.loc)
library(coda   , lib.loc = lib.loc)
library(GpGp   , lib.loc = lib.loc)
library(Bidart , lib.loc = lib.loc)


###################
# Simulating data #
###################


matern_range_vals = .05*sqrt(2)^seq(0, 4)
seeds = seq(30)
range_design = as.matrix(expand.grid(matern_range_vals,seeds))
saveRDS(range_design, "range_design.RDS")

for(i in seq(nrow(range_design))){
  if (!(paste("range_exp_started_", i, sep = "") %in% list.files())){
    saveRDS(1, paste("range_exp_started_", i, sep = ""))
    seed = range_design[i,2]
    matern_range = range_design[i,1]
    set.seed(seed)
    # spatial locations
    locs = cbind(runif(11000), runif(11000))
    # heteroscedastic noise variance
    range_var = 1
    range_range = .1
    log_range = (GpGp::fast_Gp_sim(c(.3, range_range, 1.5, 0), locs = rbind(locs), m=10) -  4)
    # adding noise to get observations
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
    z = (w[seq(nrow(locs))] + .5 * rnorm(nrow(locs)))
    # splitting train test
    predicted_locs = locs[-seq(10000),]
    predicted_z = z[-seq(10000)]
    locs = locs[seq(10000),]
    z = z[seq(10000)]
    log_range = log_range[seq(10000)]
    
    # visualizing data
    plot_pointillist = function(locs, x, main)plot(locs, xlab = "", ylab = "", main = main, pch  =15, col = heat.colors(n = 100)[as.numeric(cut(x, breaks = seq(min(x), max(x), length.out = 98)))], axes = F)
    plot_pointillist(locs, w, "noiseless field")
    plot_pointillist(locs, log_range, "log range of the field")
    plot_pointillist(locs, z, "observed noisy field")
    
    t1 = Sys.time()
    PP = Bidart::get_PP(observed_locs = locs, matern_range = matern_range, n_PP = 150)
    Bidart::compare_PP_NNGP(PP)
    mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
      observed_locs = locs, observed_field = z, 
      range_PP = T, PP = PP, range_log_scale_prior = c(-4, 2.5)
    )
    
    
    for(j in seq(15))mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(mcmc_nngp_list, num_threads_per_chain = 5, thinning = 1, lib.loc = lib.loc)
    
    estimation = Bidart::estimate_parameters(mcmc_nngp_list, get_samples = T, burn_in = .2)
    
    
    estimated_log_range_NNGP = 
      Bidart::X_PP_mult_right(
        X = matrix(1, nrow(locs)), 
        PP = PP, use_PP = T, 
        Y = estimation$summaries$range_beta[1,,]) 
    par(mfrow = c(1,2))
    plot_pointillist(locs, estimated_log_range_NNGP, main = "estimated log range")
    plot_pointillist(locs, log_range, main = "true log range")
    
    
    fixed_effects = Bidart::predict_fixed_effects(mcmc_nngp_list, burn_in = .2, lib.loc = lib.loc)
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
      num_threads_per_chain = 5, lib.loc = lib.loc, 
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
    
    
    res = c(
      w_smooth_MSE = mean((estimation$summaries$beta[1,1,1] + estimation$summaries$field_at_observed_locs[1,,] - w[seq(nrow(locs))])^2),
      w_pred_MSE =   mean((estimation$summaries$beta[1,1,1] + NNGP_pred$summaries$field[1,,] - w[-seq(nrow(locs))])^2), 
      log_range_MSE = mean((log_range- estimated_log_range_NNGP)^2),
      log_range_scale =  estimation$summaries$range_log_scale[1], 
      obs_coverage = coverage, 
      pred_coverage = pred_coverage,
      obs_elpd = NNGP_obs_elpd, 
      pred_elpd = NNGP_pred_elpd,
      time = as.numeric(Sys.time()-t1, unit = "mins"), 
      seed = seed
    )
    saveRDS(res, paste("range_exp_ended_", i, sep = ""))
  }
}
