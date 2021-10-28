run_stat_validation = readRDS("Heavy_metals/run_stat_validation.RDS")
run_full_validation = readRDS("Heavy_metals/run_nr_basis_validation.RDS")
run_hetero_validation = readRDS("Heavy_metals/run_n_validation.RDS")
test_data_set = readRDS("Heavy_metals/test_data_set.RDS")
KL_basis_test = as.matrix(test_data_set$test_KL)[,seq(25)]

###################################
# Estimating and predicting stuff #
###################################

# full nonstat 

# observed locations 
run_full_estimation              = Bidart::estimate_parameters  (mcmc_nngp_list = run_full_validation, get_samples = T)
run_full_fixed_effect_smooth     = Bidart::predict_fixed_effects(mcmc_nngp_list = run_full_validation, run_full_validation$data$covariates$X$arg, individual_fixed_effects = NULL)
run_full_noise_smooth            = Bidart::predict_noise        (mcmc_nngp_list = run_full_validation, X_noise_pred = run_full_validation$data$covariates$noise_X$arg , individual_fixed_effects = NULL)

# predicted locations
run_full_noise_prediction        = Bidart::predict_noise        (mcmc_nngp_list = run_full_validation, X_noise_pred = cbind(KL_basis_test, test_data_set$test_X), individual_fixed_effects = NULL)
run_full_fixed_effect_prediction = Bidart::predict_fixed_effects(mcmc_nngp_list = run_full_validation, test_data_set$test_X, individual_fixed_effects = NULL)
run_full_latent_field_prediction = Bidart::predict_latent_field (mcmc_nngp_list = run_full_validation, predicted_locs = test_data_set$test_locs, n_cores = 3, 
  X_range_pred = cbind(KL_basis_test, test_data_set$test_X[, c("gcarb", "globedem", "twi")]), predict_range = T)

# noise only nonstat model

# observed locations 
run_hetero_estimation              = Bidart::estimate_parameters  (mcmc_nngp_list = run_hetero_validation, get_samples = T)
run_hetero_fixed_effect_smooth     = Bidart::predict_fixed_effects(mcmc_nngp_list = run_hetero_validation, run_hetero_validation$data$covariates$X$arg, individual_fixed_effects = NULL)
run_hetero_noise_smooth            = Bidart::predict_noise        (mcmc_nngp_list = run_hetero_validation, X_noise_pred = run_hetero_validation$data$covariates$noise_X$arg , individual_fixed_effects = NULL)

# predicted locations
run_hetero_noise_prediction        = Bidart::predict_noise        (mcmc_nngp_list = run_hetero_validation, X_noise_pred = cbind(KL_basis_test, test_data_set$test_X), individual_fixed_effects = NULL)
run_hetero_fixed_effect_prediction = Bidart::predict_fixed_effects(mcmc_nngp_list = run_hetero_validation, test_data_set$test_X, individual_fixed_effects = NULL)
run_hetero_latent_field_prediction = Bidart::predict_latent_field (mcmc_nngp_list = run_hetero_validation, predicted_locs = test_data_set$test_locs, n_cores = 3)


# stat model

# observed locations 
run_stat_estimation              = Bidart::estimate_parameters  (mcmc_nngp_list = run_stat_validation, get_samples = T)
run_stat_fixed_effect_smooth     = Bidart::predict_fixed_effects(mcmc_nngp_list = run_stat_validation, run_stat_validation$data$covariates$X$arg, individual_fixed_effects = NULL)
run_stat_noise_smooth            = Bidart::predict_noise        (mcmc_nngp_list = run_stat_validation, individual_fixed_effects = NULL)

# predicted locations
run_stat_noise_prediction        = Bidart::predict_noise        (mcmc_nngp_list = run_stat_validation, individual_fixed_effects = NULL)
run_stat_fixed_effect_prediction = Bidart::predict_fixed_effects(mcmc_nngp_list = run_stat_validation, test_data_set$test_X, individual_fixed_effects = NULL)
run_stat_latent_field_prediction = Bidart::predict_latent_field (mcmc_nngp_list = run_stat_validation, predicted_locs = test_data_set$test_locs, n_cores = 3)



####################
# Comparing models #
####################

#sum(dnorm(
#  x = test_data_set$test_obs,
#  mean = run_full_latent_field_prediction$summaries$field[1,] + run_full_fixed_effect_prediction$summaries$total_linear_effects[1,,], 
#  sd = exp(.5 * run_full_noise_prediction$summaries$total_linear_effects[1,,]), 
#  log = T
#))
#
#sum(dnorm(
#  x = test_data_set$test_obs,
#  mean = run_stat_latent_field_prediction$summaries$field[1,] + run_stat_fixed_effect_prediction$summaries$total_linear_effects[1,,], 
#  sd = exp(.5 * run_stat_noise_prediction$summaries$total_linear_effects[1,,]), 
#  log = T
#))


run_stat_pred_log_score = Bidart::log_score_Gaussian(
  observed_field         = test_data_set$test_obs, 
  latent_field_samples   = run_stat_latent_field_prediction$predicted_samples$field[,], 
  fixed_effects_samples  = run_stat_fixed_effect_prediction$predicted_samples$total_linear_effects[,,], 
  log_noise_samples      = run_stat_noise_prediction$predicted_samples$total_linear_effects[,,]
)

run_full_pred_log_score = Bidart::log_score_Gaussian(
  observed_field =        test_data_set$test_obs, 
  latent_field_samples =  run_full_latent_field_prediction$predicted_samples$field[,], 
  fixed_effects_samples = run_full_fixed_effect_prediction$predicted_samples$total_linear_effects[,,], 
  log_noise_samples =     run_full_noise_prediction$predicted_samples$total_linear_effects[,,]
)
run_hetero_pred_log_score = Bidart::log_score_Gaussian(
  observed_field =        test_data_set$test_obs, 
  latent_field_samples =  run_hetero_latent_field_prediction$predicted_samples$field[,], 
  fixed_effects_samples = run_hetero_fixed_effect_prediction$predicted_samples$total_linear_effects[,,], 
  log_noise_samples =     run_hetero_noise_prediction       $predicted_samples$total_linear_effects[,,]
)



run_stat_smooth_log_score = Bidart::log_score_Gaussian(
  observed_field =        run_stat_validation$data$observed_field, 
  latent_field_samples =  run_stat_estimation$samples$field_at_observed_locs[,,], 
  fixed_effects_samples = run_stat_fixed_effect_smooth$predicted_samples$total_linear_effects[,,], 
  log_noise_samples =     run_stat_noise_smooth       $predicted_samples$total_linear_effects[,,]
)
run_full_smooth_log_score = Bidart::log_score_Gaussian(
  observed_field =        run_full_validation$data$observed_field, 
  latent_field_samples =  run_full_estimation$samples$field_at_observed_locs[,,], 
  fixed_effects_samples = run_full_fixed_effect_smooth$predicted_samples$total_linear_effects[,,], 
  log_noise_samples =     run_full_noise_smooth       $predicted_samples$total_linear_effects[,,]
)
run_hetero_smooth_log_score = Bidart::log_score_Gaussian(
  observed_field =        run_hetero_validation$data$observed_field, 
  latent_field_samples =  run_hetero_estimation$samples$field_at_observed_locs[,,], 
  fixed_effects_samples = run_hetero_fixed_effect_smooth$predicted_samples$total_linear_effects[,,], 
  log_noise_samples =     run_hetero_noise_smooth       $predicted_samples$total_linear_effects[,,]
)



tests_result = matrix(
  c(
    run_full_pred_log_score$total,
    run_hetero_pred_log_score$total,
    run_stat_pred_log_score$total,
    run_full_smooth_log_score$total,
    run_hetero_smooth_log_score$total,
    run_stat_smooth_log_score$total,
    Bidart::DIC(run_full_validation),
    Bidart::DIC(run_hetero_validation),
    Bidart::DIC(run_stat_validation)
  ), 3
)


row.names(tests_result) = c("range + noise", "noise",  "stationary")
colnames(tests_result) = c("pred. log score", "smooth. log score", "DIC")




plot(test_data_set$test_locs, col = 1+(run_full_pred_log_score$per_obs>run_stat_pred_log_score$per_obs), pch = 16, cex = .5)
plot(run_full_validation$data$observed_locs, col = 1+(run_full_smooth_log_score$per_obs>run_stat_smooth_log_score$per_obs), pch = 16, cex = .5)

plot(test_data_set$test_locs, col = 1+(run_hetero_pred_log_score$per_obs>run_stat_pred_log_score$per_obs), pch = 16, cex = .5)

plot(test_data_set$test_locs, col = 1+(run_full_pred_log_score$per_obs>run_hetero_pred_log_score$per_obs), pch = 16, cex = .5)

color_scale = run_full_latent_field_prediction$summaries$log_range[1,]
color_scale = 1 + floor((color_scale - min(color_scale))/(max(color_scale)-min(color_scale))*100)
plot(test_data_set$test_locs, col = heat.colors(100)[color_scale], pch = 16)
points(test_data_set$test_locs, col = 1+(run_full_pred_log_score$per_obs>run_stat_pred_log_score$per_obs), pch = 16, cex = .5)

color_scale = run_full_noise_smooth$summaries$total_linear_effects[1,,]
color_scale = 1 + floor((color_scale - min(color_scale))/(max(color_scale)-min(color_scale))*100)
plot(run_full_validation$data$observed_locs, col = heat.colors(100)[color_scale], pch = 16)
points(run_full_validation$data$observed_locs, col = 1+(run_full_smooth_log_score$per_obs>run_stat_smooth_log_score$per_obs), pch = 16, cex = .5)
plot(run_full_validation$data$observed_locs, col = 1+(run_full_smooth_log_score$per_obs>run_stat_smooth_log_score$per_obs), pch = 16, cex = .5)

smooth_win = (run_full_smooth_log_score$per_obs>run_stat_smooth_log_score$per_obs)
smooth_noise = scale(run_full_noise_smooth$summaries$total_linear_effects[1,,])
smooth_range = scale(as.matrix(cbind(1, run_full_validation$data$covariates$range_X$arg)) %*% run_full_estimation$summaries$range_beta[1,,])
summary(glm(smooth_win~smooth_noise+smooth_range, family = "binomial"))

pred_win =  (run_full_pred_log_score$per_obs>run_stat_pred_log_score$per_obs)
pred_noise = scale(run_full_noise_prediction$summaries$total_linear_effects[1,,])
pred_range = scale(run_full_latent_field_prediction$summaries$log_range[1,])
summary(glm(smooth_win~smooth_noise+smooth_range, family = "binomial"))


par(mfrow = c(3, 3))
for(i in sample(seq(length(test_data_set$test_obs)), 9))
{
  if(mean(run_full_pred_log_score$samples[i,])>mean(run_stat_pred_log_score$samples[i,]))main = paste(i, "nonstat wins")
  if(mean(run_full_pred_log_score$samples[i,])<mean(run_stat_pred_log_score$samples[i,]))main = paste(i, "stat wins")
  plot(density(run_full_pred_log_score$samples[i,]), 
       xlim = c(min(density(run_full_pred_log_score$samples[i,])$x, density(run_stat_pred_log_score$samples[i,])$x), max(density(run_full_pred_log_score$samples[i,])$x, density(run_stat_pred_log_score$samples[i,])$x)),
       ylim = c(min(density(run_full_pred_log_score$samples[i,])$y, density(run_stat_pred_log_score$samples[i,])$y), max(density(run_full_pred_log_score$samples[i,])$y, density(run_stat_pred_log_score$samples[i,])$y)),
       main = main
  )
  lines(density(run_stat_pred_log_score$samples[i,]), col = 2)
  legend("topleft", legend = c("stat", "nonstat"), fill = c(2, 1))
}
par(mfrow = c(1, 1))

