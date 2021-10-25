
run_stat_validation = readRDS("Heavy_metals/run_stat_validation.RDS")
run_noise_range_validation = readRDS("Heavy_metals/run_nr_basis_validation.RDS")
run_noise_validation = readRDS("Heavy_metals/run_n_validation.RDS")
test_data_set = readRDS("Heavy_metals/test_data_set.RDS")

KL_basis = as.matrix(test_data_set$test_KL)

# full nonstat model
noise_range_estimation              = Bidart::estimate_parameters(run_noise_range_validation, get_samples = T)
noise_range_noise_prediction        = Bidart::predict_noise(mcmc_nngp_list = run_noise_range_validation, X_noise_pred = cbind(test_data_set$test_X, KL_basis), individual_fixed_effects = NULL)
noise_range_fixed_effect_prediction = Bidart::predict_fixed_effects(run_noise_range_validation, test_data_set$test_X, individual_fixed_effects = NULL)
noise_range_latent_field_prediction = Bidart::predict_latent_field(mcmc_nngp_list = run_noise_range_validation, predicted_locs = test_data_set$test_locs[!duplicated(test_data_set$test_locs),], n_cores = 3, X_range_pred = cbind(test_data_set$test_X, as.matrix(KL_basis))[!duplicated(test_data_set$test_locs),], X_scale_pred = cbind(test_data_set$test_X, as.matrix(KL_basis))[!duplicated(test_data_set$test_locs),], predict_range = T)

# noise only nonstat model
noise_estimation              = Bidart::estimate_parameters(run_noise_range_validation)
noise_noise_prediction        = Bidart::predict_noise(mcmc_nngp_list = run_noise_range_validation, X_noise_pred = cbind(test_data_set$test_X, KL_basis), individual_fixed_effects = NULL)
noise_fixed_effect_prediction = Bidart::predict_fixed_effects(run_noise_range_validation, test_data_set$test_X, individual_fixed_effects = NULL)
noise_latent_field_prediction = Bidart::predict_latent_field(mcmc_nngp_list = run_noise_range_validation, predicted_locs = test_data_set$test_locs[!duplicated(test_data_set$test_locs),], n_cores = 3, X_range_pred = cbind(test_data_set$test_X, as.matrix(KL_basis))[!duplicated(test_data_set$test_locs),], X_scale_pred = cbind(test_data_set$test_X, as.matrix(KL_basis))[!duplicated(test_data_set$test_locs),], predict_range = T)

# stat model
stat_noise_prediction        = Bidart::predict_noise(run_stat_validation, individual_fixed_effects = NULL)
stat_estimation              = Bidart::estimate_parameters(run_stat_validation)
stat_fixed_effect_prediction = Bidart::predict_fixed_effects(run_stat_validation, test_data_set$test_X, individual_fixed_effects = NULL)
stat_latent_field_prediction = Bidart::predict_latent_field(mcmc_nngp_list = run_stat_validation, predicted_locs = test_data_set$test_locs[!duplicated(test_data_set$test_locs),], n_cores = 3)



sum(dnorm(
  test_data_set$test_obs, 
  mean = stat_fixed_effect_prediction$summaries$total_linear_effects[1,,] + stat_latent_field_prediction$summaries$field[1,idx,], 
  sd = exp(.5 * stat_noise_prediction$summaries$total_linear_effects[1,,]),
  log = T))

sum(dnorm(
  test_data_set$test_obs, 
  mean = nonstat_fixed_effect_prediction$summaries$total_linear_effects[1,,] + nonstat_latent_field_prediction$summaries$field[1,idx,], 
  sd = exp(.5 * nonstat_noise_prediction$summaries$total_linear_effects[1,,]),
  log = T))


Bidart::log_score_Gaussian(
  observed_field = test_data_set$test_obs, 
  latent_field_samples = stat_latent_field_prediction$summaries$field[1,idx,], 
  fixed_effects_samples = stat_fixed_effect_prediction$summaries$total_linear_effects[1,,], 
  log_noise_samples = stat_noise_prediction$summaries$total_linear_effects[1,,])
Bidart::log_score_Gaussian(
  observed_field =        test_data_set$test_obs, 
  latent_field_samples =  nonstat_latent_field_prediction$summaries$field[1,idx,], 
  fixed_effects_samples = nonstat_fixed_effect_prediction$summaries$total_linear_effects[1,,], 
  log_noise_samples =     nonstat_noise_prediction$summaries$total_linear_effects[1,,])

test = Bidart::log_score_Gaussian(
  observed_field = test_data_set$test_obs, 
  latent_field_samples = 
    do.call(cbind, list(
      nonstat_latent_field_prediction$predicted_samples$chain_1$field[idx,,], 
      nonstat_latent_field_prediction$predicted_samples$chain_2$field[idx,,], 
      nonstat_latent_field_prediction$predicted_samples$chain_3$field[idx,,]
      )),
  log_noise_samples = 
    do.call(cbind, list(
      nonstat_noise_prediction$predicted_samples$chain_1$total_linear_effects[,1,], 
      nonstat_noise_prediction$predicted_samples$chain_2$total_linear_effects[,1,], 
      nonstat_noise_prediction$predicted_samples$chain_3$total_linear_effects[,1,]
      )),
  fixed_effects_samples =  
    do.call(cbind, list(
      nonstat_fixed_effect_prediction$predicted_samples$chain_1$total_linear_effects[,1,], 
      nonstat_fixed_effect_prediction$predicted_samples$chain_2$total_linear_effects[,1,], 
      nonstat_fixed_effect_prediction$predicted_samples$chain_3$total_linear_effects[,1,]
    ))
)

test_ = Bidart::log_score_Gaussian(
  observed_field = test_data_set$test_obs, 
  latent_field_samples = 
    do.call(cbind, list(
      stat_latent_field_prediction$predicted_samples$chain_1$field[idx,,], 
      stat_latent_field_prediction$predicted_samples$chain_2$field[idx,,], 
      stat_latent_field_prediction$predicted_samples$chain_3$field[idx,,]
    )),
  log_noise_samples = 
    do.call(cbind, list(
      stat_noise_prediction$predicted_samples$chain_1$total_linear_effects[,1,], 
      stat_noise_prediction$predicted_samples$chain_2$total_linear_effects[,1,], 
      stat_noise_prediction$predicted_samples$chain_3$total_linear_effects[,1,]
    )),
  fixed_effects_samples =  
    do.call(cbind, list(
      stat_fixed_effect_prediction$predicted_samples$chain_1$total_linear_effects[,1,], 
      stat_fixed_effect_prediction$predicted_samples$chain_2$total_linear_effects[,1,], 
      stat_fixed_effect_prediction$predicted_samples$chain_3$total_linear_effects[,1,]
    ))
)


test$total
test_$total
par(mfrow = c(3, 3))
for(i in sample(seq(length(test_data_set$test_obs)), 9))
{
  if(mean(test$samples[i,])>mean(test_$samples[i,]))main = paste(i, "nonstat wins")
  if(mean(test$samples[i,])<mean(test_$samples[i,]))main = paste(i, "stat wins")
  plot(density(test$samples[i,]), 
       xlim = c(min(density(test$samples[i,])$x, density(test_$samples[i,])$x), max(density(test$samples[i,])$x, density(test_$samples[i,])$x)),
       ylim = c(min(density(test$samples[i,])$y, density(test_$samples[i,])$y), max(density(test$samples[i,])$y, density(test_$samples[i,])$y)),
       main = main
  )
  lines(density(test_$samples[i,]), col = 2)
  legend("topleft", legend = c("stat", "nonstat"), fill = c(2, 1))
}
mean(test$samples[3181,])
mean(test_$samples[3181,])
dev.off()
pdf("Heavy_metals/correct_predictions.pdf")
plot(test_data_set$test_locs, col = 1 + (apply(test$samples, 1, mean)>apply(test_$samples, 1, mean)), pch = 16, cex = .5)
legend("bottomleft", legend = c("nonstat wins", "stat wins"), fill = c(2, 1))
dev.off()


Bidart::plot_pointillist_painting(test_data_set$test_locs, nonstat_latent_field_prediction$summaries$log_range[1,idx,])
points(test_data_set$test_locs [test$per_obs>test_$per_obs,], pch = 16, cex = .1)



Dist = as.matrix(dist(test_data_set$test_locs))
hist(Dist)
stoopid_kernel = apply(Dist<1, 1, mean)
Bidart::plot_pointillist_painting(test_data_set$test_locs, stoopid_kernel)
pred_range = nonstat_latent_field_prediction$summaries$log_range[1,idx,]
GLM = glm(as.numeric(test$per_obs>test_$per_obs)~pred_range+stoopid_kernel)
summary(GLM)
table(test$per_obs>test_$per_obs, GLM$fitted.values>.5)
sum(test$per_obs>test_$per_obs)
sum(test$per_obs<test_$per_obs)

Bidart::plot_pointillist_painting(test_data_set$test_locs, GLM$fitted.values * (GLM$fitted.values>.5))
Bidart::plot_pointillist_painting(test_data_set$test_locs, pred_range)

plot(stoopid_kernel,  GLM$fitted.values)
plot(pred_range,  GLM$fitted.values)
abline(h = .5)

mean(GLM$fitted.values>.5)


