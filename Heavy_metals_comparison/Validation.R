
# load the runs
runs = list(
stat_exp = readRDS("Heavy_metals_comparison/run_stat_exp_validation.RDS"),
stat_15 = readRDS("Heavy_metals_comparison/run_stat_15_validation.RDS"),
hetero_exp = readRDS("Heavy_metals_comparison/run_n_exp_validation.RDS"),
hetero_15 = readRDS("Heavy_metals_comparison/run_n_15_validation.RDS"),
noise_range_exp = readRDS("Heavy_metals_comparison/run_nr_exp_basis_validation.RDS"),
noise_range_15 = readRDS("Heavy_metals_comparison/run_nr_15_basis_validation.RDS"),
full_exp = readRDS("Heavy_metals_comparison/run_nsr_exp_basis_validation.RDS"),
full_15 = readRDS("Heavy_metals_comparison/run_nsr_15_basis_validation.RDS")
)

# load validation data set
test_data_set = readRDS("Heavy_metals_comparison/test_data_set.RDS")

# load KL basisat validation locations
KL_basis_test = as.matrix(test_data_set$test_KL)[,seq(25)]


#######
# DIC #
#######

DICs = parallel::mclapply(runs, Bidart::DIC, mc.cores = length(runs), burn_in = .25)

##################################
# Smoothing (observed locations) #
###################################

# estimate the parameters (in particular the latent field at the observed locations)
params_estimation = parallel::mclapply(runs, Bidart::estimate_parameters, 
                                       get_samples = T, burn_in = .25,# MCMC samples returned by the function
                                       mc.cores = 3)
gc()
# predict the fixed effects at observed locations
fixed_effects_smooth = parallel::mclapply(runs,
                                          function(x)Bidart::predict_fixed_effects(
                                            mcmc_nngp_list = x, # MCMC samples
                                            X_pred = x$data$covariates$X$arg, # covariates used for prediction, here observe covariates
                                            individual_fixed_effects = NULL, burn_in = .25
                                            ), 
                                          mc.cores = 3)
gc()

# predict the noise at observed locations
noise_smooth = parallel::mclapply(
  runs, function(x)
    Bidart::predict_noise(
    mcmc_nngp_list = x, 
    X_noise_pred = x$data$covariates$noise_X$arg, 
    individual_fixed_effects = NULL, burn_in = .25
    )
  , mc.cores = 2
)
gc()

smooth_log_score = parallel::mcmapply(FUN = Bidart::log_score_Gaussian, 
       observed_field = lapply(runs, function(x)x$data$observed_field), 
       latent_field_sample = lapply(params_estimation, function(x)x$samples$field_at_observed_locs[,,]), 
       fixed_effects_samples = lapply(fixed_effects_smooth, function(x)x$predicted_samples$total_linear_effects[,,]), 
       log_noise_samples = lapply(noise_smooth, function(x)x$predicted_samples$total_linear_effects[,,]), 
       mc.cores = length(runs), SIMPLIFY = F
)
gc()

lapply(smooth_log_score, function(x)x$total)

#####################################
# Predicting (Validation locations) #
#####################################

# estimate the parameters (in particular the latent field at the observed locations)
latent_field_preds = lapply(runs, 
                            function(x) 
                            {
                              X_range_pred = NULL
                              if(x$data$covariates$range_X$arg!="No covariates were provided")X_range_pred = cbind(KL_basis_test, test_data_set$test_X[, c("gcarb", "globedem", "twi")])
                              X_scale_pred = NULL
                              if(x$data$covariates$scale_X$arg!="No covariates were provided")X_scale_pred = cbind(KL_basis_test, test_data_set$test_X[, c("gcarb", "globedem", "twi")])
                              Bidart::predict_latent_field (
                                mcmc_nngp_list = x, predicted_locs = test_data_set$test_locs, n_cores = 3, 
                                X_range_pred = X_range_pred, X_scale_pred =  X_scale_pred,
                                predict_range = T, predict_scale = T, burn_in = .25
                                )
                            })

gc()
                              
# predict the fixed effects at observed locations
fixed_effects_pred = parallel::mclapply(runs,
                                        function(x)Bidart::predict_fixed_effects(
                                          mcmc_nngp_list = x, # MCMC samples
                                          X_pred = test_data_set$test_X, # covariates used for prediction, here observe covariates
                                            individual_fixed_effects = NULL, burn_in = .25
                                            ),  mc.cores = length(runs))
gc()

# predict the noise at observed locations
noise_pred = parallel::mclapply(runs,
                                function(x)
                                {
                                  X_noise_pred = "No covariates were provided"
                                  if(x$data$covariates$noise_X$arg!="No covariates were provided")X_noise_pred = cbind(KL_basis_test, test_data_set$test_X)
                                  Bidart::predict_noise(
                                  mcmc_nngp_list = x, # MCMC samples
                                  X_noise_pred = X_noise_pred, # covariates used for prediction, here observed covariates
                                  individual_fixed_effects = NULL, burn_in = .25
                                  )
                                }
                                , mc.cores = length(runs))
gc()


pred_log_score = parallel::mcmapply(FUN = Bidart::log_score_Gaussian, 
                                      observed_field = list(test_data_set$test_obs), 
                                      latent_field_sample = lapply(latent_field_preds, function(x) x$predicted_samples$field[,]), 
                                      fixed_effects_samples = lapply(fixed_effects_pred, function(x)x$predicted_samples$total_linear_effects[,,]), 
                                      log_noise_samples = lapply(noise_pred, function(x)x$predicted_samples$total_linear_effects[,,]), 
                                      mc.cores = length(runs), SIMPLIFY = F
)
gc()

lapply(pred_log_score, function(x)x$total)

saveRDS(
cbind(
#model = rep(c("stat", "full", "hetero"), each = 2),
#smoothness = rep(c(".5", "1.5"), 3),
pred = round(sapply(pred_log_score, function(x)x$total)), 
smooth = round(sapply(smooth_log_score, function(x)x$total)), 
DIC = round(unlist(DICs))
), "Heavy_metals_comparison/NNGP_comparison.RDS"
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



### # re run with only nr and stat in order to see which model wins and where
### 
### 
### run_full_pred_log_score = pred_log_score[[2]]
### run_full_smooth_log_score = smooth_log_score[[2]]
### run_stat_pred_log_score = pred_log_score[[1]]
### run_stat_smooth_log_score = smooth_log_score[[1]]
### 
### smooth_win = (run_full_smooth_log_score$per_obs>run_stat_smooth_log_score$per_obs)
### smooth_noise = scale(noise_smooth[[2]]$summaries$total_linear_effects[1,,])
### smooth_range = scale(as.matrix(cbind(1, runs[[2]]$data$covariates$range_X$arg)) %*% params_estimation[[2]]$summaries$range_beta[1,,])
### summary(glm(smooth_win~smooth_noise+smooth_range, family = "binomial"))
### 
### pred_win =  (run_full_pred_log_score$per_obs>run_stat_pred_log_score$per_obs)
### pred_noise = scale(noise_pred[[2]]$summaries$total_linear_effects[1,,])
### pred_range = scale(latent_field_preds[[2]]$summaries$log_range[1,])
### summary(glm(pred_win~pred_noise+pred_range, family = "binomial"))

