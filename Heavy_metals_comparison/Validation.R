
# load the runs
runs = list()
i = 1
TF = c(T, F)
for(nonstat_noise in TF)
{
for(nonstat_range in TF)
{
  for(nonstat_scale in TF)
    {
      runs[[i]] = readRDS(paste("Heavy_metals_comparison/run_noise", nonstat_noise, "scale", nonstat_scale, "range", nonstat_range, ".RDS", sep = "_"))
      i = i+1
    }
  }
}

# load validation data set
test_data_set = readRDS("Heavy_metals_comparison/test_data_set.RDS")
train_data_set = readRDS("Heavy_metals_comparison/validation_train.RDS")

#######
# DIC #
#######

DICs = parallel::mclapply(runs, Bidart::DIC, mc.cores = length(runs), burn_in = .33)

##################################
# Smoothing (observed locations) #
###################################

# estimate the parameters (in particular the latent field at the observed locations)
params_estimation = parallel::mclapply(runs, Bidart::estimate_parameters, 
                                       get_samples = T, burn_in = .33,# MCMC samples returned by the function
                                       mc.cores = 4)
gc()
# predict the fixed effects at observed locations
fixed_effects_smooth = parallel::mclapply(runs,
                                          function(x)Bidart::predict_fixed_effects(
                                            mcmc_nngp_list = x, # MCMC samples
                                            X_pred = x$data$covariates$X$arg, # covariates used for prediction, here observe covariates
                                            burn_in = .33
                                            ), 
                                          mc.cores = 4)
gc()

# predict the noise at observed locations
noise_smooth = lapply(
  runs, function(x)
    Bidart::predict_noise(
    mcmc_nngp_list = x, 
    X_noise_pred = x$data$covariates$noise_X$arg, 
    predicted_locs = x$data$observed_locs,
    burn_in = .33
    )
)
gc()

smooth_log_score = parallel::mcmapply(FUN = Bidart::log_score_Gaussian, 
       observed_field = lapply(runs, function(x)x$data$observed_field), 
       latent_field_sample = lapply(params_estimation, function(x)x$samples$field_at_observed_locs[,,]), 
       fixed_effects_samples = lapply(fixed_effects_smooth, function(x)x$predicted_samples[,,]), 
       log_noise_samples = lapply(noise_smooth, function(x)x$predicted_samples), 
       mc.cores = length(runs), SIMPLIFY = F
)
gc()

lapply(smooth_log_score, function(x)x$total)

#####################################
# Predicting (Validation locations) #
#####################################


# estimate the parameters (in particular the latent field at the observed locations)
latent_field_preds = lapply(seq(length(runs)), 
                            function(x) 
                            {
                              print(x)
                              x = runs[[x]]
                              X_range_pred = NULL
                              X_scale_pred = NULL
                              if(!all(x$data$covariates$range_X$arg=="No covariates were provided"))X_range_pred = test_data_set$test_X[, c("gcarb", "globedem", "twi")]                             
                              if(!all(x$data$covariates$scale_X$arg=="No covariates were provided"))X_scale_pred = test_data_set$test_X[, c("gcarb", "globedem", "twi")]                             
                              Bidart::predict_latent_field (
                                mcmc_nngp_list = x, predicted_locs = test_data_set$test_locs, n_cores = 8, 
                                X_range_pred = X_range_pred,
                                X_scale_pred = X_scale_pred,
                                burn_in = .33
                                )
                            })

gc()
                              
# predict the fixed effects at observed locations
fixed_effects_pred = parallel::mclapply(runs,
                                        function(x)Bidart::predict_fixed_effects(
                                          mcmc_nngp_list = x, 
                                          X_pred = test_data_set$test_X, 
                                          burn_in = .33
                                            ),  mc.cores = length(runs))
gc()

# predict the noise at observed locations
noise_pred = lapply(runs,
                                function(x)
                                {
                                  X_noise_pred = test_data_set$test_X
                                  if(x$data$covariates$noise_X$arg == "No covariates were provided") X_noise_pred = NULL
                                  if(x$data$covariates$noise_X$arg != "No covariates were provided") if(ncol(x$data$covariates$noise_X$arg)>12) X_noise_pred = cbind(as.matrix(test_data_set$test_X), Bidart::predict_KL_basis(predicted_locs = test_data_set$test_locs, KL_basis = train_data_set$train_KL_decomposition)%*%diag(train_data_set$train_KL_decomposition$KL_decomposition$d))
                                  Bidart::predict_noise(
                                  mcmc_nngp_list = x, # MCMC samples
                                  predicted_locs = test_data_set$test_locs,
                                  X_noise_pred = X_noise_pred, # covariates used for prediction, here observed covariates
                                  burn_in = .33
                                  )
                                })
gc()


pred_log_score = parallel::mcmapply(FUN = Bidart::log_score_Gaussian, 
                                      observed_field = list(test_data_set$test_obs), 
                                      latent_field_sample = lapply(latent_field_preds, function(x) x$predicted_samples$field[,]), 
                                      fixed_effects_samples = lapply(fixed_effects_pred, function(x)x$predicted_samples[,,]), 
                                      log_noise_samples = lapply(noise_pred, function(x)x$predicted_samples[,,]), 
                                      mc.cores = length(runs), SIMPLIFY = F
)
gc()

lapply(pred_log_score, function(x)x$total)
barplot(sapply(pred_log_score, function(x)x$total))


saveRDS(
cbind(
model = c("NRS", "NR", "NS", "N", "RS", "R", "S", "stat"),
pred = round(sapply(pred_log_score, function(x)x$total)), 
smooth = round(sapply(smooth_log_score, function(x)x$total)), 
DIC = round(unlist(DICs)),
time= unname(sapply(runs, function(x)round(x$iterations$checkpoints[nrow(x$iterations$checkpoints), 2])))
), "Heavy_metals_comparison/NNGP_comparison.RDS"
)

