#################################
# Options to set before the run #
#################################

setwd("NDVI")
plot_stuff = T


########################
# Loading the data set #
########################

predicted_locs = readRDS("test_locs.RDS")
# the problem is that the rounding merges locations and my script does not allow to predict @ observed locations
test_field = readRDS("test_field.RDS")


run_names = list.files("res/")[grep(".RDS", list.files("res/"))]
run_names = run_names[grep("run", run_names)]
run_names = run_names[-1]

run_name = run_names[1]

for(run_name in run_names){
  mcmc_nngp_list = readRDS(paste("res/", run_name, sep=""))
  print(run_name)
   ## just checking that train and test are matching 
   #plot(
   #  mcmc_nngp_list$data$observed_locs, 
   #  cex = .05
   #)
   #points(
   #  predicted_locs, 
   #  cex = .05, col=2
   #)
   #legend("topleft", legend = c("train", "test"), fill = c(1,2))
   
   ###############
   # Predictions #
   ###############
   # This part predicts the process at unobserved locations. 
   # It will be used for prediction error. 
   
   # Latent field 
   preds_field = Bidart::predict_latent_field(
     # mcmc samples 
     mcmc_nngp_list = mcmc_nngp_list, 
     # locations where predictions must be done, must be unobserved locations
     predicted_locs = predicted_locs, 
     # Those are the covariates used to predict the nonstationary covariance parameters
     # Since there are no covariates, the NULL input is used
     X_range_pred = NULL, 
     X_scale_pred = NULL, 
     # MCMC burn-in
     burn_in = .1, 
     # number of OMP process per socket
     num_threads_per_chain = 10)
   
   # Fixed effects 
   preds_fixed = Bidart::predict_fixed_effects(
     mcmc_nngp_list = mcmc_nngp_list, 
     # covariates at new locations
     X_pred = as.data.frame(predicted_locs), 
     burn_in = .1)
   
   # Noise 
   preds_noise = Bidart::predict_noise(
     mcmc_nngp_list = mcmc_nngp_list, 
     # covariates at new locations (none)
     X_noise_pred = NULL, 
     # new locations (used for GP prior)
     predicted_locs =  predicted_locs,
     burn_in = .1)
   
   ##############
   # Estimation #
   ##############
   # This part estimates the process at observed locations. 
   # It will be used for estimation/training error. 
   
   # Estimate the parameters, including the latent field at observed locations
   estimated_parameters = Bidart::estimate_parameters(
     mcmc_nngp_list = mcmc_nngp_list, 
     burn_in = .1, 
     get_samples = T
   )
   
   # "Predict" the fixed effects at observed locations
   estimated_fixed = Bidart::predict_fixed_effects(
     mcmc_nngp_list = mcmc_nngp_list, 
     # covariates at new locations
     X_pred = mcmc_nngp_list$data$covariates$X$arg, 
     burn_in = .1)
   
   # "Predict" the noise variance at observed locations 
   estimation_noise = Bidart::predict_noise(
     mcmc_nngp_list = mcmc_nngp_list, 
     # covariates at new locations (none)
     X_noise_pred = mcmc_nngp_list$data$covariates$noise_X$arg, 
     # new locations (used for GP prior)
     predicted_locs =  mcmc_nngp_list$data$observed_locs,
     burn_in = .1)
   
   
   ##########
   # Scores # 
   ##########
   model_DIC = Bidart::DIC(mcmc_nngp_list = mcmc_nngp_list, burn_in = .1)
   
   train_log_dens = Bidart::log_score_Gaussian(
     observed_field = mcmc_nngp_list$data$observed_field, 
     latent_field_samples = estimated_parameters$samples$field_at_observed_locs, 
     log_noise_samples = estimation_noise$predicted_samples, 
     fixed_effects_samples = estimated_fixed$predicted_samples
   )
   
   test_log_dens = Bidart::log_score_Gaussian(
     observed_field = test_field, 
     latent_field_samples = preds_field$predicted_samples$field, 
     log_noise_samples = preds_noise$predicted_samples, 
     fixed_effects_samples = preds_fixed$predicted_samples
   )
   
  
  sink("prediction_comparison.txt", append = T)
  print(c(
    run_name, 
    round(model_DIC), 
    round(train_log_dens$total/mcmc_nngp_list$vecchia_approx$n_obs, 2),
    round(test_log_dens$total/nrow(predicted_locs),2), 
    round(min(unlist(Bidart::ESS(mcmc_nngp_list = mcmc_nngp_list, burn_in = .1))))
    ))

  sink()
  
  
}
