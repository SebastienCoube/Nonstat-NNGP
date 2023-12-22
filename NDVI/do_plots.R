#################################
# Options to set before the run #
#################################

setwd("NDVI")


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


mcmc_nngp_list = readRDS(paste("res/", run_name, sep=""))
print(run_name)

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

saveRDS(object = list(
  estimated_fixed = estimated_fixed, 
  estimated_parameters = estimated_parameters, 
  estimation_noise = estimation_noise, 
  preds_field = preds_field, 
  preds_fixed = preds_fixed, 
  preds_noise = preds_noise
  ), "preds_for_map.RDS")

############
# Plotting #
############

list2env(readRDS("preds_for_map.RDS"), environment())
PP = mcmc_nngp_list$hierarchical_model$PP

load(file ="data_cleaned_small_expanded.RData")
response_variable = data_cleaned_small$NDVI
response_variable = as.vector(scale(response_variable))
locs = cbind(data_cleaned_small$scaled_x, data_cleaned_small$scaled_y)
pdf("observed.pdf")
Bidart::plot_pointillist_painting(
  locs,
  response_variable, 
  cex = .1
)
dev.off()


pdf("noise.pdf")
Bidart::plot_pointillist_painting(
  rbind(mcmc_nngp_list$data$observed_locs, predicted_locs),
  c(estimation_noise$summaries["mean",,], preds_noise$summaries["mean",,]), 
  cex = .1
)
dev.off()

pdf("log_sd.pdf")
Bidart::plot_pointillist_painting(
  rbind(mcmc_nngp_list$data$locs, predicted_locs),
  log(c(estimated_parameters$summaries$field["sd",,], preds_field$summaries$field["sd",,])), 
  cex = .1
)
dev.off()

pdf("ellipses.pdf")
Bidart::plot_pointillist_painting(
  rbind(mcmc_nngp_list$data$locs, predicted_locs),
  c(estimated_parameters$summaries$field["mean",,], preds_field$summaries$field["mean",,]), 
  cex = .1
)

locs_4_ellipses = as.matrix(expand.grid(seq(.02, 1.08, .05), seq(.02, 1.08, .05)))
locs_4_ellipses=  split(locs_4_ellipses, row(locs_4_ellipses))
dists = sapply(locs_4_ellipses, function(x)min(fields::rdist(predicted_locs, matrix(x, 1))))
locs_4_ellipses = locs_4_ellipses[which(dists<.01)] 
locs_4_ellipses = do.call(rbind, locs_4_ellipses)
locs_4_ellipses_idx = apply(locs_4_ellipses, 1, function(x)which.min(fields::rdist(predicted_locs, matrix(x, 1))))
Bidart::plot_ellipses(
  locs_4_ellipses, 
  preds_field$summaries$log_range["mean",locs_4_ellipses_idx,], 
  add=T, shrink = .02
)

dev.off()
