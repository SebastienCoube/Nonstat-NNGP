#################################
# Options to set before the run #
#################################

## --------- IMPORTANT ---------- ##
# For some reason you have to set manually num_threads_per_chain in the function mcmc_nngp_run

# working directory
# setwd("NDVI_r2")

# rounding the coordinates allows to save computation effort at
# the expense of some accuracy of the coordinates
# this parameter controls the rounding of the locations. the smaller the more rounding
rounding_factor = 1/400
# Does some informative plots prior to MCMC run
plot_stuff = F

# load
.libPaths("../R_packages/")
library("Bidart")


########################
# Loading the data set #
########################
# load the data set
#load(file = "data_cleaned_small_expanded.RData")
load(file ="data_cleaned_small_expanded.RData")
response_variable = data_cleaned_small$NDVI
response_variable = as.vector(scale(response_variable))
locs = cbind(data_cleaned_small$scaled_x, data_cleaned_small$scaled_y)

########################
# rounding coordinates #
########################
locs  = rounding_factor * round(locs/rounding_factor)
print("total number of observations")
print(nrow(locs))
print("number of rounded locations")
print(sum(!duplicated(locs)))

#########################
# train-test separation #
#########################

set.seed(1)
# matching unique locs with all locs
unique_locs = do.call(rbind, unique(split(locs, row(locs))))
locs_match = match(split(locs, row(locs)), split(unique_locs, row(unique_locs)))
hctam_scol = lapply(seq(nrow(unique_locs)), function(i)c())
for(i in seq(length(locs_match)))hctam_scol[[locs_match[i]]] = c(hctam_scol[[locs_match[i]]], i)
  
# lumps 
km = kmeans(unique_locs, 6000, iter.max = 30)
plot(unique_locs, cex = .2, pch = 16, col = km$cluster, main = "KL clustering")
lumps_out = GpGp::order_maxmin(km$centers)[seq(500)]
lumps_centers_idx = apply(fields::rdist(km$centers[lumps_out,], unique_locs), 1, which.min)
lumps_out_idx = which(km$cluster %in% lumps_out)
# loo
loo_idx = seq(nrow(unique_locs))[-lumps_out_idx][GpGp::order_maxmin(unique_locs[-lumps_out_idx,])[seq(2000)]]


pdf("train_test_split.pdf")
# plotting train locs
plot(unique_locs[-c(lumps_out_idx, loo_idx),], cex = .2, pch = 16, col = 1 , xlab = "lon", ylab = "lat")
# removing lumps
points(unique_locs[lumps_out_idx,], cex = .2, pch = 16, col = 2)
# adding lumps centers
points(unique_locs[lumps_centers_idx,], cex = .3, pch = 16, col = 4)
# adding loo
points(unique_locs[loo_idx,], cex = .3, pch = 16, col = 7)
legend("topleft", fill = c(1,2,4, 7), legend = c("train", "removed", "test lump", "test pseudo-loo"))
dev.off()

# indexing data sets
loo_idx_expanded = unlist(hctam_scol[loo_idx])
lump_out_idx_expanded = unlist(hctam_scol[lumps_out_idx])
lump_center_idx_expanded = unlist(hctam_scol[lumps_centers_idx])
train_locs_idx_expanded = seq(nrow(locs))[-c(loo_idx_expanded, lump_out_idx_expanded)]

# saving
train_locs = locs[train_locs_idx_expanded,]
saveRDS(train_locs, "train_locs.RDS")
plot(train_locs, cex = .2, pch = 16)
train_field = response_variable[train_locs_idx_expanded]
saveRDS(train_field, "train_field.RDS")
loo_locs = locs[loo_idx_expanded,]
points(loo_locs, cex = .2, pch = 16, col=7)
saveRDS(loo_locs, "loo_locs.RDS")
saveRDS(response_variable[loo_idx_expanded],  "loo_field.RDS",)
lump_center_locs = locs[lump_center_idx_expanded,]
points(lump_center_locs, cex = .2, pch = 16, col=2)
saveRDS(lump_center_locs, "lump_center_locs.RDS")
saveRDS(response_variable[lump_center_idx_expanded],  "lump_center_field.RDS",)

####################
# Getting PP basis #
####################

PP = Bidart::get_PP(
  observed_locs = train_locs, # spatial sites
  matern_range = .2,
  n_PP = 50, # number of knots
  m = 15 # number of NNGP parents
  )
saveRDS(PP, "PP.RDS")

if(plot_stuff){
  # comparison between PP and NNGP
  Bidart::compare_PP_NNGP(PP)
}

#######
# Run #
#######

runs_design = data.frame(
  "nickname" = c("NNGP_full_monty", "NNGP_vanilla", "NNGP_heterosk", "NNGP_aniso_heterosk", "NNGP_lociso_heterosk"),
  "range_PP" = c(T, F, F, F, T),
  "aniso" =    c(T, F, F, T, F), 
  "noise_PP" = c(T, F, T, T, T)
) 

i_model=2
for(i_model in seq(nrow(runs_design))){
if(length(grep(paste(runs_design[i_model,1], "_lock", sep=""), list.files()))==0){
  t1 = Sys.time()
  # saving to lock slot
  saveRDS(1, paste(runs_design[i_model,1], "_lock", sep=""))
  ###########################
  # initializing the bousin #
  ###########################
  mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
    observed_locs = train_locs, 
    observed_field = train_field, 
    X = as.data.frame(train_locs), 
    m = 10, 
    nu = 1.5, 
    range_PP = runs_design[i_model,2], 
    noise_PP = runs_design[i_model,4], 
    anisotropic = runs_design[i_model,3], 
    sphere = F, 
    n_chains = 2, 
    PP = PP, 
    noise_log_scale_prior = c(-3, 5)
  )
  ######################
  # running the bousin #
  ######################
  for(i_pass in seq(20)){
    mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(
      mcmc_nngp_list,
      n_cores = 2, 
      num_threads_per_chain = ceiling(parallel::detectCores()/(2)), 
      plot_diags = T, 
      lib.loc = "../R_packages"
    )
  }
  
  # saving diagnostic plots of all the bousin
  pdf(paste(runs_design[i_model,1], ".pdf", sep=""))
  Bidart::diagnostic_plots(mcmc_nngp_list, starting_proportion = .4, burn_in = .4)
  dev.off()
  
  #############################
  # Predictions of the bousin #
  #############################
  # This part predicts the process at unobserved locations. 
  # It will be used for prediction error. 
  predicted_locs = rbind(lump_center_locs, loo_locs)
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
    burn_in = .4, 
    # number of OMP process per socket
    num_threads_per_chain = ceiling(parallel::detectCores()/(2)), 
    parallel = T, lib.loc = "../R_packages"
    )
  
  # Fixed effects 
  preds_fixed = Bidart::predict_fixed_effects(
    mcmc_nngp_list = mcmc_nngp_list, 
    # covariates at new locations
    X_pred = as.data.frame(predicted_locs), 
    burn_in = .4, lib.loc = "../R_packages"
  )
    
  # Noise 
  preds_noise = Bidart::predict_noise(
    mcmc_nngp_list = mcmc_nngp_list, 
    # covariates at new locations (none)
    X_noise_pred = NULL, 
    # new locations (used for GP prior)
    predicted_locs =  predicted_locs, 
    burn_in = .4)
  
  # Estimate the parameters, including the latent field at observed locations
  estimated_parameters = Bidart::estimate_parameters(
    mcmc_nngp_list = mcmc_nngp_list, 
    burn_in = .4, 
    get_samples = T, lib.loc = "../R_packages"
  )
  
  # "Predict" the fixed effects at observed locations
  estimated_fixed = Bidart::predict_fixed_effects(
    mcmc_nngp_list = mcmc_nngp_list, 
    # covariates at new locations
    X_pred = mcmc_nngp_list$data$covariates$X$arg, lib.loc = "../R_packages",
    burn_in = .4)
  
  # "Predict" the noise variance at observed locations 
  estimation_noise = Bidart::predict_noise(
    mcmc_nngp_list = mcmc_nngp_list, 
    # covariates at new locations (none)
    X_noise_pred = mcmc_nngp_list$data$covariates$noise_X$arg, 
    # new locations (used for GP prior)
    predicted_locs =  mcmc_nngp_list$data$observed_locs, 
    burn_in = .4)
  
  
  ##########
  # Scores # 
  ##########
  
  scores = c(
    "train" = Bidart::log_score_Gaussian(
      observed_field = train_field,
      latent_field_samples = estimated_parameters$samples$field_at_observed_locs, 
      log_noise_samples = estimation_noise$predicted_samples,
      fixed_effects_samples =  estimated_fixed$predicted_samples
    )$total / nrow(train_locs), 
    "loo" = Bidart::log_score_Gaussian(
      observed_field = response_variable[loo_idx_expanded],
      latent_field_samples = preds_field$predicted_samples$field[-seq(nrow(lump_center_locs)),,, drop=F], 
      log_noise_samples = preds_noise$predicted_samples[-seq(nrow(lump_center_locs)),,, drop=F],
      fixed_effects_samples =  preds_fixed$predicted_samples[-seq(nrow(lump_center_locs)),,, drop=F]
    )$total / length(loo_idx_expanded), 
    "lump" = Bidart::log_score_Gaussian(
      observed_field = response_variable[lump_center_idx_expanded],
      latent_field_samples = preds_field$predicted_samples$field[seq(nrow(lump_center_locs)),,, drop=F], 
      log_noise_samples = preds_noise$predicted_samples[seq(nrow(lump_center_locs)),,, drop=F],
      fixed_effects_samples =  preds_fixed$predicted_samples[seq(nrow(lump_center_locs)),,, drop=F]
    )$total / length(lump_center_idx_expanded), 
    "time" = Sys.time() - t1
  )

  saveRDS(scores, paste(runs_design[i_model,1], "_score", ".RDS", sep=""))
  pdf(paste(runs_design[i_model,1], ".pdf", sep=""))
  Bidart::diagnostic_plots(mcmc_nngp_list, burn_in = .4, starting_proportion = .4)
  dev.off()
  saveRDS(
    object = 
      list(  
        mcmc_nngp_list,
        preds_field,
        preds_fixed,
        preds_noise,
        estimated_parameters,
        estimated_fixed,
        estimation_noise
      ), 
    file = paste(runs_design[i_model,1], ".RDS", sep = "")
  )
  }
}

