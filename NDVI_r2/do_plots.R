#setwd("NDVI_r2")
#mcmc_nngp_list = readRDS("res/NNGP_full_monty.RDS")

.libPaths("../R_packages/")
mcmc_nngp_list = readRDS("NNGP_full_monty.RDS")


mcmc_nngp_list = mcmc_nngp_list[[1]]

##########################
# getting pred locations #
##########################

med_dist = median(FNN::knn.dist(data = mcmc_nngp_list$data$locs, k=1))
predicted_locs = as.matrix(
expand.grid(
  seq(
    min(mcmc_nngp_list$data$locs[,1])-.1,
    max(mcmc_nngp_list$data$locs[,1])+.1, 
    med_dist
  ), seq(
    min(mcmc_nngp_list$data$locs[,2])-.1,
    max(mcmc_nngp_list$data$locs[,2])+.1, 
    med_dist
  )
)
)
dist_predicted_locs = FNN::knnx.dist(data = mcmc_nngp_list$data$locs, query = predicted_locs, k=1)
predicted_locs = predicted_locs[which(dist_predicted_locs>(med_dist)),]
dist_predicted_locs = FNN::knnx.dist(data = mcmc_nngp_list$data$locs, query = predicted_locs, k=1)
predicted_locs = predicted_locs[which(dist_predicted_locs<.1),]
plot(predicted_locs, cex = .05, pch = 16)




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
 burn_in = .4, lib.loc = "../R_packages/",
 # number of OMP process per socket
 num_threads_per_chain = 10)
saveRDS(preds_field, "maps/preds_field.RDS")
# Fixed effects 
preds_fixed = Bidart::predict_fixed_effects(
 mcmc_nngp_list = mcmc_nngp_list, 
 # covariates at new locations
 X_pred = as.data.frame(predicted_locs),  lib.loc = "../R_packages/",
 burn_in = .4)
saveRDS(preds_fixed, "maps/preds_fixed.RDS")

# Noise 
preds_noise = Bidart::predict_noise(
 mcmc_nngp_list = mcmc_nngp_list, 
 # covariates at new locations (none)
 X_noise_pred = NULL, 
 # new locations (used for GP prior)
 predicted_locs =  predicted_locs,
 burn_in = .4)
saveRDS(preds_noise, "maps/preds_noise.RDS")

##############
# Estimation #
##############
# This part estimates the process at observed locations. 
# It will be used for estimation/training error. 

# Estimate the parameters, including the latent field at observed locations
estimated_parameters = Bidart::estimate_parameters(
 mcmc_nngp_list = mcmc_nngp_list,  lib.loc = "../R_packages/",
 burn_in = .4, 
 get_samples = T
)
saveRDS(estimated_parameters, "maps/estimated_parameters.RDS")

# "Predict" the fixed effects at observed locations
estimated_fixed = Bidart::predict_fixed_effects(
 mcmc_nngp_list = mcmc_nngp_list,  lib.loc = "../R_packages/",
 # covariates at new locations
 X_pred = mcmc_nngp_list$data$covariates$X$arg, 
 burn_in = .1)
saveRDS(estimated_fixed, "maps/estimated_fixed.RDS")

# "Predict" the noise variance at observed locations 
estimation_noise = Bidart::predict_noise(
 mcmc_nngp_list = mcmc_nngp_list, 
 # covariates at new locations (none)
 X_noise_pred = mcmc_nngp_list$data$covariates$noise_X$arg, 
 # new locations (used for GP prior)
 predicted_locs =  mcmc_nngp_list$data$observed_locs,
 burn_in = .1)
saveRDS(estimation_noise, "maps/estimation_noise.RDS")


############
# Ellipses #
############
locs_4_ellipses = as.matrix(expand.grid(seq(-.2, 1.3, .05), seq(-.2, 1.3, .05)))
locs_4_ellipses = locs_4_ellipses + rnorm(1, 0, .00001)
dist_locs_4_ellipses = FNN::knnx.dist(data = mcmc_nngp_list$data$locs, query = locs_4_ellipses, k=1)
locs_4_ellipses = locs_4_ellipses[which(dist_locs_4_ellipses<.5),]
ellipses_pred = Bidart::predict_latent_field(
  mcmc_nngp_list = mcmc_nngp_list, 
  predicted_locs = locs_4_ellipses,  
  lib.loc = "../R_packages/",
  burn_in = .4, num_threads_per_chain = 10
)
saveRDS(ellipses_pred, "maps/ellipses_pred.RDS")


setwd("NDVI_r2")
############################
# plotting range log scale #
############################
pdf("range_var.pdf")
Bidart::plot_log_scale(
  log_scale_arrays = list(
     mcmc_nngp_list$records$chain_1$range_log_scale, 
     mcmc_nngp_list$records$chain_2$range_log_scale
     ), 
  iterations = mcmc_nngp_list$iterations$thinning, 
  starting_proportion = .4, varname = "Variance matrix of the range PP"
)
dev.off()

png("range_var.png")
Bidart::plot_log_scale(
  log_scale_arrays = list(
    mcmc_nngp_list$records$chain_1$range_log_scale, 
    mcmc_nngp_list$records$chain_2$range_log_scale
  ), 
  iterations = mcmc_nngp_list$iterations$thinning, 
  starting_proportion = .4, varname = "Variance matrix of the range PP"
)
dev.off()

############
# Plotting #
############

setwd("NDVI_r2/")
mcmc_nngp_list = readRDS("res/NNGP_full_monty.RDS")
mcmc_nngp_list = mcmc_nngp_list[[1]]

estimated_fixed = readRDS("res/maps/estimated_fixed.RDS")
estimated_parameters = readRDS("res/maps/estimated_parameters.RDS")
estimation_noise = readRDS("res/maps/estimation_noise.RDS")
preds_field = readRDS("res/maps/preds_field.RDS")
preds_fixed = readRDS("res/maps/preds_fixed.RDS")
preds_noise = readRDS("res/maps/preds_noise.RDS")
ellipse_preds = readRDS("res/maps/ellipses_pred.RDS")
PP = mcmc_nngp_list$hierarchical_model$PP
gc()


png("maps/observed.png", width = 4000, height = 4000)
par(mar = c(0,0,0,0))
par(bg="blue")
plot(0,0, type= "n", xlim = c(-.1, 1.2), ylim = c(-.1, 1.2), bty="n", 
     xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", 
)
Bidart::plot_pointillist_painting(
  mcmc_nngp_list$data$observed_locs,
  mcmc_nngp_list$data$observed_field, 
  cex = 3, add = T
)
dev.off()


png("maps/noise.png", width = 4000, height = 4000)
par(mar = c(0,0,0,0))
par(bg="white")
plot(0,0, type= "n", xlim = c(-.1, 1.2), ylim = c(-.1, 1.2), bty="n", 
     xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", 
)
Bidart::plot_pointillist_painting(
  rbind(mcmc_nngp_list$data$observed_locs, preds_noise$predicted_locs),
  c(estimation_noise$summaries["mean",,], preds_noise$summaries["mean",,]), 
  cex = 3, add = T
)
dev.off()


png("maps/log_sd.png", width = 4000, height = 4000)
par(mar = c(0,0,0,0))
par(bg="white")
plot(0,0, type= "n", xlim = c(-.1, 1.2), ylim = c(-.1, 1.2), bty="n", 
     xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", 
)
Bidart::plot_pointillist_painting(
  rbind(mcmc_nngp_list$data$locs, preds_noise$predicted_locs),
  log(c(estimated_parameters$summaries$field["sd",,], preds_field$summaries$field["sd",,])), 
  cex = 3, add = T
)
dev.off()



png("maps/ellipses.png", width = 4000, height = 4000)
par(mar = c(0,0,0,0))
par(bg="white", lwd = 10)
plot(0,0, type= "n", xlim = c(-.1, 1.2), ylim = c(-.1, 1.2), bty="n", 
     xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", 
)
Bidart::plot_pointillist_painting(
  rbind(mcmc_nngp_list$data$locs, preds_noise$predicted_locs),
  (c(estimated_parameters$summaries$field["mean",,], preds_field$summaries$field["mean",,])), 
  cex = 3, add = T)

locs_4_ellipses = as.matrix(expand.grid(seq(-.2, 1.3, .05), seq(-.2, 1.3, .05)))
dist_locs_4_ellipses = FNN::knnx.dist(data = mcmc_nngp_list$data$locs, query = locs_4_ellipses, k=1)
locs_4_ellipses = locs_4_ellipses[which(dist_locs_4_ellipses<.5),]
dist_locs_4_ellipses = FNN::knnx.dist(data = mcmc_nngp_list$data$locs, query = locs_4_ellipses, k=1)
ellipse_idx = which((dist_locs_4_ellipses<.03)&(apply((100*locs_4_ellipses)%%10, 1, sum)<1))

plot_ellipses(
  locs_4_ellipses[ellipse_idx,], 
  ellipse_preds$summaries$log_range["mean",ellipse_idx,],
  add=T, shrink = sqrt(8*mcmc_nngp_list$hierarchical_model$nu)
)
dev.off()

