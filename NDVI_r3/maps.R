############
# Plotting #
############
setwd("NDVI_r3/")

setwd("res")
mcmc_nngp_list = readRDS("NNGP_full_monty.RDS")
mcmc_nngp_list = mcmc_nngp_list[[1]]

estimated_fixed = readRDS("maps/estimated_fixed.RDS")
estimated_parameters = readRDS("maps/estimated_parameters.RDS")
estimation_noise = readRDS("maps/estimation_noise.RDS")
preds_field = readRDS("maps/preds_field.RDS")
preds_fixed = readRDS("maps/preds_fixed.RDS")
preds_noise = readRDS("maps/preds_noise.RDS")
ellipse_preds = readRDS("maps/ellipses_pred.RDS")
PP = mcmc_nngp_list$hierarchical_model$PP
gc()
setwd("../maps_for_JCGS")

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




png("observed.png", width = 4000, height = 4000)
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

png("observed_scale.png", width = 100, height = 400)
Bidart::pointillist_colorscale(mcmc_nngp_list$data$observed_field)
dev.off()

png("noise.png", width = 4000, height = 4000)
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

png("noise_scale.png", width = 100, height = 400)
par(mar = c(0,0,0,0))
Bidart::pointillist_colorscale(
  c(estimation_noise$summaries["mean",,], preds_noise$summaries["mean",,])
)
dev.off()


png("log_sd.png", width = 4000, height = 4000)
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
png("log_sd_scale.png", width = 100, height = 400)
  Bidart::pointillist_colorscale(log(c(estimated_parameters$summaries$field["sd",,], preds_field$summaries$field["sd",,])))
dev.off()



png("ellipses.png", width = 4000, height = 4000)
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

Bidart::plot_ellipses(
  locs_4_ellipses[ellipse_idx,], 
  ellipse_preds$summaries$log_range["mean",ellipse_idx,],
  add=T, shrink = sqrt(8*mcmc_nngp_list$hierarchical_model$nu)
)
dev.off()



png("mean_scale.png", width = 100, height = 400)
Bidart::pointillist_colorscale((c(estimated_parameters$summaries$field["mean",,], preds_field$summaries$field["mean",,])))
dev.off()


mean_log_range =  Bidart::X_PP_mult_right(
  X = mcmc_nngp_list$data$covariates$range_X$X_locs, 
  PP = mcmc_nngp_list$hierarchical_model$PP, use_PP = T, 
  locs_idx = mcmc_nngp_list$vecchia_approx$hctam_scol_1, Y = estimated_parameters$summaries$range_beta[1,,]
)[,1]
Bidart::plot_pointillist_painting(mcmc_nngp_list$data$locs, mean_log_range)
