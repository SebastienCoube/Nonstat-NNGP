

# getting US map
#remove(list = ls());gc()
library(magrittr)
usa <- maps::map("state", fill = T, plot = F)
IDs <- sapply(strsplit(usa$names, ":"), function(x) x[1])
usa <- maptools::map2SpatialPolygons(usa, IDs=IDs, proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
# getting 5Km spatial grid and lowering resolution to 10 Km
grid.list <- c("dairp.asc", "dmino.asc", "dquksig.asc", "dTRI.asc", "gcarb.asc",
               "geomap.asc", "globedem.asc", "minotype.asc", "nlights03.asc", "sdroads.asc",
               "twi.asc", "vsky.asc", "winde.asc", "glwd31.asc")
gridmaps <- rgdal::readGDAL("Heavy_metals/usgrids5km/dairp.asc")
names(gridmaps)[1] <- sub(".asc", "", grid.list[1])
for(i in grid.list[-1]) {
  gridmaps@data[sub(".asc", "", i[1])] <- rgdal::readGDAL(paste("Heavy_metals/usgrids5km/",i, sep = ""))$band1
}
# projection coords
AEA <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0
+ellps=GRS80 +datum=NAD83 +units=m +no_defs"
sp::proj4string(gridmaps) = sp::CRS(AEA)
# overlaying map on us map
usa_aea = sp::spTransform(usa, sp::CRS(AEA))
gridmaps_overlay = sp::over(gridmaps, usa_aea)
# transforming gridmaps into WGS84 data
gridmaps_WGS84 = sp::spTransform(gridmaps, sp::CRS("+proj=longlat +datum=WGS84"))
predicted_coords = sp::coordinates(gridmaps_WGS84)
# getting X_pred
X_pred = as.matrix(as.data.frame((gridmaps_WGS84)))
X_pred = X_pred[,-c(6, 8, 14, 15, 16)]
# getting na indices
non_na_indices = (!is.na(gridmaps_overlay))&(!apply(X_pred, 1, anyNA))
# removing missing values in X_pred
predicted_coords = predicted_coords[non_na_indices,]
X_pred = X_pred[non_na_indices,]
# normalizing X_pred
X_mean = readRDS("Heavy_metals/processed_data.RDS")$X_locs_mean
X_sd = readRDS("Heavy_metals/processed_data.RDS")$X_locs_sd
for(i in seq(ncol(X_pred)))
{
  X_pred[,i] = X_pred[,i]-X_mean[i]
  X_pred[,i] = X_pred[,i]/X_sd[i]
}
apply(X_pred, 2, mean)
apply(X_pred, 2, sd)
# adding spatial basis
train_data_set = readRDS("Heavy_metals/validation_train.RDS")
PP = train_data_set$train_PP
# rounding
round_idx = !duplicated(round(predicted_coords*5, 0)/5)
round_idx_rev = match(split(round(predicted_coords*5, 0)/5, row(predicted_coords)), split((round(predicted_coords*5, 0)/5)[round_idx,], row(predicted_coords[round_idx,])))
predicted_coords = predicted_coords[round_idx,]
X_pred = X_pred[round_idx,]

# predicting and estimating nonstat model
run = readRDS("Heavy_metals/run_noise_TRUE_scale_TRUE_.RDS")
prediction = Bidart::predict_latent_field(mcmc_nngp_list = run, predicted_locs = predicted_coords, X_scale_pred = X_pred[, c("gcarb", "globedem", "twi")], num_threads_per_chain = 10, parallel = T, burn_in = .15)
saveRDS(prediction, "Heavy_metals/prediction_field_nonstat.RDS")
estimates = Bidart::estimate_parameters(run, burn_in = .15)
saveRDS(estimates, "Heavy_metals/estimates_nonstat.RDS")
noise_prediction = Bidart::predict_noise(mcmc_nngp_list = run, predicted_locs = predicted_coords, X_noise_pred =  X_pred, burn_in = .15)
saveRDS(noise_prediction, "Heavy_metals/prediction_noise.RDS")
rm(prediction) ; gc()
# predicting stat model
run = readRDS("Heavy_metals/run_noise_FALSE_scale_FALSE_.RDS")
prediction = Bidart::predict_latent_field(mcmc_nngp_list = run, predicted_locs = predicted_coords, X_range_pred = NULL, X_scale_pred = NULL, num_threads_per_chain = 10, parallel = T, burn_in = .15)
saveRDS(prediction, "Heavy_metals/prediction_field_stat.RDS")
estimates = Bidart::estimate_parameters(run, burn_in = .15)
saveRDS(estimates, "Heavy_metals/estimates_stat.RDS")

# loading predictions
prediction_ns = readRDS("Heavy_metals/prediction_field_nonstat.RDS")
prediction_stat = readRDS("Heavy_metals/prediction_field_stat.RDS")
prediction_ns$predicted_samples = NULL
prediction_stat$predicted_samples = NULL
estimates_ns = readRDS("Heavy_metals/estimates_nonstat.RDS")
estimates_stat = readRDS("Heavy_metals/estimates_stat.RDS")
noise_prediction = readRDS("Heavy_metals/prediction_noise.RDS")
##gc()

# adding stat estimates
for(name in names(prediction_stat$summaries))
{
  gridmaps[[paste(name, "stat_mean", sep = "_")]] = NA
  gridmaps[[paste(name, "stat_mean", sep = "_")]][non_na_indices] = prediction_stat$summaries[[name]][1,round_idx_rev,]
}
gridmaps[[paste("field_stat_sd")]] = NA
gridmaps[[paste("field_stat_sd")]][non_na_indices] = prediction_stat$summaries[["field"]][5,round_idx_rev,]
# adding stat estimates
for(name in names(prediction_stat$summaries))
{
  gridmaps[[paste(name, "ns_mean", sep = "_")]] = NA
  gridmaps[[paste(name, "ns_mean", sep = "_")]][non_na_indices] = prediction_ns$summaries[[name]][1,round_idx_rev,]
}
gridmaps[[paste("field_ns_sd")]] = NA
gridmaps[[paste("field_ns_sd")]][non_na_indices] = prediction_ns$summaries[["field"]][5,round_idx_rev,]


# plotting mean
pdf("Heavy_metals/field_mean.pdf", width = 13, height = 4)
zlim = c(min(gridmaps$field_ns_mean, gridmaps$field_stat_mean, na.rm = T), max(gridmaps$field_ns_mean, gridmaps$field_stat_mean, na.rm = T))
layout(matrix(1:3, 1, 3), widths = c(5,5,1))
sp::plot(gridmaps["field_ns_mean"], what = "image", zlim = zlim, main = "nonstationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat_mean"], what = "image", zlim = zlim, main = "stationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat_mean"], what = "scale", zlim = zlim )
dev.off()
# plotting mean
layout(matrix(1:3, 1, 3), widths = c(5,5,1))
zlim = c(min(gridmaps$field_ns_mean, gridmaps$field_stat_mean, na.rm = T), max(gridmaps$field_ns_mean, gridmaps$field_stat_mean, na.rm = T))
pdf("Heavy_metals/field_mean_ns.pdf", width = 5, height = 5)
sp::plot(gridmaps["field_ns_mean"], what = "image", zlim = zlim, main = "nonstationary" )
sp::plot(usa_aea, add = T)
dev.off()
pdf("Heavy_metals/field_mean_stat.pdf", width = 5, height = 5)
sp::plot(gridmaps["field_stat_mean"], what = "image", zlim = zlim, main = "stationary" )
sp::plot(usa_aea, add = T)
dev.off()
pdf("Heavy_metals/field_mean_scale.pdf", width = 2, height = 5)
sp::plot(gridmaps["field_stat_mean"], what = "scale", zlim = zlim )
dev.off()


# plotting sd
pdf("Heavy_metals/field_sd.pdf", width = 13, height = 4)
zlim = c(min(gridmaps$field_ns_sd, gridmaps$field_stat_sd, na.rm = T), max(gridmaps$field_ns_sd, gridmaps$field_stat_sd, na.rm = T))
layout(matrix(1:3, 1, 3), widths = c(5,5,1), heights = c(5,5,5))
sp::plot(gridmaps["field_ns_sd"], what = "image", zlim = zlim, main = "nonstationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat_sd"], what = "image", zlim = zlim, main = "stationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat_sd"], what = "scale", zlim = zlim )
dev.off()


# scale
pdf("Heavy_metals/scale.pdf", height = 5, width = 10)
sp::plot(gridmaps["log_scale_ns_mean"], main = "")
sp::plot(usa_aea, add = T)
dev.off()
pdf("Heavy_metals/scale_title.pdf", height = 5, width = 10)
sp::plot(gridmaps["log_scale_ns_mean"], main = "Mean log-variance of the latent field")
sp::plot(usa_aea, add = T)
dev.off()





# noise 
gridmaps[["noise"]][non_na_indices] = noise_prediction$summaries[1,round_idx_rev,]
pdf("Heavy_metals/noise.pdf", height = 5, width = 10)
sp::plot(gridmaps["noise"], main = "")
sp::plot(usa_aea, add = T)
dev.off()
gridmaps[["noise"]][non_na_indices] = noise_prediction$summaries[1,round_idx_rev,]
pdf("Heavy_metals/noise_title.pdf", height = 5, width = 10)
sp::plot(gridmaps["noise"], main = "Mean log-variance of the noise")
sp::plot(usa_aea, add = T)
dev.off()

