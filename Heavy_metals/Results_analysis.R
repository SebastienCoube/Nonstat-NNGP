



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
    # KL_decomposition = readRDS("Heavy_metals/KL_decomposition.RDS") # loading karhunen loeve decomposition
    # KL_basis = Bidart::predict_KL_basis(predicted_locs = predicted_coords, KL_basis = KL_decomposition)
    # saveRDS(KL_basis, "KL_basis_pred.RDS")
KL_basis = readRDS("KL_basis_pred.RDS")
X_range_pred = as.matrix(cbind(KL_basis, X_pred[, c("gcarb", "globedem", "twi")]))
X_pred_basis = as.matrix(cbind(KL_basis, X_pred));colnames(X_pred_basis)[seq(25)] = as.character(seq(25))


### predicting and estimating nonstat model
###run = readRDS("Heavy_metals/run_nr")
###source("Bidart/R/mcmc_nngp_predict_nonstationary.R")
###prediction = Bidart::predict_latent_field(mcmc_nngp_list = run, predicted_locs = predicted_coords, X_range_pred = X_range_pred, n_cores = 3)
###saveRDS(prediction, "Heavy_metals/prediction_field_nr_basis.RDS")
###estimates = Bidart::estimate_parameters(readRDS("Heavy_metals/run_nr"))
###saveRDS(estimates, "Heavy_metals/estimates_nr_basis.RDS")
###rm(prediction) ; gc()
###
##### predicting stat model
###run = readRDS("Heavy_metals/run_stat")
###prediction = Bidart::predict_latent_field(mcmc_nngp_list = run, predicted_locs = predicted_coords, X_range_pred = NULL, X_scale_pred = NULL, n_cores = 3)
###saveRDS(prediction, "Heavy_metals/prediction_field_stat.RDS")
###estimates = Bidart::estimate_parameters(run)
###saveRDS(estimates, "Heavy_metals/estimates_nr_basis.RDS")

# loading predictions
prediction_stat = readRDS("Heavy_metals/prediction_field_stat.RDS")
prediction_nr = readRDS("Heavy_metals/prediction_field_nr_basis.RDS")
prediction_nr$predicted_samples = NULL
prediction_stat$predicted_samples = NULL
estimates_nr = readRDS("Heavy_metals/estimates_nr_basis.RDS")
gc()

# adding stat estimates
for(name in names(prediction_stat$summaries))
{
  gridmaps[[paste(name, "stat", sep = "_")]] = NA
  gridmaps[[paste(name, "stat", sep = "_")]][non_na_indices] = prediction_stat$summaries[[name]][1,]
}
gridmaps[[paste("field_stat_sd")]] = NA
gridmaps[[paste("field_stat_sd")]][non_na_indices] = prediction_stat$summaries[["field"]][5,]
# adding nonstat estimates
for(name in names(prediction_nr$summaries))
{
  gridmaps[[paste(name, "nr", sep = "_")]] = NA
  gridmaps[[paste(name, "nr", sep = "_")]][non_na_indices] = prediction_nr$summaries[[name]][1,]
}
gridmaps[[paste("field_nr_sd")]] = NA
gridmaps[[paste("field_nr_sd")]][non_na_indices] = prediction_nr$summaries[["field"]][5,]



# plotting mean
pdf("Heavy_metals/field_mean.pdf", width = 13, height = 4)
zlim = c(-2, 4)
layout(matrix(1:3, 1, 3), widths = c(5,5,1))
sp::plot(gridmaps["field_nr"], what = "image", zlim = zlim, main = "nonstationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat"], what = "image", zlim = zlim, main = "stationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat"], what = "scale", zlim = zlim )
dev.off()
# plotting mean
zlim = c(-2, 4)
layout(matrix(1:3, 1, 3), widths = c(5,5,1))
pdf("Heavy_metals/field_mean_nr.pdf", width = 5, height = 5)
sp::plot(gridmaps["field_nr"], what = "image", zlim = zlim, main = "nonstationary" )
sp::plot(usa_aea, add = T)
dev.off()
pdf("Heavy_metals/field_mean_stat.pdf", width = 5, height = 5)
sp::plot(gridmaps["field_stat"], what = "image", zlim = zlim, main = "stationary" )
sp::plot(usa_aea, add = T)
dev.off()
pdf("Heavy_metals/field_mean_scale.pdf", width = 2, height = 5)
sp::plot(gridmaps["field_stat"], what = "scale", zlim = zlim )
dev.off()


# plotting sd
pdf("Heavy_metals/field_sd.pdf", width = 13, height = 4)
zlim = c(0, 1.15)
layout(matrix(1:3, 1, 3), widths = c(5,5,1), heights = c(5,5,5))
sp::plot(gridmaps["field_nr_sd"], what = "image", zlim = zlim, main = "nonstationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat_sd"], what = "image", zlim = zlim, main = "stationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat_sd"], what = "scale", zlim = zlim )
dev.off()


# range
gridmaps[["range"]] = NA
gridmaps[["range"]][non_na_indices] = cbind(1, X_range_pred) %*% estimates_nr$summaries$range_beta[1,,]
pdf("Heavy_metals/range.pdf", height = 5, width = 10)
sp::plot(gridmaps["range"], main = "")
sp::plot(usa_aea, add = T)
dev.off()

gridmaps[["range_covariates"]] = NA
gridmaps[["range_covariates"]][non_na_indices] = cbind(1, X_pred[, c("gcarb", "globedem", "twi")]) %*% estimates_nr$summaries$range_beta[1, c("(Intercept)", "gcarb", "globedem", "twi"),]
pdf("Heavy_metals/range_covariates", height = 5, width = 10)
sp::plot(gridmaps["range_covariates"], main = "")
sp::plot(usa_aea, add = T)
dev.off()

gridmaps[["range_basis"]] = NA
gridmaps[["range_basis"]][non_na_indices] = cbind(1, X_pred_basis[,seq(25)]) %*% estimates_nr$summaries$range_beta[1,seq(26),]
pdf("Heavy_metals/range_basis", height = 5, width = 10)
sp::plot(gridmaps["range_basis"], main = "")
sp::plot(usa_aea, add = T)
dev.off()


pdf("Heavy_metals/range_all.pdf", height = 4, width = 16)
layout(matrix(1:4, 1, 4), widths = c(5,5,5,1))
zlim = c(
  min(c(gridmaps$range, gridmaps$range_covariates, gridmaps$range_basis), na.rm = T), 
  max(c(gridmaps$range, gridmaps$range_covariates, gridmaps$range_basis), na.rm = T)
)
sp::plot(gridmaps["range"],  zlim = zlim, what = "image", main = "Total mean log range")
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["range_covariates"], zlim = zlim, what = "image", main = "Mean log range (covariates + intercept)")
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["range_basis"], zlim = zlim, what = "image", main = "Mean log range (basis + intercept)")
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["range_basis"], what = "scale", zlim = zlim )
dev.off()




# noise 

gridmaps[["noise"]] = NA
gridmaps[["noise"]][non_na_indices] = cbind(1, X_pred_basis) %*% estimates_nr$summaries$noise_beta[1,,]
pdf("Heavy_metals/noise.pdf", height = 5, width = 10)
sp::plot(gridmaps["noise"], main = "")
sp::plot(usa_aea, add = T)
dev.off()

gridmaps[["noise_covariates"]] = NA
gridmaps[["noise_covariates"]][non_na_indices] = cbind(1, X_pred) %*% estimates_nr$summaries$noise_beta[1,-seq(2, 26),]
pdf("Heavy_metals/noise_covariates", height = 5, width = 10)
sp::plot(gridmaps["noise_covariates"], main = "")
sp::plot(usa_aea, add = T)
dev.off()

gridmaps[["noise_basis"]] = NA
gridmaps[["noise_basis"]][non_na_indices] = cbind(1, X_pred_basis[,seq(25)]) %*% estimates_nr$summaries$noise_beta[1,seq(26),]
pdf("Heavy_metals/noise_basis", height = 5, width = 10)
sp::plot(gridmaps["noise_basis"], main = "")
sp::plot(usa_aea, add = T)
dev.off()


pdf("Heavy_metals/noise_all.pdf", height = 3, width = 12)
layout(matrix(1:4, 1, 4), widths = c(5,5,5,1))
zlim = c(
  min(c(gridmaps$noise, gridmaps$noise_covariates, gridmaps$noise_basis), na.rm = T), 
  max(c(gridmaps$noise, gridmaps$noise_covariates, gridmaps$noise_basis), na.rm = T)
)
sp::plot(gridmaps["noise"],  zlim = zlim, what = "image", main = "Total mean log noise")
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["noise_covariates"], zlim = zlim, what = "image", main = "Mean log noise (covariates + intercept)")
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["noise_basis"], zlim = zlim, what = "image", main = "Mean log noise (basis + intercept)")
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["noise_basis"], what = "scale", zlim = zlim )
dev.off()




run = readRDS("Heavy_metals/run_nr")
estimation_with_samples = Bidart::estimate_parameters(run, get_samples = T)


samples = do.call(rbind, 
                  list(
                    estimation_with_samples$samples$range_beta[,,], 
                    estimation_with_samples$samples$noise_beta[,,], 
                    estimation_with_samples$samples$scale_beta,
                    estimation_with_samples$samples$beta[,,]
                  )
                  )

row.names(samples) = c(
  sapply(colnames(estimation_with_samples$summaries$range_beta), function(x)paste("range_", x, sep = "")) ,
  sapply(colnames(estimation_with_samples$summaries$noise_beta), function(x)paste("noise_", x, sep = "")) ,
  sapply(colnames(estimation_with_samples$summaries$scale_beta), function(x)paste("scale_", x, sep = "")) ,
  sapply(colnames(estimation_with_samples$summaries$beta)      , function(x)paste("beta_", x, sep = ""))   
)

png("Heavy_metals/corrplot.png", width = 1500, height = 1500)
corrplot::corrplot(cor(as.data.frame(t(samples))))
dev.off()
behp = mnt::test.BHEP(t(samples), MC.rep = 500)
ACP = FactoMineR::PCA(t(samples), ncp = 86, graph = F)

par(mfrow = c(3, 3))
for(i in seq(86))
{
  qqnorm(ACP$ind$coord[,i])
}

par(mfrow = c(1, 1))

qqnorm(rnorm(nrow(samples))%*%samples)

