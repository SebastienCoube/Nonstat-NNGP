

# full model wins ! 

remove(list = ls());gc()
run = readRDS("Heavy_metals/run_nsr")
Bidart::DIC(run)
run = readRDS("Heavy_metals/run_nsr_basis")
Bidart::DIC(run)
remove(list = ls());gc()
run = readRDS("Heavy_metals/run_ns")
Bidart::DIC(run)
run = readRDS("Heavy_metals/run_stat")
Bidart::DIC(run)


# getting US map
remove(list = ls());gc()
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
#gridmaps = raster::raster(gridmaps) %>% raster::aggregate(2) %>% as("SpatialGridDataFrame")
for(i in grid.list[-1]) {
  #gridmaps@data[sub(".asc", "", i[1])] <- (rgdal::readGDAL(paste("Heavy_metals/usgrids5km/",i, sep = "")) %>% raster::raster() %>% raster::aggregate(2) %>% as("SpatialGridDataFrame"))$band1
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
# adding spatial basis
KL_decomposition = readRDS("Heavy_metals/KL_decomposition.RDS") # loading karhunen loeve decomposition
locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs # getting PP basis
locs_ = unique(locs) 
set.seed(1)
locs_ = locs_[GpGp::order_maxmin(locs_),]
locs_ = rbind(locs_, predicted_coords)
NNarray = GpGp::find_ordered_nn(locs_, 5)
PP_basis = Matrix::solve(
  Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], 
                       j = NNarray[!is.na(NNarray)], 
                       x = GpGp::vecchia_Linv(c(1, 3, 1, 0), "matern_isotropic", locs_, NNarray)[!is.na(NNarray)], 
                       triangular = T
  ), 
  diag(1, nrow(locs_), 1000)
)
KL_basis = (PP_basis %*% KL_decomposition$v[,seq(20)] %*% diag(1/KL_decomposition$d[seq(20)]))[,seq(20)]
rm(PP_basis); rm(KL_decomposition); gc()

# normalizing
X_mean = readRDS("Heavy_metals/processed_data.RDS")$X_locs_mean
X_sd = readRDS("Heavy_metals/processed_data.RDS")$X_locs_sd
for(i in seq(ncol(X_pred)))
{
  X_pred[,i] = X_pred[,i]-X_mean[i]
  X_pred[,i] = X_pred[,i]/X_sd[i]
}
apply(X_pred, 2, mean)
apply(X_pred, 2, sd)

#NNarray = GpGp::find_ordered_nn(predicted_coords, 5)
#estimates_nsr = readRDS("Heavy_metals/estimates_nsr.RDS")
#t1 = Sys.time()
#test = cbind(1, X_pred) %*% estimates_nsr$range_beta[1,,]
#test = Bidart::compute_sparse_chol(covfun_name = "nonstationary_exponential_isotropic", range_beta = estimates_nsr$range_beta[1,,], NNarray = NNarray, locs = predicted_coords, range_X = cbind(1, X_pred), compute_derivative = F)
#Sys.time()-t1

## predicting nonstat model
run = readRDS("Heavy_metals/run_nsr_basis")

par(mfrow = c(2, 1))
Bidart::plot_pointillist_painting(run$data$locs, run$data$covariates$range_X$X_locs[,13])
Bidart::plot_pointillist_painting(locs_, KL_basis[,1])
Bidart::plot_pointillist_painting(run$data$locs, run$data$covariates$range_X$X_locs[,14])
Bidart::plot_pointillist_painting(locs_, KL_basis[,2])



prediction = Bidart::predict_latent_field(mcmc_nngp_list = run, predicted_locs = predicted_coords, X_range_pred = X_pred, X_scale_pred = X_pred, n_cores = 3)
saveRDS(prediction, "Heavy_metals/prediction_field_nsr_basis.RDS")
estimates = estimate_parameters(run)
saveRDS(estimates, "Heavy_metals/estimates_nsr_basis.RDS")
rm(prediction) ; gc()

## predicting stat model
run = readRDS("Heavy_metals/run_stat")
prediction = predict_latent_field(mcmc_nngp_list = run, predicted_locs = predicted_coords, X_range_pred = NULL, X_scale_pred = NULL, n_cores = 3)
saveRDS(prediction, "Heavy_metals/prediction_field_stat.RDS")

# loading predictions
prediction_stat = readRDS("Heavy_metals/prediction_field_stat.RDS")
prediction_nsr = readRDS("Heavy_metals/prediction_field_nsr.RDS")
prediction_nsr$predicted_samples = NULL
prediction_nsr$predicted_samples = NULL
estimates_nsr = readRDS("Heavy_metals/estimates_nsr.RDS")
gc()

# adding stat estimates
for(name in names(prediction_stat$summaries))
{
  gridmaps[[paste(name, "stat", sep = "_")]] = NA
  gridmaps[[paste(name, "stat", sep = "_")]][non_na_indices] = prediction_stat$summaries[[name]][1,,]
}
gridmaps[[paste("field_stat_sd")]] = NA
gridmaps[[paste("field_stat_sd")]][non_na_indices] = prediction_stat$summaries[["field"]][5,,]
# adding nonstat estimates
for(name in names(prediction_nsr$summaries))
{
  gridmaps[[paste(name, "nsr", sep = "_")]] = NA
  gridmaps[[paste(name, "nsr", sep = "_")]][non_na_indices] = prediction_nsr$summaries[[name]][1,,]
}
gridmaps[[paste("field_nsr_sd")]] = NA
gridmaps[[paste("field_nsr_sd")]][non_na_indices] = prediction_nsr$summaries[["field"]][5,,]
for(name in colnames(X_pred))
{
  name_ = paste("range", name, sep = "_")
  gridmaps[[name_]] = NA
  gridmaps[[name_]][non_na_indices] = X_pred[,name] * estimates_nsr$range_beta[1,name,]
  name_ = paste("scale", name, sep = "_")
  gridmaps[[name_]] = NA
  gridmaps[[name_]][non_na_indices] = X_pred[,name] * estimates_nsr$scale_beta[1,name,]
  name_ = paste("noise", name, sep = "_")
  gridmaps[[name_]] = NA
  gridmaps[[name_]][non_na_indices] = X_pred[,name] * estimates_nsr$noise_beta[1,name,]
}



# plotting mean
pdf("Heavy_metals/field_mean.pdf", width = 13, height = 4)
zlim = c(-2, 4)
layout(matrix(1:3, 1, 3), widths = c(5,5,1))
sp::plot(gridmaps["field_nsr"], what = "image", zlim = zlim, main = "nonstationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat"], what = "image", zlim = zlim, main = "stationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat"], what = "scale", zlim = zlim )
dev.off()
# plotting mean
zlim = c(-2, 4)
#layout(matrix(1:3, 1, 3), widths = c(5,5,1))
pdf("Heavy_metals/field_mean_nsr.pdf", width = 5, height = 5)
sp::plot(gridmaps["field_nsr"], what = "image", zlim = zlim, main = "nonstationary" )
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
sp::plot(gridmaps["field_nsr_sd"], what = "image", zlim = zlim, main = "nonstationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat_sd"], what = "image", zlim = zlim, main = "stationary" )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["field_stat_sd"], what = "scale", zlim = zlim )
dev.off()


# PCA of parameters
covparms_beta = 
  cbind(
    estimates_nsr$range_beta[1, -1,],
    estimates_nsr$scale_beta[1, -1,],
    estimates_nsr$noise_beta[1, -1,]
  )
colnames(covparms_beta) = c("range", "scale", "noise")
row.names(covparms_beta) = c(
  "air_pollution_dens", 
  "mining_dens", 
  "earthquake_dens", 
  "toxic_release_dens", 
  "green_biomass", 
  "elevation", 
  "night_lights_dens", 
  "roads_density", 
  "wetness", 
  "visible_sky", 
  "wind_effect" 
)
ACP = FactoMineR::PCA(covparms_beta)
pdf("Heavy_metals/pca_variables.pdf")
FactoMineR::plot.PCA(ACP, choix = "var")
dev.off()
pdf("Heavy_metals/pca_individuals.pdf", width = 5, height = 5)
FactoMineR::plot.PCA(ACP, choix = "ind")
dev.off()



xtable::xtable(digits = 3,
  cbind(
    t(estimates_nsr$range_beta[c(1, 2, 4), ,]),
    t(estimates_nsr$scale_beta[c(1, 2, 4), ,]),
    t(estimates_nsr$noise_beta[c(1, 2, 4), ,])
  ) 
)


# plotting PCA first component
gridmaps[["field_incoherence"]] = NA
gridmaps[["field_incoherence"]][non_na_indices] = c(X_pred%*%ACP$ind$coord[,1])
pdf("Heavy_metals/PC1.pdf", height = 5, width = 10)
sp::plot(gridmaps["field_incoherence"][,], main = "", )
sp::plot(usa_aea, add = T)
dev.off()



# plotting PCA second component
gridmaps[["range_noise"]] = NA
gridmaps[["range_noise"]][non_na_indices] = c(X_pred%*%ACP$ind$coord[,2])
pdf("Heavy_metals/PC2.pdf", height = 5, width = 10)
sp::plot(gridmaps["range_noise"], main = "")
sp::plot(usa_aea, add = T)
dev.off()






# range
gridmaps[["range"]] = NA
gridmaps[["range"]][non_na_indices] = cbind(1, X_pred) %*% estimates_nsr$range_beta[1,,]
pdf("Heavy_metals/range.pdf", height = 5, width = 10)
sp::plot(gridmaps["range"], main = "")
sp::plot(usa_aea, add = T)
dev.off()
# scale
gridmaps[["scale"]] = NA
gridmaps[["scale"]][non_na_indices] = cbind(1, X_pred) %*% estimates_nsr$scale_beta[1,,]
pdf("Heavy_metals/scale.pdf", height = 5, width = 10)
sp::plot(gridmaps["scale"], main = "")
sp::plot(usa_aea, add = T)
dev.off()
# noise
gridmaps[["noise"]] = NA
gridmaps[["noise"]][non_na_indices] = cbind(1, X_pred) %*% estimates_nsr$noise_beta[1,,]
pdf("Heavy_metals/noise.pdf", height = 5, width = 10)
sp::plot(gridmaps["noise"], main = "")
sp::plot(usa_aea, add = T)
dev.off()




sp::plot(gridmaps["log_range"], main = "", )
sp::plot(gridmaps["log_scale"], main = "", )
sp::plot(gridmaps["field"], main = "", )
sp::plot(gridmaps["field"], main = "", )



sp::plot(gridmaps["range_twi"], main = "", )
sp::plot(gridmaps["range_globedem"], main = "", )
sp::plot(gridmaps["range_gcarb"], main = "", )
sp::plot(gridmaps["range_nlights03"], main = "", )
sp::plot(gridmaps["range_sdroads"], main = "", )
sp::plot(gridmaps["range_dairp"], main = "", )
sp::plot(gridmaps["range_dmino"], main = "", )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["noise_sdroads"], main = "", )
sp::plot(usa_aea, add = T)


sp::plot(gridmaps["field"], main = "", )
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["log_range"], main = "", )
sp::plot(usa_aea, add = T)


run = readRDS("Heavy_metals/run_nsr")

samples = do.call(rbind, 
  lapply(
    names(run$records$chain_1)[grep("beta", names(run$records$chain_1))], 
    function(name)
    {
      res = do.call(abind::abind, lapply(run$records, function(x)x[[name]][,,-seq(100)]))
      row.names(res) = paste(name, 
                             c(
                               "(Intercept)",
                               "air_pollution_dens", 
                               "mining_dens", 
                               "earthquake_dens", 
                               "toxic_release_dens", 
                               "green_biomass", 
                               "elevation", 
                               "night_lights_dens", 
                               "roads_density", 
                               "wetness", 
                               "visible_sky", 
                               "wind_effect" 
                             )
                             , sep = "_")
      res
    }
  ))
colnames(samples) = NULL

pdf("Heavy_metals/corrplot.pdf", width = 10, height = 10)
corrplot::corrplot(cor(as.data.frame(t(samples))))
dev.off()
behp = mnt::test.BHEP(t(samples), MC.rep = 500)
ACP = FactoMineR::PCA(t(samples), ncp = 11)
qqnorm(ACP$ind$coord[,11])
dev.off()



gridmaps[["field_diff_stat_nonstat"]] =NA 
gridmaps[["field_diff_stat_nonstat"]] =gridmaps[["field_nsr"]] - gridmaps[["field_stat"]]

sp::plot(gridmaps["field_diff_stat_nonstat"])
sp::plot(usa_aea, add = T)
