

# full model wins ! 

##DICs = c()
##run = readRDS("Heavy_metals/run_nsr")
##DICs = c(DICs, Bidart::DIC(run))
##run = readRDS("Heavy_metals/run_nsr_basis")
##DICs = c(DICs, Bidart::DIC(run))
##run = readRDS("Heavy_metals/run_nsr_basis_only")
##DICs = c(DICs, Bidart::DIC(run))
##run = readRDS("Heavy_metals/run_ns")
##DICs = c(DICs, Bidart::DIC(run))
##run = readRDS("Heavy_metals/run_ns_basis")
##DICs = c(DICs, Bidart::DIC(run))
##run = readRDS("Heavy_metals/run_stat")
##DICs = c(DICs, Bidart::DIC(run))
##remove(run);gc()


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
####KL_decomposition = readRDS("Heavy_metals/KL_decomposition.RDS") # loading karhunen loeve decomposition
####locs = readRDS("Heavy_metals/processed_data.RDS")$observed_locs # getting PP basis
####locs_ = unique(locs) 
####set.seed(1)
####locs_ = locs_[GpGp::order_maxmin(locs_),]
####locs_ = rbind(locs_, predicted_coords)
####NNarray = GpGp::find_ordered_nn(locs_, 5)
####PP_basis = Matrix::solve(
####  Matrix::sparseMatrix(i = row(NNarray)[!is.na(NNarray)], 
####                       j = NNarray[!is.na(NNarray)], 
####                       x = GpGp::vecchia_Linv(c(1, 3, 1, 0), "matern_isotropic", locs_, NNarray)[!is.na(NNarray)], 
####                       triangular = T
####  ), 
####  diag(1, nrow(locs_), 1000)
####)
####KL_basis = (PP_basis %*% KL_decomposition$v[,seq(20)] %*% diag(1/KL_decomposition$d[seq(20)]))[,seq(20)]
####colnames(KL_basis) = seq(20)
####saveRDS(KL_basis, "Heavy_metals/KL_basis.RDS")
####rm(PP_basis); rm(KL_decomposition); gc()
KL_basis = readRDS("Heavy_metals/KL_basis.RDS")
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


run = readRDS("Heavy_metals/run_nsr_basis_only")
X_pred_basis = cbind(X_pred, KL_basis[-seq(run$vecchia_approx$n_locs),])
basis = KL_basis[-seq(run$vecchia_approx$n_locs),]

#### predicting and estimating nonstat model
##run = readRDS("Heavy_metals/run_nsr_basis")
#prediction = Bidart::predict_latent_field(mcmc_nngp_list = run, predicted_locs = predicted_coords, X_range_pred = basis, X_scale_pred = basis, n_cores = 3)
#saveRDS(prediction, "Heavy_metals/prediction_field_nsr_basis.RDS")
#estimates = Bidart::estimate_parameters(run)
#saveRDS(estimates, "Heavy_metals/estimates_nsr_basis.RDS")
#rm(prediction) ; gc()

### predicting stat model
#run = readRDS("Heavy_metals/run_stat")
#prediction = Bidart::predict_latent_field(mcmc_nngp_list = run, predicted_locs = predicted_coords, X_range_pred = NULL, X_scale_pred = NULL, n_cores = 3)
#saveRDS(prediction, "Heavy_metals/prediction_field_stat.RDS")

# loading predictions
prediction_stat = readRDS("Heavy_metals/prediction_field_stat.RDS")
prediction_nsr = readRDS("Heavy_metals/prediction_field_nsr_basis.RDS")
prediction_nsr$predicted_samples = NULL
estimates_nsr = readRDS("Heavy_metals/estimates_nsr_basis.RDS")
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
for(i in seq(ncol(X_pred_basis)))
{
  name_ = paste("noise", name, sep = "_")
  gridmaps[[name_]] = NA
  gridmaps[[name_]][non_na_indices] = X_pred_basis[,i] * estimates_nsr$noise_beta[1,i,]
}
for(i in seq(ncol(basis)))
{
  name = colnames(X_pred_basis)[i]
  name_ = paste("range", name, sep = "_")
  gridmaps[[name_]] = NA
  gridmaps[[name_]][non_na_indices] = X_pred_basis[,i] * estimates_nsr$range_beta[1,i,]
  name_ = paste("scale", name, sep = "_")
  gridmaps[[name_]] = NA
  gridmaps[[name_]][non_na_indices] = X_pred_basis[,i] * estimates_nsr$scale_beta[1,i,]
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


# range
gridmaps[["range"]] = NA
#gridmaps[["range"]][non_na_indices] = cbind(1, X_pred_basis) %*% estimates_nsr$range_beta[1,,]
gridmaps[["range"]][non_na_indices] = cbind(1, basis) %*% estimates_nsr$range_beta[1,,]
pdf("Heavy_metals/range.pdf", height = 5, width = 10)
sp::plot(gridmaps["range"], main = "")
sp::plot(usa_aea, add = T)
dev.off()

##gridmaps[["range_covariates"]] = NA
##gridmaps[["range_covariates"]][non_na_indices] = X_pred %*% estimates_nsr$range_beta[1,seq(2, 12),]
##pdf("Heavy_metals/range_covariates", height = 5, width = 10)
##sp::plot(gridmaps["range_covariates"], main = "")
##sp::plot(usa_aea, add = T)
##dev.off()
##
##gridmaps[["range_basis"]] = NA
##gridmaps[["range_basis"]][non_na_indices] = cbind(1, X_pred_basis[,-seq(11)]) %*% estimates_nsr$range_beta[1,-seq(2, 12),]
##pdf("Heavy_metals/range_basis", height = 5, width = 10)
##sp::plot(gridmaps["range_basis"], main = "")
##sp::plot(usa_aea, add = T)
##dev.off()
##
##
##pdf("Heavy_metals/range_all.pdf", height = 4, width = 16)
##layout(matrix(1:4, 1, 4), widths = c(5,5,5,1))
##zlim = c(
##  min(c(gridmaps$range, gridmaps$range_covariates, gridmaps$range_basis), na.rm = T), 
##  max(c(gridmaps$range, gridmaps$range_covariates, gridmaps$range_basis), na.rm = T)
##)
##sp::plot(gridmaps["range"],  zlim = zlim, what = "image", main = "Total mean log range")
##sp::plot(usa_aea, add = T)
##sp::plot(gridmaps["range_covariates"], zlim = zlim, what = "image", main = "Mean log range explained by the environmental covariates")
##sp::plot(usa_aea, add = T)
##sp::plot(gridmaps["range_basis"], zlim = zlim, what = "image", main = "Mean log range explained by the spatial basis and the intercept")
##sp::plot(usa_aea, add = T)
##sp::plot(gridmaps["range_basis"], what = "scale", zlim = zlim )
##dev.off()



# scale
gridmaps[["scale"]] = NA
#gridmaps[["scale"]][non_na_indices] = cbind(1, X_pred_basis) %*% estimates_nsr$scale_beta[1,,]
gridmaps[["scale"]][non_na_indices] = cbind(1, basis) %*% estimates_nsr$scale_beta[1,,]
pdf("Heavy_metals/scale.pdf", height = 5, width = 10)
sp::plot(gridmaps["scale"], main = "")
sp::plot(usa_aea, add = T)
dev.off()

##gridmaps[["scale_covariates"]] = NA
##gridmaps[["scale_covariates"]][non_na_indices] = X_pred %*% estimates_nsr$scale_beta[1,seq(2, 12),]
##pdf("Heavy_metals/scale_covariates", height = 5, width = 10)
##sp::plot(gridmaps["scale_covariates"], main = "")
##sp::plot(usa_aea, add = T)
##dev.off()
##
##gridmaps[["scale_basis"]] = NA
##gridmaps[["scale_basis"]][non_na_indices] = cbind(1, X_pred_basis[,-seq(11)]) %*% estimates_nsr$scale_beta[1,-seq(2, 12),]
##pdf("Heavy_metals/scale_basis", height = 5, width = 10)
##sp::plot(gridmaps["scale_basis"], main = "")
##sp::plot(usa_aea, add = T)
##dev.off()
##
##
##pdf("Heavy_metals/scale_all.pdf", height = 4, width = 16)
##layout(matrix(1:4, 1, 4), widths = c(5,5,5,1))
##zlim = c(
##  min(c(gridmaps$scale, gridmaps$scale_covariates, gridmaps$scale_basis), na.rm = T), 
##  max(c(gridmaps$scale, gridmaps$scale_covariates, gridmaps$scale_basis), na.rm = T)
##)
##sp::plot(gridmaps["scale"],  zlim = zlim, what = "image", main = "Total mean log scale")
##sp::plot(usa_aea, add = T)
##sp::plot(gridmaps["scale_covariates"], zlim = zlim, what = "image", main = "Mean log scale explained by the environmental covariates")
##sp::plot(usa_aea, add = T)
##sp::plot(gridmaps["scale_basis"], zlim = zlim, what = "image", main = "Mean log scale explained by the spatial basis and the intercept")
##sp::plot(usa_aea, add = T)
##sp::plot(gridmaps["scale_basis"], what = "scale", zlim = zlim )
##dev.off()

# noise 

gridmaps[["noise"]] = NA
gridmaps[["noise"]][non_na_indices] = cbind(1, X_pred_basis) %*% estimates_nsr$noise_beta[1,,]
pdf("Heavy_metals/noise.pdf", height = 5, width = 10)
sp::plot(gridmaps["noise"], main = "")
sp::plot(usa_aea, add = T)
dev.off()

gridmaps[["noise_covariates"]] = NA
gridmaps[["noise_covariates"]][non_na_indices] = X_pred %*% estimates_nsr$noise_beta[1,seq(2, 12),]
pdf("Heavy_metals/noise_covariates", height = 5, width = 10)
sp::plot(gridmaps["noise_covariates"], main = "")
sp::plot(usa_aea, add = T)
dev.off()

gridmaps[["noise_basis"]] = NA
gridmaps[["noise_basis"]][non_na_indices] = cbind(1, X_pred_basis[,-seq(11)]) %*% estimates_nsr$noise_beta[1,-seq(2, 12),]
pdf("Heavy_metals/noise_basis", height = 5, width = 10)
sp::plot(gridmaps["noise_basis"], main = "")
sp::plot(usa_aea, add = T)
dev.off()


pdf("Heavy_metals/noise_all.pdf", height = 4, width = 16)
layout(matrix(1:4, 1, 4), widths = c(5,5,5,1))
zlim = c(
  min(c(gridmaps$noise, gridmaps$noise_covariates, gridmaps$noise_basis), na.rm = T), 
  max(c(gridmaps$noise, gridmaps$noise_covariates, gridmaps$noise_basis), na.rm = T)
)
sp::plot(gridmaps["noise"],  zlim = zlim, what = "image", main = "Total mean log noise")
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["noise_covariates"], zlim = zlim, what = "image", main = "Mean log noise explained by the environmental covariates")
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["noise_basis"], zlim = zlim, what = "image", main = "Mean log noise explained by the spatial basis and the intercept")
sp::plot(usa_aea, add = T)
sp::plot(gridmaps["noise_basis"], what = "scale", zlim = zlim )
dev.off()



#sp::plot(gridmaps["range_twi"], main = "", )
#sp::plot(gridmaps["range_globedem"], main = "", )
#sp::plot(gridmaps["range_gcarb"], main = "", )
#sp::plot(gridmaps["range_nlights03"], main = "", )
#sp::plot(gridmaps["range_sdroads"], main = "", )
#sp::plot(gridmaps["range_dairp"], main = "", )
#sp::plot(gridmaps["range_dmino"], main = "", )
#sp::plot(usa_aea, add = T)
#sp::plot(gridmaps["noise_sdroads"], main = "", )
#sp::plot(usa_aea, add = T)



#run = readRDS("Heavy_metals/run_nsr_basis")
run = readRDS("Heavy_metals/run_nsr_basis_only")

samples = do.call(rbind, 
  lapply(
    names(run$records$chain_1)[grep("beta", names(run$records$chain_1))], 
    function(name)
    {
      res = do.call(abind::abind, lapply(run$records, function(x)x[[name]][,,-seq(100)]))
      row.names(res) = paste(name, 
                             row.names(run$states$chain_1$params[[name]])
                             , sep = "_")
      res
    }
  ))
colnames(samples) = NULL

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

dev.off()

