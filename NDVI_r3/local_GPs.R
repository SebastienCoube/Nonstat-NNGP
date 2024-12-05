#############
# load data #
#############

.libPaths("../R_packages/")
#setwd("NDVI_r3")
observed_field = readRDS("train_field.RDS")
coords = readRDS("train_locs.RDS")
loo_field = readRDS("loo_field.RDS")
loo_coords = readRDS("loo_locs.RDS")
lump_field = readRDS("lump_center_field.RDS")
lump_coords = readRDS("lump_center_locs.RDS")

idx = seq(nrow(coords))
#idx = seq(nrow(coords))[order(runif(nrow(coords)))[seq(100000)]]
coords  = coords[idx,]
observed_field  = observed_field[idx]
Bidart::plot_pointillist_painting(coords, observed_field)

library(hetGP); library(lhs)
library(liGP); library(laGP)

model_noise = c(F, T, T)
model_range = c(F, T, T)

res = matrix(0, 3, 4)
res[,1] = c("lGP_stat", "lGP_het", "lGP_het_range")

unique_coords = coords[!duplicated(coords),]
coords_match = match(split(coords, row(coords)), split(unique_coords, row(unique_coords)))
unique_loo_coords = loo_coords[!duplicated(loo_coords),]
loo_coords_match = match(split(loo_coords, row(loo_coords)), split(unique_loo_coords, row(unique_loo_coords)))
unique_lump_coords = lump_coords[!duplicated(lump_coords),]
lump_coords_match = match(split(lump_coords, row(lump_coords)), split(unique_lump_coords, row(unique_lump_coords)))

t1 = Sys.time()
lhs_design <- randomLHS(10,2)
n <- 80
Xmt <- scale_ipTemplate(coords, n, space_fill_design=lhs_design, method='qnorm')$Xm.t

g_prior <- garg(list(mle=TRUE), observed_field)
theta_prior <- darg(NULL, coords)


liGP_smooth <- liGP(XX=unique_coords, X=coords, Y=observed_field, Xm=Xmt, N=n, 
                    theta = .01, 
                    g = g_prior, 
                    epsK=1e-5, num_thread = 30, 
                    reps = F, nu = 1.5)
liGP_loo <- liGP(XX=unique_loo_coords, X=coords, Y=observed_field, Xm=Xmt, N=n, 
                 theta = NULL, 
                 g = g_prior, 
                 epsK=1e-5, num_thread = 30, 
                 reps = F, nu = 1.5)
liGP_lump <- liGP(XX=unique_lump_coords, X=coords, Y=observed_field, Xm=Xmt, N=n, 
                  theta = NULL, 
                  g = g_prior, 
                  epsK=1e-5, num_thread = 30, 
                  reps = F, nu = 1.5)
## a few plots to make sure predictions are well matched
# plot(observed_field, liGP_smooth$mean[coords_match])
# abline(a = 0, b = 1)
# plot(loo_field,  liGP_loo$mean[loo_coords_match])
# abline(a = 0, b = 1)
# plot(lump_field,  liGP_lump$mean[lump_coords_match])
# abline(a = 0, b = 1)
smooth_sd = sqrt(liGP_smooth$nu * liGP_smooth$var)[coords_match]
loo_sd = sqrt(liGP_loo$nu * liGP_loo$var)[loo_coords_match]
lump_sd = sqrt(liGP_lump$nu * liGP_lump$var)[lump_coords_match]

scores = list(
  smooth_elpd = 
    mean(dnorm(observed_field, liGP_smooth$mean[coords_match], smooth_sd, log=T)),
  loo_elpd = 
    mean(dnorm(loo_field, liGP_loo$mean[loo_coords_match], loo_sd, log=T)), 
  lump_elpd = 
    mean(dnorm(lump_field, liGP_lump$mean[lump_coords_match], lump_sd, log=T)), 
  train_cover = mean((observed_field > liGP_smooth$mean[coords_match] + smooth_sd * qnorm(.025))&(observed_field < liGP_smooth$mean[coords_match] + smooth_sd * qnorm(.975))),
  loo_cover = mean((loo_field > liGP_loo$mean[loo_coords_match] + loo_sd * qnorm(.025))&(loo_field < liGP_loo$mean[loo_coords_match] + loo_sd * qnorm(.975))),
  lump_cover = mean((lump_field > liGP_lump$mean[lump_coords_match] + lump_sd * qnorm(.025))&(lump_field < liGP_lump$mean[lump_coords_match] + lump_sd * qnorm(.975))),
  loo_MSE = mean((loo_field - liGP_loo$mean[loo_coords_match])^2),
  lump_MSE = mean((lump_field - liGP_lump$mean[lump_coords_match])^2),
  time = as.numeric(Sys.time()-t1, units = "hours")
)
saveRDS(scores, "local_GPs_score.RDS")
