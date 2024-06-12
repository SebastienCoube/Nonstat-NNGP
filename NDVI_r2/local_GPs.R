#############
# load data #
#############

.libPaths("../R_packages/")
#setwd("NDVI_r2")
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


scores = list(
  smooth = 
    mean(dnorm(observed_field, liGP_smooth$mean[coords_match], sqrt(liGP_smooth$var + liGP_smooth$nu * liGP_smooth$var)[coords_match], log=T)),
  loo = 
    mean(dnorm(loo_field, liGP_loo$mean[loo_coords_match], sqrt(liGP_loo$var + liGP_loo$nu * liGP_loo$var)[loo_coords_match], log=T)), 
  lump = 
    mean(dnorm(lump_field, liGP_lump$mean[lump_coords_match], sqrt(liGP_lump$var + liGP_lump$nu * liGP_lump$var)[lump_coords_match], log=T)), 
  time = as.numeric(Sys.time()-t1, units = "hours")
)
saveRDS(scores, "local_GPs_score.RDS")
