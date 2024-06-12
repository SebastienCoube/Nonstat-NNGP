

.libPaths("../R_packages/")
library(GpGp);library(FNN);
library(Bidart);library(abind);
library(fields);library(parallel);
library(Matrix);
library(sp);
library(Bidart)
library(INLA)

inla.binary.install()
1
library(INLA)

#############
# load data #
#############
observed_field = readRDS("train_field.RDS")
coords = readRDS("train_locs.RDS")
coords_loo = readRDS("loo_locs.RDS")
coords_lump = readRDS("lump_center_locs.RDS")

#################
# get inla mesh #
#################

t1 = Sys.time()

prmesh = INLA::inla.mesh.2d(coords)
prmesh$n
plot(prmesh)

##############
# Covariates #
##############

# train covariates
X_split = split(coords, col(coords)); names(X_split) = c("coord_1", "coord_2")
X_split_loo = split(coords_loo, col(coords_loo)); names(X_split_loo) = c("coord_1", "coord_2")
X_split_lump = split(coords_lump, col(coords_lump)); names(X_split_lump) = c("coord_1", "coord_2")

########################################
# stationary spatially coherent effect #
########################################
spde_spatially_coherent = 
  INLA::inla.spde2.matern(
    mesh = prmesh, alpha=2 
    )
s.index = INLA::inla.spde.make.index(
  name = "spatial.field",
  spde_spatially_coherent$n.spde)


#######
# run #
#######
stk <- INLA::inla.stack(
  data = list(observed_field = observed_field),
  A = list(
    INLA::inla.spde.make.A(
      prmesh, # SPDE mesh
      loc=coords # coordinates of the observations of the interest variable
    ), 1
  ),
  effects=list(
    s.index,  #the spatial field index
    # the length of the spatial index is equal to the number of mesh nodes
    X_split #the covariates
  ))

stk_loo <- INLA::inla.stack(
  data = list(observed_field = NA),
  A = list(
    INLA::inla.spde.make.A(
      prmesh, # SPDE mesh
      loc=coords_loo # coordinates of the observations of the interest variable
    ), 1
  ),
  effects=list(
    s.index,  #the spatial field index
    # the length of the spatial index is equal to the number of mesh nodes
    X_split_loo #the covariates
  ))
stk_lump <- INLA::inla.stack(
  data = list(observed_field = NA),
  A = list(
    INLA::inla.spde.make.A(
      prmesh, # SPDE mesh
      loc=coords_lump # coordinates of the observations of the interest variable
    ), 1
  ),
  effects=list(
    s.index,  #the spatial field index
    # the length of the spatial index is equal to the number of mesh nodes
    X_split_lump #the covariates
  ))
 
stk = inla.stack(stk, stk_lump, stk_loo)

formula = observed_field ~
  coord_1 + coord_2 +  # fixed effects
  f(spatial.field, model = spde_spatially_coherent) # spatially coherent random effect

res <- INLA::inla(
  formula, data = INLA::inla.stack.data(stk), 
  family = "gaussian",
  control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(stk)),
  verbose = T, control.compute = list(cpo = TRUE, dic = TRUE)
)

saveRDS(res, "INLA_stat.RDS")

# Bidart::plot_pointillist_painting(
# prmesh$loc,
# res$summary.random$spatial.field$mean, 
# cex = .2
# )

# Plotting predicted values at observed
Bidart::plot_pointillist_painting(
  coords,
  res$summary.fitted.values$mean[seq(nrow(coords))],
  cex = .2
)
# predicting predicted values at lump centers
Bidart::plot_pointillist_painting(
  coords_lump,
  res$summary.fitted.values$mean[seq(nrow(coords)+1, nrow(coords)+nrow(coords_lump))],
  cex = .2
)
# predicting predicted values at loo locations
Bidart::plot_pointillist_painting(
  coords_loo,
  res$summary.fitted.values$mean[seq(nrow(coords)+nrow(coords_lump)+1, nrow(coords)+nrow(coords_lump)+nrow(coords_loo))],
  cex = .2
)



##################################
# log-score normal approximation #
##################################


score = list()
# observations 

score$train = mean(
  dnorm(
    observed_field,   
    res$summary.fitted.values$mean[seq(nrow(coords))], 
    sqrt( res$summary.fitted.values$sd[seq(nrow(coords))]^2 + 1/ res$summary.hyperpar["Precision for the Gaussian observations", "mode"]), 
    log = T
  )
)

# lump centers 
score$lump = mean(
  dnorm(
    readRDS("lump_center_field.RDS"),   
    res$summary.fitted.values$mean[seq(nrow(coords)+1, nrow(coords)+nrow(coords_lump))], 
    sqrt( res$summary.fitted.values$sd[seq(nrow(coords)+1, nrow(coords)+nrow(coords_lump))]^2 + 1/ res$summary.hyperpar["Precision for the Gaussian observations", "mode"]), 
    log = T
  )
)
# loo centers 

score$loo = mean(
  dnorm(
    readRDS("loo_field.RDS"),   
    res$summary.fitted.values$mean[seq(nrow(coords)+nrow(coords_lump)+1, nrow(coords)+nrow(coords_lump)+nrow(coords_loo))], 
    sqrt( res$summary.fitted.values$sd[seq(nrow(coords)+nrow(coords_lump)+1, nrow(coords)+nrow(coords_lump)+nrow(coords_loo))]^2 + 1/ res$summary.hyperpar["Precision for the Gaussian observations", "mode"]), 
    log = T
  )
)

score$time = as.numeric(Sys.time()-t1, units = "mins")
saveRDS(score, "inla_stat_score.RDS")