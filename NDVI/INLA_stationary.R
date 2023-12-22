
library(GpGp);library(FNN);
library(Bidart);library(abind);
library(fields);library(parallel);
library(Matrix);
library(sp);
library(INLA)
library(Bidart)
INLA::inla.binary.install()

#############
# load data #
#############
observed_field = readRDS("train_field.RDS")
coords = readRDS("train_locs.RDS")
coords_test = readRDS("test_locs.RDS")

#################
# get inla mesh #
#################

lonlat_CRS = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
prmesh = INLA::inla.mesh.2d(coords)
prmesh$n
plot(prmesh)

##############
# Covariates #
##############

# train covariates
X = coords
X_split = split(X, col(X)); names(X_split) = c("coord_1", "coord_2")
X_split_test = split(coords_test, col(coords_test)); names(X_split_test) = c("coord_1", "coord_2")

########################################
# stationary spatially coherent effect #
########################################
A_field = INLA::inla.spde.make.A(prmesh, # SPDE mesh
                                 loc=coords # coordinates of the observations of the interest variable
)
spde_spatially_coherent = 
  INLA::inla.spde2.matern(
    mesh = prmesh, alpha=2 
    )
s.index = INLA::inla.spde.make.index(
  name = "spatial.field",
  spde_spatially_coherent$n.spde)


##############################
# run  with stationary noise #
##############################
 stk <- INLA::inla.stack(
   data = list(observed_field = observed_field),
   A = list(
            A_field, 
            1
   ),
   effects=list(
     s.index,  #the spatial field index
     # the length of the spatial index is equal to the number of mesh nodes
     X_split #the covariates
   ))
 
 formula = observed_field ~
   - 1 + # removing superfluous intercept
   coord_1 + coord_2 +  # fixed effects
   f(spatial.field, model = spde_spatially_coherent) # spatially coherent random effect
 
 
 
 clik <- list(hyper = list(prec = list(initial = 20,
                                       fixed = TRUE)))
 res <- INLA::inla(formula, data = INLA::inla.stack.data(stk), family = "gaussian",
                   control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(stk)), 
                   verbose = T, control.compute = list(cpo = TRUE, dic = TRUE))
 
 saveRDS(res, "runINLA_stat.RDS")

