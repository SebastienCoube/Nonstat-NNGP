


#############
# load data #
#############

observed_field = readRDS("train_field.RDS")
coords = readRDS("train_locs.RDS")

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

# get PP at INLA mesh nodes #
PP = Bidart::get_PP(
  observed_locs = prmesh$loc, # spatial sites
  matern_range = .2,
  n_PP = 50, # number of knots
  m = 15 # number of NNGP parents
)

PP_basis_prmesh = as.matrix(Matrix::solve(PP$sparse_chol, diag(1, nrow(PP$sparse_chol), PP$n_PP))[-seq(PP$n_PP),][PP$idx,])
Bidart::plot_pointillist_painting(prmesh$loc, PP_basis_prmesh[,1])
PP_basis_obs = PP_basis_prmesh[match(split(coords, row(coords)), split(prmesh$loc[,-3], row(prmesh$loc[,-3]))),]



#######################################
# nonstationary heteroscedastic noise #
#######################################

X_noise_split = split(PP_basis_obs, col(PP_basis_obs)); names(X_noise_split) = paste("PP_", seq(ncol(PP_basis_obs)), sep="")

# PP_basis_obs corresponds to observations of the covariates
# ***at the coordinates of the observations***
A_noise = Matrix::Diagonal(n=nrow(PP_basis_obs), x=1)
# A noise is an identity matrix with
# as many rows and columns as the number of observations
spde_noise <- INLA::inla.spde2.generic(
  M0 = Matrix::Diagonal(n=nrow(PP_basis_obs), x=1.0),
  M1 = Matrix::Diagonal(n=nrow(PP_basis_obs), x=0.0),
  M2 = Matrix::Diagonal(n=nrow(PP_basis_obs), x=0.0),
  B0 = cbind(0, PP_basis_obs), B1 = matrix(0, 1, 1+ ncol(PP_basis_obs)),
  B2 = matrix(0, 1, 1+ ncol(PP_basis_obs)),
  theta.mu = rep(0,ncol(PP_basis_obs)),
  theta.Q = diag(rep(.00001,ncol(PP_basis_obs))), # very weak prior
  transform = "identity")  


########################################
# stationary spatially coherent effect #
########################################
A_field = INLA::inla.spde.make.A(prmesh, # SPDE mesh
                                 loc=coords # coordinates of the observations of the interest variable
)
spde_spatially_coherent = INLA::inla.spde2.matern(mesh = prmesh, alpha=2)


#################################################
# combining noise and spatially coherent effect #
#################################################


stk <- INLA::inla.stack(
  data = list(y = observed_field),
  A = list(A_noise,
           A_field,
           1
  ),
  effects=list(
    INLA::inla.spde.make.index(
      name = "heteroscedastic.noise",
      spde_noise$n.spde), #the noise index
    # the length of the noise index is equal to the number of observations
    INLA::inla.spde.make.index(
      name = "spatial.field",
      spde_spatially_coherent$n.spde),  #the spatial field index
    # the length of the spatial index is equal to the number of mesh nodes
    X_split #the covariates
  ))

formula = observed_field ~
  - 1 + # removing superfluous intercept
  coord_1 + coord_2 +  # fixed effects
  f(heteroscedastic.noise, model = spde_noise)+ # heteroskedastic Gaussian noise
  f(spatial.field, model = spde_spatially_coherent) # spatially coherent random effect


#################################
# run  with nonstationary noise #
#################################

clik <- list(hyper = list(prec = list(initial = 20,
                                      fixed = TRUE)))
res <- INLA::inla(formula, data = INLA::inla.stack.data(stk), control.family = clik,
                  control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(stk)), verbose = T)

saveRDS(res, "runINLA.RDS")