
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
t1 = Sys.time()


#idx = seq(nrow(coords))[order(runif(nrow(coords)))[seq(10000)]]
#coords  = coords[idx,]
#observed_field = observed_field[idx]

#################
# get inla mesh #
#################

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

##############
# Make noise #
##############

# making noise
PP = Bidart::get_PP(observed_locs = rbind(coords, coords_lump, coords_loo), matern_range = .2, 
                    n_PP = 50)
Bidart::compare_PP_NNGP(PP)
X_noise = cbind(1, as.matrix(Matrix::solve(PP$sparse_chol, diag(1, nrow(PP$sparse_chol), PP$n_PP))[-seq(PP$n_PP),][PP$idx,][seq(nrow(coords)),]))
Bidart::plot_pointillist_painting( rbind(coords), X_noise[,2], main = "PP 1")


t1 = Sys.time()

# A noise is an identity matrix with
# as many rows and columns as the number of observations
spde_noise <- INLA::inla.spde2.generic(
  M0 = Matrix::Diagonal(n=nrow(X_noise), x=1.0),
  M1 = Matrix::Diagonal(n=nrow(X_noise), x=0.0),
  M2 = Matrix::Diagonal(n=nrow(X_noise), x=0.0),
  B0 = cbind(0, X_noise), B1 = matrix(0, 1, 1+ ncol(X_noise)),
  B2 = matrix(0, 1, 1+ ncol(X_noise)),
  theta.mu = rep(0,ncol(X_noise)),
  theta.Q = diag(rep(.01,ncol(X_noise))),# weak prior
  transform = "identity")  
n.index = INLA::inla.spde.make.index(
  name = "heteroscedastic.noise",
  spde_noise$n.spde)



#########
# Stack #
#########
stk <- INLA::inla.stack(
  data = list(observed_field = observed_field),
  A = list(
    # noise
    Matrix::Diagonal(n=nrow(X_noise), x=1),
    # spatial
    INLA::inla.spde.make.A(
      prmesh, # SPDE mesh
      loc=coords # coordinates of the observations of the interest variable
    ),
    # fixed
    1
  ),
  effects=list(
    n.index, #the noise index
    s.index,  #the spatial field index
    # the length of the spatial index is equal to the number of mesh nodes
    X_split #the covariates
  )
  )

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

###################
# formula and run #
###################

formula = observed_field ~
  coord_1 + coord_2 +  # fixed effects
  f(heteroscedastic.noise, model = spde_noise)+
  f(spatial.field, model = spde_spatially_coherent) # spatially coherent random effect

res <- INLA::inla(
  formula, data = INLA::inla.stack.data(stk), 
  control.family = list(hyper = list(prec = list(initial = 20, fixed = TRUE))),
  control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(stk)),
  verbose = T, control.compute = list(cpo = TRUE, dic = TRUE), 
  control.inla=list(int.strategy="eb")
)


saveRDS(res, "runINLA_heterosk.RDS")

##########
# Scores #
##########
score = list()

noise_var = 
  exp(
    -2*
    # X noise
    cbind(1, as.matrix(Matrix::solve(PP$sparse_chol, diag(1, nrow(PP$sparse_chol), PP$n_PP))[-seq(PP$n_PP),][PP$idx,])) %*%
    # reg coeffs
    res$summary.hyperpar[grep("heteroscedastic.noise", dimnames(res$summary.hyperpar)[[1]]), "mode"]
  )
pdf("INLA_heterosk_var.pdf")
Bidart::plot_pointillist_painting(rbind(coords, coords_lump, coords_loo), noise_var)
dev.off()  
Bidart::plot_pointillist_painting(coords, X_noise )
# observations log-score
score$train =  mean(
    dnorm(
      observed_field,   
      res$summary.fitted.values$mean[seq(nrow(coords))], 
      (res$summary.fitted.values$sd[seq(nrow(coords))]^2 + noise_var[seq(nrow(coords))])^.5 , 
      log = T
    )
  )
  # lump centers log-score
score$lump =    mean(
    dnorm(
      readRDS("lump_center_field.RDS"),   
      res$summary.fitted.values$mean[seq(nrow(coords)+1, nrow(coords)+nrow(coords_lump))], 
      (res$summary.fitted.values$sd[seq(nrow(coords)+1, nrow(coords)+nrow(coords_lump))]^2 + noise_var[seq(nrow(coords)+1, nrow(coords)+nrow(coords_lump))])^.5 , 
      log = T
    )
  )
  # loo
score$loo = mean(
    dnorm(
      readRDS("loo_field.RDS"),   
      res$summary.fitted.values$mean[seq(nrow(coords)+nrow(coords_lump)+1, nrow(coords)+nrow(coords_lump)+nrow(coords_loo))], 
      (res$summary.fitted.values$sd[seq(nrow(coords)+nrow(coords_lump)+1, nrow(coords)+nrow(coords_lump)+nrow(coords_loo))]^2 + noise_var[seq(nrow(coords)+nrow(coords_lump)+1, nrow(coords)+nrow(coords_lump)+nrow(coords_loo))])^.5 , 
      log = T
    )
  )

score$time = as.numeric(Sys.time()-t1, units = "mins")
saveRDS(score, "inla_heterosk_score.RDS")


