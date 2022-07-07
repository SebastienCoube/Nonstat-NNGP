install.packages(pkgs = "/home/user/s/scoube/Bidart_1.0.tar.gz", lib = "/home/user/s/scoube/R_packages/", repos = NULL, type = "source")
library(GpGp, lib.loc = "/home/user/s/scoube/R_packages/")
library(Bidart, lib.loc = "/home/user/s/scoube/R_packages/")
library(expm, lib.loc = "/home/user/s/scoube/R_packages/")
library(FNN, lib.loc = "/home/user/s/scoube/R_packages/")
library(abind, lib.loc = "/home/user/s/scoube/R_packages/")
library(parallel, lib.loc = "/home/user/s/scoube/R_packages/")
library(Matrix, lib.loc = "/home/user/s/scoube/R_packages/")


KL_nonstat_GP_NNGP = function(covmat, NNGP_precision, full_GP_det)
{
 
  H = 
    .5 * (nrow(covmat) * log(2*pi*exp(1)) + full_GP_det)
  KL = 
    .5 * (
      sum(NNGP_precision * covmat) 
      - nrow(covmat)
      - 2 *sum(log(vecchia_linv[,1]))
      - full_GP_det
    )
  c(H, KL)
}


seed = seq(30)
ordering = c("maxmin", "random", "coord", "middleout")
nparents = c(1, 5, 10, 20)
range_scale = c(.1, .3,   .5)
res = NULL
n = 10000
for(s in seed)
{
  for(sc in range_scale)
  {
    set.seed(s)
    locs = 5 * cbind(runif(n), runif(n))
    log_range = do.call(cbind, list(
      GpGp::fast_Gp_sim(covparms = c(sc, .5, 1, 0), locs =  locs, m = 50),
      GpGp::fast_Gp_sim(covparms = c(sc, .5, 1, 0), locs =  locs, m = 50), 
      GpGp::fast_Gp_sim(covparms = c(sc, .5, 1, 0), locs =  locs, m = 50)))
    log_range[,c(1, 2)] = log_range[,c(1, 2)] + log(.1)/sqrt(2)
    covmat = Bidart::nonstat_covmat(log_range = matrix(log_range, nrow = n), covfun_name = "nonstationary_exponential_anisotropic", locs = locs)
    full_GP_det = c(determinant(covmat)$mod)
    for(o in ordering)
    {
      for(m in nparents)
      {
        print(s)
        print(sc)
        print(o)
        print(m)
        print("")
        ordeuru = seq(n)
        if(o == "maxmin")oreduru = GpGp::order_maxmin(locs)
        if(o == "coord")oreduru = GpGp::order_coordinate(locs)
        if(o == "middleout")oreduru = GpGp::order_middleout(locs)
        locs_ = locs[oreduru,]
        log_range_ = log_range[oreduru,]
        NNarray = GpGp::find_ordered_nn(locs_, m = m)
        vecchia_linv = Bidart::nonstat_vecchia_Linv(log_range = log_range_, covfun_name = "nonstationary_exponential_anisotropic", locs = locs_, NNarray = NNarray, sphere = F, compute_derivative = F)[[1]]
        NNGP_precision = Matrix::crossprod(Matrix::sparseMatrix(
          i = row(NNarray)[!is.na(NNarray)], 
          j = NNarray[!is.na(NNarray)], 
          x = vecchia_linv[!is.na(NNarray)]
        ))[order(oreduru), order(oreduru)]
        res = rbind(res, c(s, sc, o, m,  KL_nonstat_GP_NNGP(covmat = covmat, NNGP_precision = NNGP_precision, full_GP_det = full_GP_det)))
      }
    }
  }
}
boxplot(as.numeric(res[,6]) ~ res[,2]+res[,3]+ res[,4])


saveRDS(res, "res_elliptic.RDS")