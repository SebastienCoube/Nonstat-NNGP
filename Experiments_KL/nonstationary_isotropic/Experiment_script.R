

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
nparents = c(5, 10, 20)
range_scale = c(.1, .3,   .5)
res = NULL
n = 2500

for(s in seed)
{
  for(sc in range_scale)
  {
    set.seed(s)
    locs = 5 * cbind(runif(n), runif(n))
    #locs = as.matrix(expand.grid(seq(0, 5, .1), seq(0, 5, .1)))
    KL = Bidart::get_KL_basis(locs, lonlat = F, covfun_name = "matern15_isotropic", covparms = c(1, 1, .00001), n_PP = 500, n_KL = 30)
    log_range = log(.05) + sqrt(sc) * Bidart::X_KL_mult_right(Y = rnorm(31), KL = KL, X = cbind(rep(0, nrow(locs))), use_KL = T)
    covmat = Bidart::nonstat_covmat(log_range = matrix(log_range, nrow = n), covfun_name = "nonstationary_exponential_isotropic", locs = locs)
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
        log_range_ = matrix(log_range[oreduru], n)
        NNarray = GpGp::find_ordered_nn(locs_, m = m)
        vecchia_linv = Bidart::nonstat_vecchia_Linv(log_range = log_range_, covfun_name = "nonstationary_exponential_isotropic", locs = locs_, NNarray = NNarray, sphere = F, compute_derivative = F, nu = .5)[[1]]
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



#Bidart::plot_pointillist_painting(locs[oreduru,],as.vector(Matrix::solve(Matrix::sparseMatrix(
#  i = row(NNarray)[!is.na(NNarray)], 
#  j = NNarray[!is.na(NNarray)], 
#  x = vecchia_linv[!is.na(NNarray)]
#), rnorm(nrow(locs)))))
#Bidart::plot_pointillist_painting(locs[oreduru,],log_range_)

saveRDS(res, "Experiments_KL/nonstationary_isotropic/res_circ.RDS")
