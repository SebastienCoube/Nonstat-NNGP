#' @export
derivative_sandwiches = function(
  derivatives, 
  left_vector, 
  right_vector, 
  NNarray
)
{
  M = matrix(0, length(left_vector), length(derivatives))
  for( i in seq(length(derivatives)))M[,i] = Bidart::derivative_sandwich(derivatives[[i]], left_vector, right_vector, NNarray)
  # changing basis between det/aniso and canonical
  if(ncol(M)==3) M = M %*% t(matrix(
    c(1, 1, 0, 
      1,-1, 0,
      0, 0, 1), 3))
  M
}

#' @export
log_determinant_derivatives = function(sparse_chol_and_grad, NNarray)
{
  M = matrix(0, nrow(NNarray), length(sparse_chol_and_grad[[2]]))
  for( i in seq(length(sparse_chol_and_grad[[2]])))M[,i] = Bidart::log_determinant_derivative(derivative = sparse_chol_and_grad[[2]][[i]], compressed_sparse_chol = sparse_chol_and_grad[[1]], NNarray = NNarray)
  if(ncol(M)==3) M = M %*% matrix(
    c(1, 1, 0, 
      1,-1, 0,
      0, 0, 1), 3)
  M
}

#' @export
expmat = function(coords)
{
  if(length(coords)==1)logm = matrix(coords)
  if(length(coords)==3)
  {
    logm = matrix(coords [c(1, 3, 3, 2)] , 2)
  }
  if(length(coords)==6)
  {
    logm = matrix(0, 3, 3)
    diag(logm) = coords[seq(3)]
    logm[lower.tri(logm)] = coords[-seq(3)]
    logm[upper.tri(logm)] = logm[lower.tri(logm)]
  }
  expm::expm(logm)
}

#' @export
symmat = function(coords)
{
  if(length(coords)==3)
  {
    logm = matrix(coords [c(1, 3, 3, 2)], 2)
    logm[c(2, 3)] = logm[c(2, 3)]
  }
  if(length(coords)==6)
  {
    logm = matrix(0, 3, 3)
    diag(logm) = coords[seq(3)]
    logm[lower.tri(logm)] = coords[-seq(3)]
    logm[upper.tri(logm)] = logm[lower.tri(logm)]
  }
  logm
}

#' @export
variance_field = function(beta, KL = NULL, use_KL = F, X, locs_idx = NULL)
{
  as.vector(exp(X_KL_mult_right(X = X, KL = KL, use_KL = use_KL, locs_idx = locs_idx, Y = beta)))
}


# NB : computes sparse chol wrt determinant/anisotropy basis of range, 
# but gives derivatives wrt canonical
#' @export
compute_sparse_chol = function(covfun_name = covfun, range_beta, NNarray, locs, range_X = NULL, KL = NULL, use_KL = F, compute_derivative = T, nu = 1, locs_idx = NULL)
{
  if ((nu != 1) & (nu != 2)  & (nu != 1.5)) stop("nu must be equal to 1.5")
  
  # sqexp stat
  if(covfun_name=="exponential_isotropic") return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), 0), covfun_name = "exponential_isotropic", locs = locs, NNarray = NNarray)))
  if(covfun_name=="exponential_sphere")    return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), 0), covfun_name = "exponential_sphere", locs = locs, NNarray = NNarray)))
  if(covfun_name=="exponential_spacetime") return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), 0), covfun_name = "exponential_spacetime", locs = locs, NNarray = NNarray)))
  if(covfun_name=="exponential_anisotropic") return(nonstat_vecchia_Linv(log_range = matrix(1, nrow(locs), 1)%*%(range_beta*2), covfun_name = "nonstationary_exponential_anisotropic", sphere = F, locs = locs, NNarray = NNarray, compute_derivative = T, nu =1))
  
  # matern stat
  if(covfun_name=="matern_isotropic") return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), nu, 0), covfun_name = "matern_isotropic", locs = locs, NNarray = NNarray)))
  if(covfun_name=="matern15_isotropic") return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), 0), covfun_name = "matern15_isotropic", locs = locs, NNarray = NNarray)))
  if(covfun_name=="matern_sphere")    return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), nu, 0), covfun_name = "matern_sphere", locs = locs, NNarray = NNarray)))
  if(covfun_name=="matern_spacetime") return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), nu, 0), covfun_name = "matern_spacetime", locs = locs, NNarray = NNarray)))
  
  
  if(ncol(range_beta)==3)range_beta = range_beta %*% matrix(
    c(1, 1, 0, 
      1,-1, 0,
      0, 0, 1), 3)
  log_range = as.matrix(
    Bidart::X_KL_mult_right(
      X = range_X, 
      KL = KL, 
      Y = range_beta,  
      use_KL = use_KL, 
      locs_idx = locs_idx))
  #Bidart::plot_ellipses(locs, log_range)
  # exp locally isotropic
  if(covfun_name=="nonstationary_exponential_isotropic")          res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_isotropic"  , sphere = F, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative, nu = nu)
  if(covfun_name=="nonstationary_exponential_isotropic_sphere")   res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_isotropic"  , sphere = T, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative, nu = nu)
  
  # matern locally isotropic
  if(covfun_name=="nonstationary_matern_isotropic")          res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_matern_isotropic"  , sphere = F, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative, nu = nu)
  if(covfun_name=="nonstationary_matern_isotropic_sphere")   res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_matern_isotropic"  , sphere = T, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative, nu = nu)
  
  # exp locally anisotropic
  if(covfun_name=="nonstationary_exponential_anisotropic")        res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_anisotropic", sphere = F, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative, nu = nu)
  if(covfun_name=="nonstationary_exponential_anisotropic_sphere") res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_anisotropic", sphere = T, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative, nu = nu)
  
  # matern locally anisotropic
  if(covfun_name=="nonstationary_matern_anisotropic")        res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_matern_anisotropic", sphere = F, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative, nu = nu)
  if(covfun_name=="nonstationary_matern_anisotropic_sphere") res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_matern_anisotropic", sphere = T, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative, nu = nu)
  res[[2]] = lapply(res[[2]], function(x)x*2)
  return(res)
}

#' @export
get_KL_basis = function(locs, lonlat = F, covfun_name = "matern_isotropic", covparms, n_PP = 500, n_KL = 50, seed = 1, m = 5, get_basis = F)
{
  set.seed(seed)
  locs_ = locs[!duplicated(locs),]
  reordering = GpGp::order_maxmin(locs_, lonlat = lonlat)
  locs_ = locs_[reordering,]
  NNarray = GpGp::find_ordered_nn(locs_, m, lonlat = lonlat)
  sparse_chol = Matrix::sparseMatrix(
    i = row(NNarray)[!is.na(NNarray)], 
    j = NNarray[!is.na(NNarray)], 
    x = GpGp::vecchia_Linv(covparms = covparms, covfun_name = covfun_name, locs = locs_, NNarray = NNarray)[!is.na(NNarray)], 
    triangular = T
  )
  PP_basis = Matrix::solve(sparse_chol, Matrix::sparseMatrix(x = 1, i = seq(n_PP), j = seq(n_PP), dims = c(nrow(locs_), n_PP)))
  KL_decomposition = irlba::irlba(PP_basis, nu = n_KL, nv = n_KL)
  idx = match(split(locs, row(locs)), split(locs_, row(locs_)))
  basis = KL_decomposition$u[idx,]
  return(list("basis" = basis, "unique_reordered_locs" = locs_, "idx" = idx, "KL_decomposition" = KL_decomposition, "lonlat" = lonlat, "n_PP" = n_PP, "m" = m, "covparms" = covparms, "covfun_name" = covfun_name, "n_KL" = n_KL, "sparse_chol" = sparse_chol, "NNarray" = NNarray))
}

#' @export
predict_KL_basis = function(predicted_locs, KL_basis, seed = 1)
{
  set.seed(seed)
  locs_ = predicted_locs[!duplicated(predicted_locs),]
  #reordering = GpGp::order_maxmin(locs_, lonlat = KL_basis$lonlat)
  #locs_ = locs_[reordering,]
  
  locs_ = rbind(KL_basis$unique_reordered_locs, locs_)
  NNarray = GpGp::find_ordered_nn(locs_, KL_basis$m, lonlat = KL_basis$lonlat)
  
  
  sparse_chol = Matrix::sparseMatrix(
    i = row(NNarray)[!is.na(NNarray)], 
    j = NNarray[!is.na(NNarray)], 
    x = GpGp::vecchia_Linv(covparms = KL_basis$covparms, covfun_name = KL_basis$covfun_name, locs = locs_, NNarray = NNarray)[!is.na(NNarray)], 
    triangular = T
  )
  PP_basis = Matrix::solve(sparse_chol, Matrix::sparseMatrix(x = 1, i = seq(KL_basis$n_PP), j = seq(KL_basis$n_PP), dims = c(nrow(locs_), KL_basis$n_PP)))
  basis = PP_basis %*% KL_basis$KL_decomposition$v %*% Matrix::sparseMatrix(i = seq(KL_basis$n_KL), j = seq(KL_basis$n_KL), x = 1/KL_basis$KL_decomposition$d)
  basis = basis[match(split(predicted_locs, row(predicted_locs)), split(locs_, row(locs_))),]
  return(basis)
}

#' @export
beta_prior_ll = function(beta, n_KL, beta_mean, beta_precision, log_scale)
{
  if(n_KL>0) 
  {
    scale_mat = Bidart::expmat(-log_scale)
    return(
      (
          +.5 * n_KL * determinant(scale_mat, logarithm = T)$mod # determinant is changed by log scale
          -sum(.5 * (c(t(beta[seq(nrow(beta)-n_KL),, drop = F])-t(beta_mean)) %*% beta_precision) * c(t(beta[seq(nrow(beta)-n_KL),,drop = F])-t(beta_mean)))
          -sum(.5 * (beta[-seq(nrow(beta)-n_KL),,drop = F] %*% scale_mat) * beta[-seq(nrow(beta)-n_KL),,drop = F])
      )
    )
  }
  return(
    sum(
      c(
        -.5 * (c(t(beta[seq(nrow(beta)-n_KL),])-t(beta_mean)) %*% beta_precision) * c(t(beta[seq(nrow(beta)-n_KL),])-t(beta_mean))
      )
    )
  )
}


# #' @export
#beta_prior_ll_derivative = function(beta, n_KL, beta_mean, beta_precision, log_scale)
#{
#  if(n_KL>0) 
#  {
#    scale_mat = Bidart::expmat(-log_scale)
#    return(
#      matrix(
#        c(
#          -(c(t(beta[seq(nrow(beta)-n_KL),,drop=F])-t(beta_mean)) %*% beta_precision),
#          -c(t(beta[-seq(nrow(beta)-n_KL),,drop = F] %*% scale_mat))
#        ), 
#        nrow(beta)
#      )
#    )
#  }
#  return(
#    matrix(
#      c(
#        -(c(t(beta[seq(nrow(beta)-n_KL),])-t(beta_mean)) %*% beta_precision)
#      ), 
#      nrow(beta)
#    )
#  )
#}

#' @export
beta_prior_ll_derivative = function(beta, n_KL, beta_mean, beta_precision, log_scale)
{
  
  res = matrix(
    c(
      -(c(t(beta[seq(nrow(beta)-n_KL),, drop = F])-t(beta_mean)) %*% beta_precision)
    ), 
    nrow(beta)-n_KL
  )
  if(n_KL>0) 
  {
    scale_mat = Bidart::expmat(-log_scale)
    res = rbind(res, 
                -beta[-seq(nrow(beta)-n_KL),,drop = F] %*% scale_mat
    )
  }
  res
}


#beta = matrix(rnorm(10))
#n_KL = 5
#beta_mean = matrix(rnorm(5), 5, 1)
#beta_precision = diag(exp(rnorm(5)), 5, 5)
#log_scale = rnorm(1)
#
#beta_prior_ll(beta = beta, n_KL = 5, beta_mean = beta_mean, beta_precision = beta_precision, log_scale = log_scale)
#beta_ = beta; beta_[1]=beta_[1]+.0001
#beta_prior_ll_derivative(beta = beta, n_KL = 5, beta_mean = beta_mean, beta_precision = beta_precision, log_scale = log_scale)
#10000*
#  (
#  beta_prior_ll(beta = beta_, n_KL = 5, beta_mean = beta_mean, beta_precision = beta_precision, log_scale = log_scale)-
#  beta_prior_ll(beta = beta, n_KL = 5, beta_mean = beta_mean, beta_precision = beta_precision, log_scale = log_scale)
#  )


#KL$idx : match between the non redundant locations of KL and the redundant observed locations
#locs_idx : match between the redundant observed locations and those of X
X_KL_mult_right = function(X, KL = NULL, use_KL = F, locs_idx = NULL, Y)
{
  if(is.null(locs_idx))locs_idx = seq(nrow(X))
  Y = as.matrix(Y)
  res = X  %*% Y[seq(ncol(X)),]
  if(use_KL)
  {
    V = matrix(0, nrow(KL$sparse_chol), ncol(Y))
    V[seq(nrow(KL$KL_decomposition$v)),]= KL$KL_decomposition$v %*% Y[-seq(ncol(X)),, drop = F] 
    res = res + as.matrix(Matrix::solve(KL$sparse_chol, V, triangular = T))[KL$idx,,drop =F][locs_idx,,drop =F]
  }
  res
}


X_KL_crossprod = function(X, KL = NULL, use_KL = F,  Y, locs_idx = NULL)
{
  if(is.null(locs_idx))locs_idx = seq(nrow(X))
  Y = as.matrix(Y)
  res = matrix(0, ncol(X)+KL$n_KL*use_KL, ncol(Y))
  res[seq(ncol(X)),] = crossprod(x = X, y = Y)
  if(use_KL)
  {
    res[-seq(ncol(X)),] = 
      t(
        t(
          as.matrix(
            Matrix::solve(
              Matrix::t(KL$sparse_chol), 
              Matrix::sparseMatrix(x = 1, i = KL$idx, j = seq(length(KL$idx))) %*% ((Matrix::sparseMatrix(i = locs_idx, j = seq(nrow(Y)), dims = c(length(KL$idx), nrow(Y))) %*% Y))
              )
          )
        )[,1:KL$n_PP] %*% KL$KL_decomposition$v
      )
  }
  res
}


#' @export
derivative_chol_expmat = function(coords)
{
  dimres = 1
  if(length(coords)==6)dimres = 3
  res = array(data = 0, dim = c(dimres, dimres, length(coords)))
  chol_expmat = chol(Bidart::expmat(coords))
  for(i in seq(length(coords)))
  {
    coords_ = coords
    coords_[i] = coords_[i] + 0.00001
    res[,,i] = 100000 * (chol(Bidart::expmat(coords_)) - chol_expmat)
  }
  res
}

#' @export
derivative_field_wrt_scale = function(field, coords)
{
  d_chol_expmat = derivative_chol_expmat(coords)
  white_field = field %*% solve(chol(Bidart::expmat(coords)))
  res = array(0, dim = c(dim(field), length(coords)))
  for(i in seq(length(coords)))
  {
    res[,,i] = white_field %*% d_chol_expmat[,,i]
  }
  res
}
