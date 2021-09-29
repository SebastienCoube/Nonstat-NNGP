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
  M
}

#' @export
log_determinant_derivatives = function(sparse_chol_and_grad, NNarray)
{
  M = matrix(0, nrow(NNarray), length(sparse_chol_and_grad[[2]]))
  for( i in seq(length(sparse_chol_and_grad[[2]])))M[,i] = Bidart::log_determinant_derivative(derivative = sparse_chol_and_grad[[2]][[i]], compressed_sparse_chol = sparse_chol_and_grad[[1]], NNarray = NNarray)
  M
}

#' @export
expmat = function(coords)
{
  if(length(coords)==1)logm = matrix(coords)
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
variance_field = function(beta, field = NULL, X)
{
  res = X%*%matrix(beta, ncol = 1)
  if(!is.null(field)) 
  {
    res = res+field
  }
  res = exp(res)
  res = as.vector(res)
  return(res)
}

#' @export
compute_sparse_chol = function(covfun_name = covfun, range_beta, NNarray, locs, range_field = NULL, range_X = NULL, compute_derivative = T)
{
  
  if((covfun_name=="exponential_isotropic")|(covfun_name=="exponential_sphere")|(covfun_name=="exponential_spacetime"))  
  {
    res = list()
    res[[1]] = GpGp::vecchia_Linv(c(1, exp(range_beta), 0), covfun_name = covfun_name, locs = locs, NNarray = NNarray)
    res[[2]] = list()
    #for(i in seq(length(range_beta)))
    #{
    #  range_beta_ = range_beta
    #  range_beta_[i] = range_beta_[i] + 0.000001
    #  res[[2]][[i]] = (GpGp::vecchia_Linv(c(1, exp(range_beta_), 0), covfun_name = covfun_name, locs = locs, NNarray = NNarray) - res[[1]])/0.000001
    #}
    return(res)
  }
  if((covfun_name=="exponential_anisotropic2D"))  
  {
    res = list()
    #
    res[[1]] = GpGp::vecchia_Linv(c(1, solve(t(Bidart::expmat(range_beta)))[cbind(c(1, 2, 2), c(1, 1, 2))] , 0), covfun_name = covfun_name, locs = locs, NNarray = NNarray)
    res[[2]] = list()
    #for(i in seq(length(range_beta)))
    #{
    #  range_beta_ = range_beta
    #  range_beta_[i] = range_beta_[i] + 0.000001
    #  res[[2]][[i]] = (GpGp::vecchia_Linv(c(1, solve(t(Bidart::expmat(range_beta_)))[cbind(c(1, 2, 2), c(1, 1, 2))], 0), covfun_name = covfun_name, locs = locs, NNarray = NNarray) - res[[1]])/0.000001
    #}
    return(res)
  }
  if(!is.null(range_X)) log_range = range_X%*%range_beta
  if(!is.null(range_field)) log_range = log_range+range_field
  if(!is.null(log_range)) log_range = as.matrix(log_range)
  if(covfun_name=="nonstationary_exponential_isotropic")          res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_isotropic"  , sphere = F, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
  if(covfun_name=="nonstationary_exponential_isotropic_sphere")   res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_isotropic"  , sphere = T, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
  if(covfun_name=="nonstationary_exponential_anisotropic")        res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_anisotropic", sphere = F, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
  if(covfun_name=="nonstationary_exponential_anisotropic_sphere") res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_anisotropic", sphere = T, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
  res[[2]] = lapply(res[[2]], function(x)x*2)
  return(res)
}

#compute_sparse_chol = function(covfun_name = covfun, range_beta, NNarray, locs, range_field = NULL, range_X = NULL, compute_derivative = T)
#{
#  
#  if(covfun_name=="exponential_isotropic")   return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), 0), covfun_name = "exponential_isotropic", locs = locs, NNarray = NNarray)))
#  if(covfun_name=="exponential_sphere")      return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), 0), covfun_name = "exponential_sphere", locs = locs, NNarray = NNarray)))
#  if(covfun_name=="exponential_spacetime")   return(list(GpGp::vecchia_Linv(c(1, exp(range_beta), 0), covfun_name = "exponential_spacetime", locs = locs, NNarray = NNarray)))
#  if(covfun_name=="exponential_anisotropic") return(Bidart::nonstat_vecchia_Linv(log_range = matrix(1, nrow(locs), 1)%*%(range_beta*2), covfun_name = "nonstationary_exponential_anisotropic", sphere = F, locs = locs, NNarray = NNarray, compute_derivative = F))
#  if(!is.null(range_X)) log_range = range_X%*%range_beta
#  if(!is.null(range_field)) log_range = log_range+range_field
#  if(!is.null(log_range)) log_range = as.matrix(log_range)
#  if(covfun_name=="nonstationary_exponential_isotropic")          res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_isotropic"  , sphere = F, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
#  if(covfun_name=="nonstationary_exponential_isotropic_sphere")   res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_isotropic"  , sphere = T, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
#  if(covfun_name=="nonstationary_exponential_anisotropic")        res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_anisotropic", sphere = F, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
#  if(covfun_name=="nonstationary_exponential_anisotropic_sphere") res = Bidart::nonstat_vecchia_Linv(log_range = log_range*2, covfun_name = "nonstationary_exponential_anisotropic", sphere = T, locs = locs, NNarray = NNarray, compute_derivative = compute_derivative)
#  res[[2]] = lapply(res[[2]], function(x)x*2)
#  return(res)
#}

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
