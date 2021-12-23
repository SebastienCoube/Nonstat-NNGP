#' @export
get_array_summary = function(samples)
{
  if(any(is.na(samples)))message("there are NAs in the samples")
  if(is.null(samples))return(NULL)
  out = apply(samples, c(1, 2), function(x)c(mean(x, na.rm = T), quantile(x, c(0.025, 0.5, 0.975), na.rm = T), sd(x, na.rm = T)))
  dimnames(out)[[1]] = c("mean", "q0.025", "median", "q0.975", "sd")
  out
}

#' # Actually predicts more than the latent field
# 1) the fixed effects and the latent fields (if specified in the model) for the range and the scale
# 2) the latent field for the response variable
#' @export
predict_latent_field = function(mcmc_nngp_list, predicted_locs, X_range_pred = NULL, X_scale_pred = NULL, burn_in = .5, n_cores = 1, 
                                predict_range = F, predict_scale = F)
{
  # Sanity checks
  if(length(intersect(split(predicted_locs, row(predicted_locs)), split(mcmc_nngp_list$data$locs, row(mcmc_nngp_list$data$locs)) )) > 0) stop("the predicted locations must contain none of the original observed locations, use estimate instead")
  if(!is.null(X_range_pred))if(nrow(X_range_pred)!=nrow(predicted_locs))stop("X_range_pred must have the same number of rows as predicted_locs")
  if(!is.null(X_scale_pred))if(nrow(X_scale_pred)!=nrow(predicted_locs))stop("X_scale_pred must have the same number of rows as predicted_locs")
  # Dealing with duplicated data
  #if(any(duplicated(predicted_locs)))stop("there should not be duplicated locations")
  idx = match(split(predicted_locs, row(predicted_locs)), split(predicted_locs[!duplicated(predicted_locs),], row(predicted_locs[!duplicated(predicted_locs),])))
  X_range_pred     = X_range_pred[!duplicated(predicted_locs),]
  X_scale_pred     = X_scale_pred[!duplicated(predicted_locs),]
  predicted_locs = predicted_locs[!duplicated(predicted_locs),]
  # adding intercept to prediction regressors
  X_range_pred = cbind(rep(1, nrow(predicted_locs)), X_range_pred)
  X_scale_pred = cbind(rep(1, nrow(predicted_locs)), X_scale_pred)
  # binding regressors of observed and predicted locations 
  X_range = rbind(as.matrix(mcmc_nngp_list$data$covariates$range_X$X_locs), as.matrix(X_range_pred))
  if(ncol(X_range_pred)>1)X_range[-seq(mcmc_nngp_list$vecchia_approx$n_locs),-1] = X_range[-seq(mcmc_nngp_list$vecchia_approx$n_locs),-1] - matrix(rep(mcmc_nngp_list$data$covariates$range_X$X_mean[-1], each = nrow(X_range_pred)), nrow(X_range_pred))
  X_scale = rbind(as.matrix(mcmc_nngp_list$data$covariates$scale_X$X_locs), as.matrix(X_scale_pred))
  if(ncol(X_scale)>1)X_scale[-seq(mcmc_nngp_list$vecchia_approx$n_locs),-1] = X_scale[-seq(mcmc_nngp_list$vecchia_approx$n_locs),-1] - matrix(rep(mcmc_nngp_list$data$covariates$scale_X$X_mean[-1], each = nrow(X_scale_pred)), nrow(X_scale_pred))
  # constructing NNarray of observed and predicted locations
  locs = rbind(mcmc_nngp_list$data$locs, predicted_locs)
  NNarray = rbind(mcmc_nngp_list$vecchia_approx$NNarray, GpGp::find_ordered_nn(locs, m = ncol(mcmc_nngp_list$vecchia_approx$NNarray)-1)[-seq(mcmc_nngp_list$vecchia_approx$n_locs),])
  NNarray_non_NA = !is.na(NNarray)
  sparse_chol_row_idx = row(NNarray)[NNarray_non_NA]
  sparse_chol_column_idx = NNarray[NNarray_non_NA]
  # if needed hyprprior NNGP for covariance parameter fields
  if(!is.null(mcmc_nngp_list$hierarchical_model$hyperprior_covariance$scale_NNGP_prior)) 
  {
    scale_NNGP_prior = Matrix::sparseMatrix(
      i = sparse_chol_row_idx, 
      j = sparse_chol_column_idx, 
      x = GpGp::vecchia_Linv(covparms = c(1, mcmc_nngp_list$hierarchical_model$hyperprior_covariance$scale_NNGP_prior$range, 
                                          mcmc_nngp_list$hierarchical_model$hyperprior_covariance$scale_NNGP_prior$smoothness, 0), 
                             covfun_name = mcmc_nngp_list$hierarchical_model$hyperprior_covariance$log_NNGP_matern_covfun, 
                             NNarray = NNarray, 
                             locs = locs)[NNarray_non_NA]
    )
  }
  if(!is.null(mcmc_nngp_list$hierarchical_model$hyperprior_covariance$range_NNGP_prior)) 
  {
    range_NNGP_prior = Matrix::sparseMatrix(
      i = sparse_chol_row_idx, 
      j = sparse_chol_column_idx, 
      x = GpGp::vecchia_Linv(covparms = c(1, mcmc_nngp_list$hierarchical_model$hyperprior_covariance$range_NNGP_prior$range, 
                                          mcmc_nngp_list$hierarchical_model$hyperprior_covariance$range_NNGP_prior$smoothness, 0), 
                             covfun_name = mcmc_nngp_list$hierarchical_model$hyperprior_covariance$log_NNGP_matern_covfun, 
                             NNarray = NNarray, 
                             locs = locs)[NNarray_non_NA]
    )
  }
  
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  # parallelization
  n_samples = length(kept_iterations)
  cl = parallel::makeCluster(n_cores, outfile = NULL)
  predicted_samples = 
  parallel::parLapply(X = mcmc_nngp_list$records, cl = cl, fun = function(chain)
                                          {
                                            res = list()
                                            # creating arrays
                                            # fields
                                            res$scale_field = NULL 
                                            res$range_field = NULL
                                            if(!is.null(mcmc_nngp_list$hierarchical_model$hyperprior_covariance$scale_NNGP_prior))res$scale_field = array(0, dim = c(nrow(predicted_locs), 1, n_samples))
                                            if(!is.null(mcmc_nngp_list$hierarchical_model$hyperprior_covariance$range_NNGP_prior))res$range_field = array(0, dim = c(nrow(predicted_locs), dim(mcmc_nngp_list$states$chain_1$params$range_beta)[2], n_samples))
                                            res$field = array(0, dim = c(nrow(predicted_locs), 1, n_samples))
                                            # total effects
                                            if(predict_scale) res$log_scale = array(0, dim = c(nrow(predicted_locs), 1, n_samples))
                                            if(predict_range) res$log_range = array(0, dim = c(nrow(predicted_locs), dim(mcmc_nngp_list$states$chain_1$params$range_beta)[2], n_samples))
                                            # linear_effects effects
                                            if(predict_scale) res$log_scale_linear = array(0, dim = c(nrow(predicted_locs), 1, n_samples))
                                            if(predict_range) res$log_range_linear = array(0, dim = c(nrow(predicted_locs), dim(mcmc_nngp_list$states$chain_1$params$range_beta)[2], n_samples))
                                            # looping over saved observations
                                            for(i_predict in seq(n_samples))
                                            {
                                              gc()
                                              #  predicting scale latent field using log NNGP prior #######
                                              if(!is.null(mcmc_nngp_list$hierarchical_model$hyperprior_covariance$scale_NNGP_prior))
                                              {
                                                 res$scale_field[,,i_predict] = 
                                                   exp(.5 * chain$scale_log_scale[,,i_start + i_predict]) *
                                                   as.vector(Matrix::solve(
                                                     scale_NNGP_prior, 
                                                     c(
                                                       exp(-.5 * chain$scale_log_scale[,,i_start + i_predict]) * as.vector(mcmc_nngp_list$hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol %*% chain$scale_field[,,i_start + i_predict]), 
                                                       rnorm(nrow(predicted_locs))
                                                       )
                                                     ))[-seq(mcmc_nngp_list$vecchia_approx$n_locs)]
                                                 if(predict_scale)res$log_scale[,,i_predict] = res$scale_field[,,i_predict] + matrix(X_scale[-seq(mcmc_nngp_list$vecchia_approx$n_locs),], ncol = ncol(X_scale))  %*% chain$scale_beta[,,i_start + i_predict]
                                              }
                                              if(is.null(res$scale_field) & predict_scale)res$log_scale[,,i_predict] = matrix(X_scale[-seq(mcmc_nngp_list$vecchia_approx$n_locs),], ncol = ncol(X_scale)) %*% chain$scale_beta[,,i_start + i_predict]
                                              if(predict_scale)res$log_scale_linear[,,i_predict] = matrix(X_scale[-seq(mcmc_nngp_list$vecchia_approx$n_locs),], ncol = ncol(X_scale)) %*% chain$scale_beta[,,i_start + i_predict]
                                              #  predicting range latent field using log NNGP prior #######
                                              range_field_complete = NULL
                                              if(!is.null(mcmc_nngp_list$hierarchical_model$hyperprior_covariance$range_NNGP_prior))
                                              {
                                                 range_field_complete = 
                                                   (
                                                   as.matrix(Matrix::solve(
                                                     range_NNGP_prior, 
                                                     rbind(
                                                       as.matrix(mcmc_nngp_list$hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol
                                                                 %*% chain$range_field[,,i_start + i_predict]) 
                                                       %*% solve(chol(Bidart::expmat(chain$range_log_scale[,,i_start + i_predict]))), 
                                                       matrix(rnorm(nrow(predicted_locs) * ncol(chain$range_field)), ncol = ncol(chain$range_field))
                                                     )
                                                     ))
                                                   ) %*% chol(Bidart::expmat(chain$range_log_scale[,,i_start + i_predict]))
                                                 res$range_field[,,i_predict] = range_field_complete [-seq(mcmc_nngp_list$vecchia_approx$n_locs),]
                                                 if(predict_range)res$log_range[,,i_predict] = res$range_field[,,i_predict] + matrix(X_range[-seq(mcmc_nngp_list$vecchia_approx$n_locs),], ncol = ncol(X_range)) %*% chain$range_beta[,,i_start + i_predict]
                                              }
                                              if(is.null(res$range_field)& predict_range)res$log_range[,,i_predict] = matrix(X_range[-seq(mcmc_nngp_list$vecchia_approx$n_locs),], ncol = ncol(X_range)) %*% chain$range_beta[,,i_start + i_predict]
                                              if(predict_range) res$log_range_linear[,,i_predict] = matrix(X_range[-seq(mcmc_nngp_list$vecchia_approx$n_locs),], ncol = ncol(X_range)) %*% chain$range_beta[,,i_start + i_predict]
                                              # computing scale of the observations  ##########
                                              # computing NNGP factor
                                              sparse_chol = 
                                                Matrix::sparseMatrix(
                                                  i = sparse_chol_row_idx, 
                                                  j = sparse_chol_column_idx, 
                                                  x = 
                                                    Bidart::compute_sparse_chol(
                                                    covfun_name = mcmc_nngp_list$hierarchical_model$covfun, 
                                                    range_beta = chain$range_beta[,,i_start + i_predict], 
                                                    NNarray = NNarray, locs = locs, 
                                                    range_field = range_field_complete, 
                                                    range_X = X_range, 
                                                    compute_derivative = F, 
                                                    nu = mcmc_nngp_list$hierarchical_model$nu
                                                    )[[1]][NNarray_non_NA], 
                                                  triangular = T
                                                )
                                             scale_field = Bidart::variance_field(
                                               beta = chain$scale_beta[,,i_start + i_predict], 
                                               field = c(chain$scale_field[,,i_start + i_predict], res$scale_field[,,i_predict]), 
                                               X = X_scale
                                               )
                                             
                                              # predicting latent field #######
                                              res$field[,,i_predict] = 
                                                (
                                                 sqrt(scale_field) * 
                                                 as.vector(Matrix::solve(
                                                   sparse_chol, 
                                                   c(
                                                     as.vector(sparse_chol[1:mcmc_nngp_list$vecchia_approx$n_locs, 1:mcmc_nngp_list$vecchia_approx$n_locs] %*% (chain$field[,,i_start + i_predict]/sqrt(scale_field[1:mcmc_nngp_list$vecchia_approx$n_locs]))), 
                                                     rnorm(nrow(predicted_locs))
                                                   )
                                                 ))
                                                )[-seq(mcmc_nngp_list$vecchia_approx$n_locs)]
                                            }
                                            return(res)
                                          }
                      )
  summaries = list()
  predicted_samples_ = list()
  for(name in names(predicted_samples[[1]]))
  {
    predicted_samples_[[name]] = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(predicted_samples, function(y)y[[name]]))
  }
  summaries = parallel::parLapply(cl = cl, X = predicted_samples_, fun = Bidart::get_array_summary)
  for(name in names(predicted_samples[[1]]))
  {
    predicted_samples_[[name]] = predicted_samples_[[name]][idx,,]
    summaries[[name]] = summaries[[name]][,idx,]
  }
  return(list("predicted_locs" = predicted_locs,"predicted_samples" = predicted_samples_, "summaries" = summaries))
}




#' @export
predict_fixed_effects = function(mcmc_nngp_list, X_pred = NULL, burn_in = .5, individual_fixed_effects = NULL)
{
  if(!is.null(X_pred)){if(ncol(X_pred)!=ncol(mcmc_nngp_list$data$covariates$X$arg))stop("The number of provided covariates does not match")}
  # adding intercept to prediction regressors
  if(!is.null(X_pred))X_pred = cbind(rep(1, nrow(X_pred)), X_pred)
  if(is.null(X_pred))X_pred = matrix(1, 1, 1)
  colnames(X_pred)[1] = "(Intercept)"
  # which individual fixed effects are kept
  removed_fixed_effects = setdiff(c(colnames(X_pred)), individual_fixed_effects)
  if(!is.null(individual_fixed_effects))if(individual_fixed_effects == "all_fixed_effects")removed_fixed_effects= NULL
  # removing mean for prediction
  if(ncol(X_pred)>1)X_pred[,-1]  = X_pred[,-1] - matrix(rep(mcmc_nngp_list$data$covariates$X$X_mean[-1], each = nrow(X_pred)), nrow(X_pred))
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  n_samples = length(kept_iterations)
  predicted_samples = lapply(
                                          X = mcmc_nngp_list$records, 
                                          FUN = function(chain)
                                          {
                                            res = list()
                                            for(name in row.names(mcmc_nngp_list$states$chain_1$params$beta))
                                            {
                                              res[[name]] = array(0, dim = c(nrow(X_pred), 1, n_samples))
                                            }
                                            # looping over saved observations
                                            for(i_predict in seq(n_samples))
                                            {
                                             # range
                                              for(name in row.names(mcmc_nngp_list$states$chain_1$params$beta))
                                              {
                                                idx = match(name, row.names(mcmc_nngp_list$states$chain_1$params$beta))
                                                res[[name]][,,i_predict] = X_pred[,idx] * chain$beta[idx,,i_start + i_predict]
                                              }
                                            }
                                            res$total_linear_effects = Reduce("+", res)
                                            res[match(removed_fixed_effects, c(colnames(X_pred)))]=NULL
                                            return(res)
                                          })
  summaries = list()
  predicted_samples_ = list()
  for(name in names(predicted_samples[[1]]))
  {
    predicted_samples_[[name]] = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(predicted_samples, function(y)y[[name]]))
    summaries[[name]] = get_array_summary(predicted_samples_[[name]])
  }
  return(list("predicted_samples" = predicted_samples_, "summaries" = summaries))
}
#' @export
predict_noise = function(mcmc_nngp_list, X_noise_pred = NULL, burn_in = .5, individual_fixed_effects = NULL)
{
  message("only fixed effects for now")
  if(!is.null(X_noise_pred)){if(X_noise_pred =="No covariates were provided")X_noise_pred = NULL}
  if(is.null(X_noise_pred)&(mcmc_nngp_list$data$covariates$noise_X$arg!="No covariates were provided")) stop("No covariates were provided for prediction, while covariates were provided for fit")
  if(!is.null(X_noise_pred)&(mcmc_nngp_list$data$covariates$noise_X$arg=="No covariates were provided")) stop("Covariates were provided for prediction, while no covariates were provided for fit")
  if(!is.null(X_noise_pred)) if(ncol(X_noise_pred)!=ncol(mcmc_nngp_list$data$covariates$noise_X$arg))stop("The number of provided covariates does not match")
  # adding intercept to prediction regressors
  if(!is.null(X_noise_pred))X_noise_pred = cbind(rep(1, nrow(X_noise_pred)), X_noise_pred)
  if(is.null(X_noise_pred))X_noise_pred = matrix(1, 1, 1)
  colnames(X_noise_pred) = colnames(mcmc_nngp_list$data$covariates$noise_X$X)
  # which individual fixed effects are kept
  removed_fixed_effects = setdiff(colnames(X_noise_pred), individual_fixed_effects)
  if(!is.null(individual_fixed_effects))if(individual_fixed_effects == "all_fixed_effects")removed_fixed_effects= NULL
  # removing mean for prediction
  if(ncol(X_noise_pred)>1)X_noise_pred[,-1]  = X_noise_pred[,-1] - matrix(rep(mcmc_nngp_list$data$covariates$noise_X$X_mean[-1], each = nrow(X_noise_pred)), nrow(X_noise_pred))
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  # parallelization
  n_samples = length(kept_iterations)
  predicted_samples = lapply(X = mcmc_nngp_list$records, 
                             FUN = 
                             function(chain)
                             {
                               res = list()
                               # creating arrays
                               for(name in row.names(mcmc_nngp_list$states$chain_1$params$noise_beta))
                               {
                                 res[[name]] = array(0, dim = c(nrow(X_noise_pred), 1, n_samples))
                               }
                               # looping over saved observations
                               for(i_predict in seq(n_samples))
                               {
                                # range
                                 for(name in row.names(mcmc_nngp_list$states$chain_1$params$noise_beta))
                                 {
                                   idx = match(name, row.names(mcmc_nngp_list$states$chain_1$params$noise_beta))
                                   res[[name]][,,i_predict] = X_noise_pred[,idx] * chain$noise_beta[idx,,i_start + i_predict]
                                 }
                               }
                               res$total_linear_effects = Reduce("+", res)
                               res[match(removed_fixed_effects, c(colnames(X_noise_pred)))]=NULL
                               return(res)
                             })
  summaries = list()
  predicted_samples_ = list()
  for(name in names(predicted_samples[[1]]))
  {
    predicted_samples_[[name]] = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(predicted_samples, function(y)y[[name]]))
    summaries[[name]] = get_array_summary(predicted_samples_[[name]])
  }
  return(list("predicted_samples" = predicted_samples_, "summaries" = summaries))
}

#' @export
#estimate_parameters = function(mcmc_nngp_list, burn_in = .5, n_cores = 1)
#{
#  # burn in
#  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
#  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
#  summaries = list()
#  for(name in names(mcmc_nngp_list$records[[1]]))
#  {
#    if(length(grep("beta", name))==0)summaries[[name]] = get_array_summary(Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations))))))
#    # de-centering estimates of beta
#    if(name == "beta")
#    {
#      M = diag(1, ncol(mcmc_nngp_list$data$covariates$X$X), ncol(mcmc_nngp_list$data$covariates$X$X))
#      M[1, -1] = -mcmc_nngp_list$data$covariates$X$X_mean[-1]
#      summaries[[name]] = get_array_summary(array_multiply_1(Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations))))), M))
#      
#    }
#    if(name == "range_beta")
#    {
#      M = diag(1, ncol(mcmc_nngp_list$data$covariates$range_X$X_locs), ncol(mcmc_nngp_list$data$covariates$range_X$X_locs))
#      M[1, -1] = -mcmc_nngp_list$data$covariates$range_X$X_mean[-1]
#      summaries[[name]] = get_array_summary(array_multiply_1(Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations))))), M))
#    }
#    if(name == "scale_beta")
#    {
#      M = diag(1, ncol(mcmc_nngp_list$data$covariates$scale_X$X_locs), ncol(mcmc_nngp_list$data$covariates$scale_X$X_locs))
#      M[1, -1] = -mcmc_nngp_list$data$covariates$scale_X$X_mean[-1]
#      summaries[[name]] = get_array_summary(array_multiply_1(Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations))))), M))
#    }
#    if(name == "noise_beta")
#    {
#      M = diag(1, ncol(mcmc_nngp_list$data$covariates$noise_X$X), ncol(mcmc_nngp_list$data$covariates$noise_X$X))
#      M[1, -1] = -mcmc_nngp_list$data$covariates$noise_X$X_mean[-1]
#      summaries[[name]] = get_array_summary(array_multiply_1(Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations))))), M))
#    }
#    dimnames(summaries[[name]])[[2]] = dimnames(mcmc_nngp_list$states$chain_1$params[[name]])[[1]]
#  }
#  summaries
#}
estimate_parameters = function(mcmc_nngp_list, burn_in = .5, n_cores = 1, get_samples = F)
{
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  summaries = list()
  samples = list()
  for(name in names(mcmc_nngp_list$records[[1]]))
  {
    # matched latent field 
    x = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y$field[mcmc_nngp_list$vecchia_approx$locs_match,,kept_iterations], dim = c(mcmc_nngp_list$vecchia_approx$n_obs, 1 , length(kept_iterations)))))
    if(get_samples)samples$field_at_observed_locs = x
    summaries$field_at_observed_locs = get_array_summary(x)
    
    # general parameters minus beta
    if(length(grep("beta", name))==0)
    {
      x = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations)))))
      if(get_samples)samples[[name]] = x
      summaries[[name]] = get_array_summary(x)
    }
    # de-centering estimates of beta
    if(name == "beta")
    {
      M = diag(1, ncol(mcmc_nngp_list$data$covariates$X$X), ncol(mcmc_nngp_list$data$covariates$X$X))
      M[1, -1] = -mcmc_nngp_list$data$covariates$X$X_mean[-1]
      x = Bidart::array_multiply_1(Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations))))), M)
      if(get_samples)samples[[name]] = x
      summaries[[name]] = get_array_summary(x)
    }
    if(name == "range_beta")
    {
      M = diag(1, ncol(mcmc_nngp_list$data$covariates$range_X$X_locs), ncol(mcmc_nngp_list$data$covariates$range_X$X_locs))
      M[1, -1] = -mcmc_nngp_list$data$covariates$range_X$X_mean[-1]
      x = Bidart::array_multiply_1(Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations))))), M)
      if(get_samples)samples[[name]] = x
      summaries[[name]] = get_array_summary(x)
    }
    if(name == "scale_beta")
    {
      M = diag(1, ncol(mcmc_nngp_list$data$covariates$scale_X$X_locs), ncol(mcmc_nngp_list$data$covariates$scale_X$X_locs))
      M[1, -1] = -mcmc_nngp_list$data$covariates$scale_X$X_mean[-1]
      x = Bidart::array_multiply_1(Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations))))), M)
      if(get_samples)samples[[name]] = x
      summaries[[name]] = get_array_summary(x)
    }
    if(name == "noise_beta")
    {
      M = diag(1, ncol(mcmc_nngp_list$data$covariates$noise_X$X), ncol(mcmc_nngp_list$data$covariates$noise_X$X))
      M[1, -1] = -mcmc_nngp_list$data$covariates$noise_X$X_mean[-1]
      x = Bidart::array_multiply_1(Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations))))), M)
      if(get_samples)samples[[name]] = x
      summaries[[name]] = get_array_summary(x)
    }
    dimnames(summaries[[name]])[[2]] = dimnames(mcmc_nngp_list$states$chain_1$params[[name]])[[1]]
  }
  list("summaries" = summaries, "samples" = samples)
}

#test = estimate_parameters(mcmc_nngp_list)
#' @export
DIC = function(mcmc_nngp_list, burn_in = .5)
{
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  mean_field = rep(0, mcmc_nngp_list$vecchia_approx$n_locs)
  mean_noise = rep(0, mcmc_nngp_list$vecchia_approx$n_obs)
  mean_beta  = rep(0, ncol(mcmc_nngp_list$data$covariates$X$X))
  ll = c()
  for(i in kept_iterations)
  {
    for(j in seq(length(mcmc_nngp_list$states)))
    {
      noise = Bidart::variance_field(beta = mcmc_nngp_list$records[[j]]$noise_beta[,,i], field = mcmc_nngp_list$records[[j]]$noise_field[mcmc_nngp_list$vecchia_approx$locs_match,,i], X = mcmc_nngp_list$data$covariates$noise_X$X)
      field = mcmc_nngp_list$records[[j]]$field[,,i]
      beta = mcmc_nngp_list$records[[j]]$beta[,,i]
      ll = c(ll, 
             sum(dnorm(mcmc_nngp_list$data$observed_field - field[mcmc_nngp_list$vecchia_approx$locs_match]-mcmc_nngp_list$data$covariates$X$X%*%beta, 0, sqrt(noise), log = T))
             )
      mean_field = mean_field + field 
      mean_noise = mean_noise + noise 
      mean_beta  = mean_beta  + beta  
    }
  }
  mean_field = mean_field / (length(kept_iterations)*length(mcmc_nngp_list$states)) 
  mean_noise = mean_noise / (length(kept_iterations)*length(mcmc_nngp_list$states)) 
  mean_beta  = mean_beta  / (length(kept_iterations)*length(mcmc_nngp_list$states))
  ll_ = sum(dnorm(mcmc_nngp_list$data$observed_field - mean_field[mcmc_nngp_list$vecchia_approx$locs_match]-mcmc_nngp_list$data$covariates$X$X%*%mean_beta, 0, sqrt(mean_noise), log = T))
  -2 * (2 * mean(ll) - ll_)
}

#' @export
log_score_Gaussian = function(observed_field, latent_field_samples, log_noise_samples, fixed_effects_samples)
{
  if(length(dim(latent_field_samples))==2)nsamples =  ncol(latent_field_samples)
  if(length(dim(latent_field_samples))==0)nsamples = 1
  samples = dnorm(
    observed_field, 
    latent_field_samples + fixed_effects_samples, 
    exp(.5 * log_noise_samples), 
    log = T
  )
  list(
  "samples" =   samples, 
  "per_obs" =   apply(samples, 1, mean), 
  "total" = 
    sum(
      samples
    )/nsamples
  )
  
}

