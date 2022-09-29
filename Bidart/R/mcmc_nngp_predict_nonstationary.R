#' @export
get_array_summary = function(samples)
{
  if(any(is.na(samples)))message("there are NAs in the samples")
  if(is.null(samples))return(NULL)
  out = apply(samples, c(1, 2), function(x)c(mean(x, na.rm = T), quantile(x, c(0.025, 0.5, 0.975), na.rm = T), sd(x, na.rm = T)))
  dimnames(out)[[1]] = c("mean", "q0.025", "median", "q0.975", "sd")
  out
}

#' @export
predict_latent_field = function(mcmc_nngp_list, predicted_locs, X_range_pred = NULL, X_scale_pred = NULL, burn_in = .5, n_cores = 1)
{
  # Sanity checks
  if(length(intersect(split(predicted_locs, row(predicted_locs)), split(mcmc_nngp_list$data$locs, row(mcmc_nngp_list$data$locs)) )) > 0) stop("the predicted locations must contain none of the original observed locations, use estimate instead")
  if(!is.null(X_range_pred))if(nrow(X_range_pred)!=nrow(predicted_locs))stop("X_range_pred must have the same number of rows as predicted_locs")
  if(!is.null(X_scale_pred))if(nrow(X_scale_pred)!=nrow(predicted_locs))stop("X_scale_pred must have the same number of rows as predicted_locs")
  # Dealing with duplicated data
  idx = match(split(predicted_locs, row(predicted_locs)), split(predicted_locs[!duplicated(predicted_locs),], row(predicted_locs[!duplicated(predicted_locs),])))
  X_range_pred     = X_range_pred[!duplicated(predicted_locs),]
  X_scale_pred     = X_scale_pred[!duplicated(predicted_locs),]
  predicted_locs = predicted_locs[!duplicated(predicted_locs),]
  # adding intercept to prediction regressors
  X_range_pred = cbind(rep(1, nrow(predicted_locs)), X_range_pred)
  X_scale_pred = cbind(rep(1, nrow(predicted_locs)), X_scale_pred)
  # binding regressors of observed and predicted locations 
  X_range = rbind(as.matrix(mcmc_nngp_list$data$covariates$range_X$X_locs), as.matrix(X_range_pred))
  X_scale = rbind(as.matrix(mcmc_nngp_list$data$covariates$scale_X$X_locs), as.matrix(X_scale_pred))
  # constructing NNarray of observed and predicted locations
  locs = rbind(mcmc_nngp_list$data$locs, predicted_locs)
  NNarray = rbind(mcmc_nngp_list$vecchia_approx$NNarray, GpGp::find_ordered_nn(locs, m = ncol(mcmc_nngp_list$vecchia_approx$NNarray)-1)[-seq(mcmc_nngp_list$vecchia_approx$n_locs),])
  NNarray_non_NA = !is.na(NNarray)
  sparse_chol_row_idx = row(NNarray)[NNarray_non_NA]
  sparse_chol_column_idx = NNarray[NNarray_non_NA]
  # constructing NNarray of observed and predicted locations for parameter fields
  if(mcmc_nngp_list$hierarchical_model$scale_KL | mcmc_nngp_list$hierarchical_model$range_KL) 
  {
    NNarray_p = rbind(KL$NNarray, GpGp::find_ordered_nn(locs, m = ncol(KL$NNarray)-1)[-seq(nrow(KL$NNarray)),])
    NNarray_non_NA_p = !is.na(NNarray_p)
    sparse_chol_row_idx_p = row(NNarray_p)[NNarray_non_NA_p]
    sparse_chol_column_idx_p = NNarray_p[NNarray_non_NA_p]
    # if needed hyprprior NNGP for covariance parameter fields
    KL_ = mcmc_nngp_list$hierarchical_model$KL
    KL_$sparse_chol =  Matrix::sparseMatrix(
      i = sparse_chol_row_idx_p, 
      j = sparse_chol_column_idx_p, 
      x = GpGp::vecchia_Linv(covparms = KL$covparms, 
                             covfun_name = KL$covfun_name, 
                             NNarray = NNarray_p, 
                             locs = rbind(KL$unique_reordered_locs, locs[-seq(nrow(KL$unique_reordered_locs)),]))[NNarray_non_NA_p]
    )
    KL_$idx = c(KL$idx, mcmc_nngp_list$vecchia_approx$n_locs+idx)
  }
  
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  # parallelization
  n_samples = length(kept_iterations)
  cl = parallel::makeCluster(n_cores, outfile = NULL)
  predicted_samples = 
    parallel::parLapply(X = mcmc_nngp_list$records, 
                        cl = cl,
                        fun = function(chain)
                        {
                          res = list()
                          # creating arrays
                          # fields
                          res$field = array(0, dim = c(nrow(predicted_locs), 1, n_samples))
                          res$log_scale = array(0, dim = c(nrow(predicted_locs), 1, n_samples))
                          res$log_range = array(0, dim = c(nrow(predicted_locs), dim(mcmc_nngp_list$states$chain_1$params$range_beta)[2], n_samples))
                          # looping over saved observations
                          # i_predict= 1
                          for(i_predict in seq(n_samples))
                          {
                            gc()
                            # computing NNGP factor
                            sparse_chol = 
                              Matrix::sparseMatrix(
                                i = sparse_chol_row_idx, 
                                j = sparse_chol_column_idx, 
                                x = 
                                  Bidart::compute_sparse_chol(
                                    covfun_name = mcmc_nngp_list$hierarchical_model$covfun, 
                                    range_beta = chain$range_beta[,,i_start + i_predict, drop = F], 
                                    NNarray = NNarray, locs = locs, 
                                    KL = KL_, use_KL = mcmc_nngp_list$hierarchical_model$range_KL, 
                                    range_X = X_range, 
                                    compute_derivative = F, 
                                    nu = mcmc_nngp_list$hierarchical_model$nu, 
                                    locs_idx = c(mcmc_nngp_list$vecchia_approx$hctam_scol_1, seq(mcmc_nngp_list$vecchia_approx$n_obs+1, mcmc_nngp_list$vecchia_approx$n_obs+nrow(predicted_locs)))
                                  )[[1]][NNarray_non_NA], 
                                triangular = T
                              )
                            res$log_range[,,i_predict] = 
                              (
                                Bidart::X_KL_mult_right(
                                  X = X_range, 
                                  KL = KL_, use_KL = mcmc_nngp_list$hierarchical_model$range_KL, 
                                  locs_idx = c(mcmc_nngp_list$vecchia_approx$hctam_scol_1, seq(mcmc_nngp_list$vecchia_approx$n_obs+1, mcmc_nngp_list$vecchia_approx$n_obs+nrow(predicted_locs))), 
                                  Y = chain$range_beta[,,i_start + i_predict, drop = F]
                                )[-seq(mcmc_nngp_list$vecchia_approx$n_locs),]
                              )
                            scale_field = Bidart::variance_field(
                              beta = chain$scale_beta[,,i_start + i_predict], KL = KL_, use_KL = mcmc_nngp_list$hierarchical_model$scale_KL,  
                              locs_idx =  c(mcmc_nngp_list$vecchia_approx$hctam_scol_1, seq(mcmc_nngp_list$vecchia_approx$n_obs+1, mcmc_nngp_list$vecchia_approx$n_obs+nrow(predicted_locs))),
                              X = X_scale
                            )
                            res$log_scale[,,i_predict] = 
                              (
                                log(scale_field[-seq(mcmc_nngp_list$vecchia_approx$n_locs)])
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
predict_fixed_effects = function(mcmc_nngp_list, X_pred = NULL, burn_in = .5)
{
  if(!is.null(X_pred)){if(ncol(X_pred)!=ncol(mcmc_nngp_list$data$covariates$X$arg))stop("The number of provided covariates does not match")}
  # adding intercept to prediction regressors
  if(!is.null(X_pred))X_pred = cbind(rep(1, nrow(X_pred)), X_pred)
  if(is.null(X_pred))X_pred = matrix(1, 1, 1)
  colnames(X_pred)[1] = "(Intercept)"
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
predict_noise = function(mcmc_nngp_list, X_noise_pred = NULL, burn_in = .5, predicted_locs = NULL)
{
  if(mcmc_nngp_list$hierarchical_model$noise_KL & is.null(predicted_locs))stop("no predicted locs were provided, yet there is a GP prior for the noise")
  if(nrow(predicted_locs)!= nrow(X_noise_pred))stop("X_noise_pred and predicted_locs should have the same number of rows")
  if(!is.null(X_noise_pred)){if(identical(X_noise_pred, "No covariates were provided")){X_noise_pred = NULL}}
  if(is.null(X_noise_pred)&(!identical(mcmc_nngp_list$data$covariates$noise_X$arg,"No covariates were provided"))) stop("No covariates were provided for prediction, while covariates were provided for fit")
  if(!is.null(X_noise_pred)&(identical(mcmc_nngp_list$data$covariates$noise_X$arg, "No covariates were provided"))) stop("Covariates were provided for prediction, while no covariates were provided for fit")
  if(!is.null(X_noise_pred)) if(ncol(X_noise_pred)!=ncol(mcmc_nngp_list$data$covariates$noise_X$arg))stop("The number of provided covariates does not match")
  # adding intercept to prediction regressors
  if(!is.null(X_noise_pred))X_noise_pred = cbind(rep(1, nrow(X_noise_pred)), X_noise_pred)
  if(is.null(X_noise_pred))X_noise_pred = matrix(1, 1, 1)
  colnames(X_noise_pred) = colnames(mcmc_nngp_list$data$covariates$noise_X$X)
  # GP if needed
  if(mcmc_nngp_list$hierarchical_model$noise_KL) 
    {
      NNarray_p = rbind(KL$NNarray, GpGp::find_ordered_nn(locs, m = ncol(KL$NNarray)-1)[-seq(nrow(KL$NNarray)),])
      NNarray_non_NA_p = !is.na(NNarray_p)
      sparse_chol_row_idx_p = row(NNarray_p)[NNarray_non_NA_p]
      sparse_chol_column_idx_p = NNarray_p[NNarray_non_NA_p]
      # if needed hyprprior NNGP for covariance parameter fields
      KL_ = mcmc_nngp_list$hierarchical_model$KL
      KL_$sparse_chol =  Matrix::sparseMatrix(
        i = sparse_chol_row_idx_p, 
        j = sparse_chol_column_idx_p, 
        x = GpGp::vecchia_Linv(covparms = KL$covparms, 
                               covfun_name = KL$covfun_name, 
                               NNarray = NNarray_p, 
                               locs = rbind(KL$unique_reordered_locs, locs[-seq(nrow(KL$unique_reordered_locs)),]))[NNarray_non_NA_p]
      )
      KL_$idx = c(KL$idx, mcmc_nngp_list$vecchia_approx$n_locs+idx)
    }
  
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  # parallelization
  n_samples = length(kept_iterations)
  predicted_samples = lapply(
    X = mcmc_nngp_list$records, 
    FUN = 
      function(chain)
      {
        res = list(array(0, dim = c(nrow(X_noise_pred), 1, n_samples)))
        # looping over saved observations
        for(i_predict in seq(n_samples))
        {
          res[,,i] = Bidart::X_KL_mult_right(
            X = X_noise_pred, 
            KL = KL_, 
            use_KL = mcmc_nngp_list$hierarchical_model$noise_KL, 
            locs_idx = , 
            Y = chain$noise_beta[,,i_start + i_predict]
            
          )
        }
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
estimate_parameters = function(mcmc_nngp_list, burn_in = .5, n_cores = 1, get_samples = F)
{
  # burn in
  kept_iterations = seq(length(mcmc_nngp_list$iterations$thinning))[which(mcmc_nngp_list$iterations$thinning > mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in)]
  i_start = max(which(mcmc_nngp_list$iterations$thinning < mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]* burn_in))
  summaries = list()
  samples = list()
  # matched latent field 
  x = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y$field[mcmc_nngp_list$vecchia_approx$locs_match,,kept_iterations], dim = c(mcmc_nngp_list$vecchia_approx$n_obs, 1 , length(kept_iterations)))))
  if(get_samples)samples$field_at_observed_locs = x
  summaries$field_at_observed_locs = get_array_summary(x)
  for(name in names(mcmc_nngp_list$records[[1]]))
  {
    x = Reduce(f = function(x, y)abind::abind(x, y, along = 3), x = lapply(mcmc_nngp_list$records, function(y)array(y[[name]][,,kept_iterations], dim = c(dim(y[[name]])[c(1, 2)], length(kept_iterations)))))
    if(get_samples)samples[[name]] = x
    summaries[[name]] = get_array_summary(x)
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
      noise = Bidart::variance_field(
        beta = mcmc_nngp_list$records[[j]]$noise_beta[,,i], 
        KL = KL, mcmc_nngp_list$hierarchical_model$noise_KL, 
        X = mcmc_nngp_list$data$covariates$noise_X$X, 
        locs_idx = NULL
        )
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

