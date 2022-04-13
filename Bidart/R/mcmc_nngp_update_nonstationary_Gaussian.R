#' @export
update_centered_beta = function(beta, field, X, beta_prior, NNGP_prior, log_scale)
{
  centered_field = as.vector(field + X$X_locs %*% beta[X$which_locs]) # centering the field
  # including a priori
  posterior_precision = # (Q1 + Q2)
    exp(-log_scale) * X$crossprod_sparse_chol_X_locs + # Q1 precision wrt latent field info
    beta_prior$precision[X$which_locs, X$which_locs] # Q2 a priori conditional precision Q[which_locs|!which_locs] = Q[which_locs]
  posterior_mean = # (Q1 + Q2)^-1 (Q1 mu1 + Q2 mu2)
    solve(posterior_precision,  # (Q1 + Q2)^-1
          exp(-log_scale) * t (as.vector(NNGP_prior$sparse_chol %*% centered_field)  %*% X$sparse_chol_X_locs) # Q1 mu1 with mu1 = Q1^-1 XT Q_nngp w  (Q1 cancels out)
          + beta_prior$precision[X$which_locs, X$which_locs] %*% beta_prior$mean[X$which_locs] # Q2 mu2, with mu2 = (mu[which_locs]|mu[!which_locs]) = mu[which_locs] - Q[which_locs]Q[which_locs, !which_locs](beta[!which_locs]-mu[!which_locs])
          )
  if(ncol(X$X_locs)!=ncol(X$X)) 
  {
    posterior_mean = posterior_mean - 
    solve(posterior_precision, # (Q1 + Q2)
          solve(matrix(beta_prior$precision[X$which_locs, X$which_locs], length(X$which_locs)), # QAA-1
                matrix(beta_prior$precision[X$which_locs,- X$which_locs], length(X$which_locs)) %*% #QAB
                  (beta[-X$which_locs]- beta_prior$mean[-X$which_locs]))) #xB - muB
  }
  # Q2 mu2, with mu2 = (mu[which_locs]|mu[!which_locs]) = mu[which_locs] - Q[which_locs]Q[which_locs, !which_locs](beta[!which_locs]-mu[!which_locs])
  new_beta  = as.vector(posterior_mean + t(chol(solve(posterior_precision))) %*% rnorm(length(posterior_mean))) 
  new_field = centered_field - as.vector(X$X_locs %*% matrix(new_beta, ncol = 1))
  return(list(
    beta = new_beta, 
    field = new_field
  ))
}



#' @export
mcmc_nngp_update_Gaussian = function(data,
                                     hierarchical_model, vecchia_approx, # model architecture
                                     state, # model state
                                     n_iterations_update = 400, thinning = .1, iter_start = 0, seed = 1# practical settings
)
{
  #################
  # Sanity checks #
  #################
  # iterations and thinning
  if((floor(n_iterations_update)!=n_iterations_update) |  n_iterations_update<1)stop("n_iterations_update must be a positive round number")
  if((thinning<0)|(thinning>1))stop("thinning is a proportion and must be between 0 and 1")
  if((floor(iter_start)!=iter_start)|  iter_start<0)stop("iter_start must be a positive round number")
  # set seed 
  set.seed(seed)
  #########################################
  # Initializing chain storage structures #
  #########################################
  # this part re-creates a small portion of the $records objects of each chain. It fills it with chain state during the run, and then updates each chain with the new values
  params_records = lapply(state$params, function(x)array(0, c(dim(as.matrix(x)), round(thinning * n_iterations_update))))
  acceptance_records = lapply(state$transition_kernels, function(x)matrix(0, 50, length(x)))
  par(mfrow = c(2, 1))
  
  # do NOT remove or the code will bug in parallel. This is magic
  if(!is.null(state$params$range_log_scale))Bidart::expmat(state$params$range_log_scale)
  Matrix::tcrossprod(state$sparse_chol_and_stuff$sparse_chol, rep(1, vecchia_approx$n_locs))
  #################
  # Gibbs sampler #
  #################
  for(iter in seq(1, n_iterations_update))
  {
    gc()
    #cat(paste(iter, " "))
    #########
    # Range #
    #########
    #####################################
    # Stationary range beta (ancillary) #
    #####################################
    # Metropolis step
    if(length(grep("nonstat", hierarchical_model$covfun))==0)
    {
      new_range_beta_0  = state$params$range_beta + exp(.5*state$transition_kernels$range_beta_ancillary) * rnorm(length(state$params$range_beta))
      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                        range_beta = new_range_beta_0, NNarray = vecchia_approx$NNarray, 
                                                                        locs = data$locs, 
                                                                        range_field = NULL, 
                                                                        range_X = NULL, 
                                                                        compute_derivative = F, 
                                                                        nu = hierarchical_model$nu
      )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
      if(
        log(runif(1))<
        -.5 * sum((new_field         [vecchia_approx$locs_match]- state$sparse_chol_and_stuff$lm_residuals)^2/state$sparse_chol_and_stuff$noise)
        +.5 * sum((state$params$field[vecchia_approx$locs_match]- state$sparse_chol_and_stuff$lm_residuals)^2/state$sparse_chol_and_stuff$noise)
        # prior
        - .5 * (new_range_beta_0        - hierarchical_model$beta_priors$range_beta$mean)^2 * hierarchical_model$beta_priors$range_beta$precision
        + .5 * (state$params$range_beta - hierarchical_model$beta_priors$range_beta$mean)^2 * hierarchical_model$beta_priors$range_beta$precision
      )
      {
        state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
        state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
        state$params$range_beta[] = new_range_beta_0
        state$params$field = new_field
        state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
        acceptance_records$range_beta_ancillary[iter - 50*(iter %/% 50)] = acceptance_records$range_beta_ancillary[iter - 50*(iter %/% 50)]+1
      }
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_beta_ancillary)>.51)state$transition_kernels$range_beta_ancillary = state$transition_kernels$range_beta_ancillary + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_beta_ancillary)<.11)state$transition_kernels$range_beta_ancillary = state$transition_kernels$range_beta_ancillary - rnorm(1, .4, .05)
          acceptance_records$range_beta_ancillary =  0*acceptance_records$range_beta_ancillary
        }
      }
      # updating MALA kernel
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_beta_ancillary)<.11)state$transition_kernels$range_beta_ancillary = state$transition_kernels$range_beta_ancillary - rnorm(1, .2, .05)
          acceptance_records$range_beta_ancillary =  0*acceptance_records$range_beta_ancillary
        }
      }
    }
    ######################################
    # Stationary range beta (sufficient) #
    ######################################
    # Metropolis step
    if(length(grep("nonstat", hierarchical_model$covfun))==0)
    {
      new_range_beta_0  = state$params$range_beta + exp(.5*state$transition_kernels$range_beta_sufficient) * rnorm(length(state$params$range_beta))
      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                        range_beta = new_range_beta_0, NNarray = vecchia_approx$NNarray, 
                                                                        locs = data$locs, 
                                                                        range_field = NULL, 
                                                                        range_X = NULL, 
                                                                        compute_derivative = F, 
                                                                        nu = hierarchical_model$nu
      )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      
      if(log(runif(1))< 
         sum(log(new_compressed_sparse_chol_and_grad[[1]][,1])) - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
         -.5 * sum((new_sparse_chol                         %*% (state$params$field /sqrt(state$sparse_chol_and_stuff$scale)))^2)
         +.5 * sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field /sqrt(state$sparse_chol_and_stuff$scale)))^2)
         # prior
         - .5 * (new_range_beta_0        - hierarchical_model$beta_priors$range_beta$mean)^2 * hierarchical_model$beta_priors$range_beta$precision
         + .5 * (state$params$range_beta - hierarchical_model$beta_priors$range_beta$mean)^2 * hierarchical_model$beta_priors$range_beta$precision
      )
      {
        state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
        state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
        state$params$range_beta[] = new_range_beta_0
        state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
        acceptance_records$range_beta_sufficient[iter - 50*(iter %/% 50) ] = acceptance_records$range_beta_sufficient[iter - 50*(iter %/% 50) ]+1
      }
      # updating kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_beta_sufficient)>.41)state$transition_kernels$range_beta_sufficient = state$transition_kernels$range_beta_sufficient + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_beta_sufficient)<.11)state$transition_kernels$range_beta_sufficient = state$transition_kernels$range_beta_sufficient - rnorm(1, .4, .05)
          acceptance_records$range_beta_sufficient =  0*acceptance_records$range_beta_sufficient
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_beta_sufficient)<.11)state$transition_kernels$range_beta_sufficient = state$transition_kernels$range_beta_sufficient - rnorm(1, .2, .05)
          acceptance_records$range_beta_sufficient =  0*acceptance_records$range_beta_sufficient
        }
      }
    }
    ###############
    ########################################
    # Nonstationary range beta (ancillary) #
    ########################################
    if(length(grep("nonstat", hierarchical_model$covfun))==1)
    {
      q = t(solve(data$covariates$range_X$chol_solve_crossprod_X_locs)) %*% state$params$range_beta # whitening wrt covariates of the range
      current_U =
        as.vector(
          + .5*t(c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened) %*% hierarchical_model$beta_priors$range_beta$precision_whitened %*% (c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened) # normal prior 
          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        )
      # MALA whitened
      state$momenta$range_beta_ancillary = sqrt(.9) * state$momenta$range_beta_ancillary + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
      p = state$momenta$range_beta_ancillary
      # Make a half step for momentum at the beginning
      p = p - exp(state$transition_kernels$range_beta_ancillary) *
        (+ as.matrix(
          0.5* matrix(hierarchical_model$beta_priors$range_beta$precision_whitened %*% (c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened), ncol = ncol(q))
          + solve(solve(data$covariates$range_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                 # normal prior derivative                
                 t(data$covariates$range_X$X_locs) %*%   # Jacobian of range field wrt range_beta
                   (
                     # natural gradient of obs likelihood wrt range field
                     Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                   left_vector = as.vector(
                                                     Matrix::solve(
                                                       Matrix::t(state$sparse_chol_and_stuff$sparse_chol), 
                                                       - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                     ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                                       * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                     )), 
                                                   right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                   NNarray = vecchia_approx$NNarray  
                     )
                   ))) 
        )/ 2
      # Make a full step for the position
      q = q + exp(state$transition_kernels$range_beta_ancillary) * p
      new_range_beta = t(data$covariates$range_X$chol_solve_crossprod_X_locs) %*% q 
      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                        range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
                                                                        locs = data$locs, 
                                                                        range_field = state$params$range_field, 
                                                                        range_X = data$covariates$range_X$X_locs, 
                                                                        nu = hierarchical_model$nu
      )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
      # Make a half step for momentum at the end.
      p = p - exp(state$transition_kernels$range_beta_ancillary) *
        (
         + as.matrix(
           0.5* matrix(hierarchical_model$beta_priors$range_beta$precision_whitened %*% (c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened), ncol = ncol(q))
           + solve(solve(data$covariates$range_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                 # normal prior derivative                
                 t(data$covariates$range_X$X_locs) %*%   # Jacobian of range field wrt range_beta
                   (
                     # natural gradient of obs likelihood wrt range field
                     Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                   left_vector = as.vector(
                                                     Matrix::solve(
                                                       Matrix::t(new_sparse_chol), 
                                                       - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                     ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                                       * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                     )), 
                                                   right_vector = new_field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                   NNarray = vecchia_approx$NNarray  
                     )
                   ))) 
        )/ 2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$range_beta_ancillary ^2) / 2
      proposed_U =
        as.vector (
          + .5*t(c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened) %*% hierarchical_model$beta_priors$range_beta$precision_whitened %*% (c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened) # normal prior 
          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        )
      proposed_K = sum(p^2) / 2
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
        {
          state$momenta$range_beta_ancillary = p
          state$params$field = new_field
          state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
          state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
          state$params$range_beta[] = new_range_beta
          state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          acceptance_records$range_beta_ancillary[iter - 50*(iter %/% 50) ] = acceptance_records$range_beta_ancillary[iter - 50*(iter %/% 50) ] + 1
        }
      }
      
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_beta_ancillary)>.91)state$transition_kernels$range_beta_ancillary = state$transition_kernels$range_beta_ancillary + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_beta_ancillary)<.61)state$transition_kernels$range_beta_ancillary = state$transition_kernels$range_beta_ancillary - rnorm(1, .4, .05)
          acceptance_records$range_beta_ancillary =  0*acceptance_records$range_beta_ancillary
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))      
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_beta_ancillary)<.61)state$transition_kernels$range_beta_ancillary = state$transition_kernels$range_beta_ancillary - rnorm(1, .2, .05)
          acceptance_records$range_beta_ancillary =  0*acceptance_records$range_beta_ancillary
        }
      }
    }
    #######################################
    # Nonstationary range beta (centered) #
    #######################################
###    if(!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior))
###    {
###      # rotating field and beta in order to get 3 decorrelated components (each component has spatial correlation tho)
###      rotated_field = state$params$range_field  %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
###      rotated_beta  = matrix(state$params$range_beta, nrow = ncol(data$covariates$range_X$X_locs)) %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
###      # center
###      centered_field = as.matrix(rotated_field + data$covariates$range_X$X_locs %*% rotated_beta)
###      sparse_chol_X = data$covariates$range_X$sparse_chol_X_locs
###      beta_covmat = data$covariates$range_X$solve_crossprod_sparse_chol_X_locs
###      # beta mean =  w^t RtRX (XtRtRX)^{-1}(((Xt))) 
###      beta_mean =  t(
###        t(as.matrix(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol %*% centered_field))  %*% sparse_chol_X #vbw^t RtRX 
###        %*% beta_covmat # (XtRtRX)^{-1} 
###      )
###      # sampling new beta
###      rotated_beta = beta_mean +  t(chol(beta_covmat)) %*% matrix(rnorm(length(beta_mean)), nrow = nrow(beta_mean))
###      # de-centering and de-rotating
###      rotated_field = centered_field - as.matrix(data$covariates$range_X$X_locs %*% rotated_beta)
###      state$params$range_field = rotated_field %*% chol(Bidart::expmat(state$params$range_log_scale))
###      state$params$range_beta[] = rotated_beta %*% chol(Bidart::expmat(state$params$range_log_scale))
###    }
    if(!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior))
    {
      
      # rotating field and beta in order to get 3 decorrelated components (each component has spatial correlation tho)
      rotated_field = state$params$range_field  %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
      rotated_beta  = matrix(state$params$range_beta, nrow = ncol(data$covariates$range_X$X_locs)) %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
      # data knowing parameter
      centered_field = as.matrix(rotated_field + data$covariates$range_X$X_locs %*% rotated_beta)
      sparse_chol_X = data$covariates$range_X$sparse_chol_X_locs
      beta_covmat = data$covariates$range_X$solve_crossprod_sparse_chol_X_locs
      # prior (whitened wrt log scale)
      rotation_matrix = t(solve(chol(Bidart::expmat(state$params$range_log_scale)))) %x% diag(1, nrow(rotated_beta), nrow(rotated_beta))
      solve_rotation_matrix = solve(rotation_matrix)
      beta_precision_prior = t(solve_rotation_matrix) %*% hierarchical_model$beta_priors$range_beta$precision %*% solve_rotation_matrix
      # combining data knowing parameter with prior 
        # Q1 + Q2
      beta_precision_posterior = diag(1, 3, 3) %x% solve(beta_covmat) + beta_precision_prior
        # (Q1 + Q2)-1 (Q1 mu1 + Q2 mu2)
      beta_mean_posterior = solve(beta_precision_posterior, 
                                  c(
                                    t(as.matrix(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol %*% centered_field))  %*% sparse_chol_X # Q1 mu1 
                                    )
                                  + 
                                    beta_precision_prior %*% rotation_matrix %*% hierarchical_model$beta_priors$range_beta$mean
                                  )
      
      # sampling new beta
      rotated_beta = matrix(beta_mean_posterior +  t(chol(solve(beta_precision_posterior))) %*% rnorm(length(beta_mean_posterior)), nrow = nrow(state$params$range_beta))
      # de-centering and de-rotating
      rotated_field = centered_field - as.matrix(data$covariates$range_X$X_locs %*% rotated_beta)
      state$params$range_field = rotated_field %*% chol(Bidart::expmat(state$params$range_log_scale))
      state$params$range_beta[] = rotated_beta %*% chol(Bidart::expmat(state$params$range_log_scale))
    }
    #########################################
    # Nonstationary range beta (sufficient) #
    #########################################
    if(length(grep("nonstat", hierarchical_model$covfun))==1)
    {
      q = t(solve(data$covariates$range_X$chol_solve_crossprod_X_locs)) %*% state$params$range_beta # whitening wrt covariates of the range
      current_U = as.vector(
          +.5*t(c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened) %*% hierarchical_model$beta_priors$range_beta$precision_whitened %*% (c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened) # normal prior 
          + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
        )
      
      # MALA whitened
      state$momenta$range_beta_sufficient = sqrt(.9) * state$momenta$range_beta_sufficient + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
      p = state$momenta$range_beta_sufficient
      # Make a half step for momentum at the beginning
      p = p - exp(state$transition_kernels$range_beta_sufficient) *
        (+ as.matrix(
           0.5* matrix(hierarchical_model$beta_priors$range_beta$precision_whitened %*% (c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened), ncol = ncol(q))  # normal prior derivative                
           + solve(solve(data$covariates$range_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                 t(data$covariates$range_X$X_locs) %*%   # Jacobian of range field wrt range_beta
                   (# natural gradient of obs likelihood wrt range field
                     Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                   left_vector = as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                   right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                   NNarray = vecchia_approx$NNarray  
                     )
                     - Bidart::log_determinant_derivatives(sparse_chol_and_grad = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
                   ))) 
        )/ 2
      # Make a full step for the position
      q = q + exp(state$transition_kernels$range_beta_sufficient) * p
      new_range_beta = t(data$covariates$range_X$chol_solve_crossprod_X_locs) %*% q 
      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                        range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
                                                                        locs = data$locs, 
                                                                        range_field = state$params$range_field, 
                                                                        range_X = data$covariates$range_X$X_locs, 
                                                                        nu = hierarchical_model$nu
      )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      # Make a half step for momentum at the end.
      p = p - exp(state$transition_kernels$range_beta_sufficient) *
        (        
         + as.matrix(
           0.5* matrix(hierarchical_model$beta_priors$range_beta$precision_whitened %*% (c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened), ncol = ncol(q))  # normal prior derivative
           + solve (solve(data$covariates$range_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  t(data$covariates$range_X$X_locs) %*%   # Jacobian of range field wrt range_beta
                   (# natural gradient of obs likelihood wrt range field
                     Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                   left_vector = as.vector(new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                   right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                   NNarray = vecchia_approx$NNarray
                     )
                     - Bidart::log_determinant_derivatives(sparse_chol_and_grad = new_compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
                   ))) 
        )/ 2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$range_beta_sufficient ^2) / 2
      proposed_U = as.vector(
          +.5*t(c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened) %*% hierarchical_model$beta_priors$range_beta$precision_whitened %*% (c(q) - hierarchical_model$beta_priors$range_beta$mean_whitened) # normal prior 
          + .5* sum((new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
        )
      proposed_K = sum(p^2) / 2
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
        {
          state$momenta$range_beta_sufficient = p
          state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
          state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
          state$params$range_beta[] = new_range_beta
          state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          acceptance_records$range_beta_sufficient[iter - 50*(iter %/% 50) ] = acceptance_records$range_beta_sufficient[iter - 50*(iter %/% 50) ] + 1
        }
      }
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_beta_sufficient)>.91)state$transition_kernels$range_beta_sufficient = state$transition_kernels$range_beta_sufficient + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_beta_sufficient)<.61)state$transition_kernels$range_beta_sufficient = state$transition_kernels$range_beta_sufficient - rnorm(1, .4, .05)
          acceptance_records$range_beta_sufficient =  0*acceptance_records$range_beta_sufficient
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_beta_sufficient)>.91)state$transition_kernels$range_beta_sufficient = state$transition_kernels$range_beta_sufficient + rnorm(1, .2, .05)
          acceptance_records$range_beta_sufficient =  0*acceptance_records$range_beta_sufficient
        }
      }
    }
    #######################################
    # Nonstationary range beta (centered) #
    #######################################
    if(!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior))
    {
      # rotating field and beta in order to get 3 decorrelated components (each component has spatial correlation tho)
      rotated_field = state$params$range_field  %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
      rotated_beta  = matrix(state$params$range_beta, nrow = ncol(data$covariates$range_X$X_locs)) %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
      # data knowing parameter
      centered_field = as.matrix(rotated_field + data$covariates$range_X$X_locs %*% rotated_beta)
      sparse_chol_X = data$covariates$range_X$sparse_chol_X_locs
      beta_covmat = data$covariates$range_X$solve_crossprod_sparse_chol_X_locs
      # prior (whitened wrt log scale)
      rotation_matrix = t(solve(chol(Bidart::expmat(state$params$range_log_scale)))) %x% diag(1, nrow(rotated_beta), nrow(rotated_beta))
      solve_rotation_matrix = solve(rotation_matrix)
      beta_precision_prior = t(solve_rotation_matrix) %*% hierarchical_model$beta_priors$range_beta$precision %*% solve_rotation_matrix
      # combining data knowing parameter with prior 
      # Q1 + Q2
      beta_precision_posterior = diag(1, 3, 3) %x% solve(beta_covmat) + beta_precision_prior
      # (Q1 + Q2)-1 (Q1 mu1 + Q2 mu2)
      beta_mean_posterior = solve(beta_precision_posterior, 
                                  c(
                                    t(as.matrix(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol %*% centered_field))  %*% sparse_chol_X # Q1 mu1 
                                  )
                                  + 
                                    beta_precision_prior %*% rotation_matrix %*% hierarchical_model$beta_priors$range_beta$mean
      )
      
      # sampling new beta
      rotated_beta = matrix(beta_mean_posterior +  t(chol(solve(beta_precision_posterior))) %*% rnorm(length(beta_mean_posterior)), nrow = nrow(state$params$range_beta))
      # de-centering and de-rotating
      rotated_field = centered_field - as.matrix(data$covariates$range_X$X_locs %*% rotated_beta)
      state$params$range_field = rotated_field %*% chol(Bidart::expmat(state$params$range_log_scale))
      state$params$range_beta[] = rotated_beta %*% chol(Bidart::expmat(state$params$range_log_scale))
    }
    ################
    
    ################################################
    # Nonstationary range latent field (ancillary) #
    ################################################
    # MALA to update the field, using centered parametrization to update the regression coefficients (including intercept)
    if(!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior))
    {
      # Make a half step for momentum at the beginning
      q = as.matrix(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol %*% state$params$range_field) %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
      current_U =
        (
          .5 * sum(q^2) # whitenend prior 
          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        )
      # MALA whitened
      state$momenta$range_field_ancillary = sqrt(.9) * state$momenta$range_field_ancillary + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
      p = state$momenta$range_field_ancillary
      p = p - exp(state$transition_kernels$range_field_ancillary_mala) *
        (q  # white noise prior derivative
         + as.matrix(
           Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                         (# natural gradient
                           Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                         left_vector = as.vector(
                                                           Matrix::solve(
                                                             Matrix::t(state$sparse_chol_and_stuff$sparse_chol), 
                                                             - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                           ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                                             * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                           )), 
                                                         right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                         NNarray = vecchia_approx$NNarray  
                           )
                         )))
         %*% t(chol(Bidart::expmat(state$params$range_log_scale)))# whitening
        )/ 2
      # Make a full step for the position
      q = q + exp(state$transition_kernels$range_field_ancillary_mala) * p
      new_range_field = as.matrix(Matrix::solve(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol,  q)) %*%  chol(Bidart::expmat(state$params$range_log_scale))
      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                        range_beta = state$params$range_beta, NNarray = vecchia_approx$NNarray, 
                                                                        locs = data$locs, 
                                                                        range_field = new_range_field, 
                                                                        range_X = data$covariates$range_X$X_locs, 
                                                                        nu = hierarchical_model$nu
      )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
      # Make a half step for momentum at the end.
      p = p - exp(state$transition_kernels$range_field_ancillary_mala) *
        (q  # white noise prior derivative
         + as.matrix(
           Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                         (# natural gradient
                           Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                         left_vector = as.vector(
                                                           Matrix::solve(
                                                             Matrix::t(new_sparse_chol), 
                                                             - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                           ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                                             * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                           )), 
                                                         right_vector = new_field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                         NNarray = vecchia_approx$NNarray  
                           )
                         )))
         %*% t(chol(Bidart::expmat(state$params$range_log_scale)))# whitening
        )/ 2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$range_field_ancillary ^2) / 2
      proposed_U =
        (
          .5 * sum(q^2) # whitenend prior 
          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        )
      proposed_K = sum(p^2) / 2
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
        {
          state$momenta$range_field_ancillary = p
          state$params$field = new_field
          state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
          state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
          state$params$range_field = new_range_field
          state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          acceptance_records$range_field_ancillary_mala[iter - 50*(iter %/% 50) ] = acceptance_records$range_field_ancillary_mala[iter - 50*(iter %/% 50) ] + 1
        }
      }
      
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_field_ancillary_mala)>.91)state$transition_kernels$range_field_ancillary_mala = state$transition_kernels$range_field_ancillary_mala + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_field_ancillary_mala)<.61)state$transition_kernels$range_field_ancillary_mala = state$transition_kernels$range_field_ancillary_mala - rnorm(1, .4, .05)
          acceptance_records$range_field_ancillary_mala =  0*acceptance_records$range_field_ancillary_mala
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))      
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_field_ancillary_mala)<.61)state$transition_kernels$range_field_ancillary_mala = state$transition_kernels$range_field_ancillary_mala - rnorm(1, .2, .05)
          acceptance_records$range_field_ancillary_mala =  0*acceptance_records$range_field_ancillary_mala
        }
      }
    }
    ##########################################################
    # Nonstationary with log-range latent field (sufficient) #
    ##########################################################
    # MALA to update the field, using centered parametrization to update the regression coefficients (including intercept)
    if((!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior)))
    {
      q = as.matrix(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol %*% state$params$range_field) %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
      current_U =
        (
          .5 * sum(q^2) # whitened prior 
          + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
        )
      # MALA whitened
      state$momenta$range_field_sufficient = sqrt(.9) * state$momenta$range_field_sufficient + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
      p = state$momenta$range_field_sufficient
      # Make a half step for momentum at the beginning
      p = p - exp(state$transition_kernels$range_field_sufficient_mala) *
        (q  # white noise prior derivative
         + as.matrix(
           Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                         (# natural gradient
                           Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                         left_vector = as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                         right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                         NNarray = vecchia_approx$NNarray  
                           )
                           - Bidart::log_determinant_derivatives(sparse_chol_and_grad = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad,  NNarray = vecchia_approx$NNarray)# derivative of determinant
                         ))) 
         %*% t(chol(Bidart::expmat(state$params$range_log_scale)))# whitening
        ) / 2 
      # Make a full step for the position
      q = q + exp(state$transition_kernels$range_field_sufficient_mala) * p
      new_range_field = as.matrix(Matrix::solve(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol,  q)) %*%  chol(Bidart::expmat(state$params$range_log_scale))
      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                        range_beta = state$params$range_beta, NNarray = vecchia_approx$NNarray, 
                                                                        locs = data$locs, 
                                                                        range_field = new_range_field, 
                                                                        range_X = data$covariates$range_X$X_locs, 
                                                                        nu = hierarchical_model$nu
      )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      # Make a half step for momentum at the end.
      p = p - exp(state$transition_kernels$range_field_sufficient_mala) *
        (q  # white noise prior derivative
         + as.matrix(
           Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                         (# natural gradient
                           Bidart::derivative_sandwiches(derivative = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                         left_vector = as.vector(new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                         right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                         NNarray = vecchia_approx$NNarray
                           )
                           - Bidart::log_determinant_derivatives(new_compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
                         ))) 
         %*% t(chol(Bidart::expmat(state$params$range_log_scale)))# whitening
        )/ 2 
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$range_field_sufficient ^2) / 2
      proposed_U =
        (
          .5 * sum(q^2) # whitenend prior 
          + .5* sum((new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
        )
      proposed_K = sum(p^2) / 2
      if(!is.nan(current_U-proposed_U+current_K- proposed_K))
      {
        if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
        {
          state$momenta$range_field_sufficient = p
          state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
          state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
          state$params$range_field = new_range_field
          state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
          acceptance_records$range_field_sufficient_mala[iter - 50*(iter %/% 50) ] = acceptance_records$range_field_sufficient_mala[iter - 50*(iter %/% 50) ] + 1
        }
      }
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_field_sufficient_mala)>.91)state$transition_kernels$range_field_sufficient_mala = state$transition_kernels$range_field_sufficient_mala + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_field_sufficient_mala)<.61)state$transition_kernels$range_field_sufficient_mala = state$transition_kernels$range_field_sufficient_mala - rnorm(1, .4, .05)
          acceptance_records$range_field_sufficient_mala =  0*acceptance_records$range_field_sufficient_mala
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_field_sufficient_mala)<.61)state$transition_kernels$range_field_sufficient_mala = state$transition_kernels$range_field_sufficient_mala - rnorm(1, .2, .05)
          acceptance_records$range_field_sufficient_mala =  0*acceptance_records$range_field_sufficient_mala
        }
      }
    }
    ####################################################
    # Variance of the nonstationary range latent field #
    ####################################################
    if(!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior)){
      
      ###      # ancillary-sufficient
      ###      new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, exp(.5 * state$transition_kernels$range_log_scale_sufficient))
      ###      new_range_field = state$params$range_field %*% solve(chol(Bidart::expmat(state$params$range_log_scale))) %*% chol(Bidart::expmat(new_range_log_scale))
      ###      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
      ###                                                                range_beta = state$params$range_beta, NNarray = vecchia_approx$NNarray, 
      ###                                                                locs = data$locs, 
      ###                                                                range_field = new_range_field, 
      ###                                                                range_X = data$covariates$range_X$X_locs, nu = hierarchical_model$nu)
      ###      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      ###      if(
      ###        + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
      ###        - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
      ###        - .5* sum((new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
      ###        + sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
      ###        > log(runif(1))
      ###      )
      ###      {
      ###        state$params$range_log_scale = new_range_log_scale
      ###        state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
      ###        state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
      ###        state$params$range_field = new_range_field
      ###        state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
      ###        acceptance_records$range_log_scale_sufficient[iter - 50*(iter %/% 50) ] = acceptance_records$range_log_scale_sufficient[iter - 50*(iter %/% 50) ] + 1
      ###      }
      ###      # updating MALA kernel
      ###      if(iter_start + iter < 1000)
      ###      {
      ###        if(iter %/% 50 ==iter / 50)
      ###        {
      ###          if(mean(acceptance_records$range_log_scale_sufficient)>.41  & state$transition_kernels$range_log_scale_sufficient< (0))state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient + rnorm(1, .4, .05)
      ###          if(mean(acceptance_records$range_log_scale_sufficient)<.11)state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient - rnorm(1, .4, .05)
      ###          acceptance_records$range_log_scale_sufficient =  0*acceptance_records$range_log_scale_sufficient
      ###        }
      ###      }
      # ancillary-sufficient
      new_range_log_scale = state$params$range_log_scale
      new_range_field = state$params$range_field
      current_U =
        (
          0  # improper prior
          + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
        )
      # HMC
      state$momenta$range_log_scale_sufficient = sqrt(.9) * state$momenta$range_log_scale_sufficient + sqrt(.1)*rnorm(length(new_range_log_scale))
      p = state$momenta$range_log_scale_sufficient
      # Make a half step for momentum at the beginning
      # derivative of range field wrt log scale
      d_field = Bidart::derivative_field_wrt_scale(as.matrix(new_range_field), new_range_log_scale)
      # derivative of potential  wrt range field
      d_potential = 
        (
          Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                        left_vector = as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                        right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                        NNarray = vecchia_approx$NNarray 
          )
          - Bidart::log_determinant_derivatives(sparse_chol_and_grad = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad,  NNarray = vecchia_approx$NNarray)# derivative of determinant
        ) 
      # chain rule
      p = p - exp(state$transition_kernels$range_log_scale_sufficient) *
        (
          apply(d_field * array(rep(d_potential, length(p)), dim = c(dim(d_potential), length(p))), 3, sum)
          + 0 # improper prior
        )/2
      # update the position
      new_range_log_scale = new_range_log_scale + exp(state$transition_kernels$range_log_scale_sufficient) * p
      if(
        (min(eigen(Bidart::expmat(new_range_log_scale))$values)>exp(-8))&
        (max(eigen(Bidart::expmat(new_range_log_scale))$values)<exp(3))
      )
      {
        new_range_field = state$params$range_field %*% solve(chol(Bidart::expmat(state$params$range_log_scale))) %*% chol(Bidart::expmat(new_range_log_scale))
        new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                          range_beta = state$params$range_beta, NNarray = vecchia_approx$NNarray, 
                                                                          locs = data$locs, 
                                                                          range_field = new_range_field, 
                                                                          range_X = data$covariates$range_X$X_locs, 
                                                                          nu = hierarchical_model$nu)
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        # Make a half step for momentum at the end
        # derivative of range field wrt log scale
        d_field = Bidart::derivative_field_wrt_scale(new_range_field, new_range_log_scale)
        # derivative of potential  wrt range field
        d_potential = 
          (# natural gradient
            Bidart::derivative_sandwiches(derivative = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                          left_vector = as.vector(new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                          right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                          NNarray = vecchia_approx$NNarray
            )
            - Bidart::log_determinant_derivatives(new_compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
          )
        # chain rule
        p = p - exp(state$transition_kernels$range_log_scale_sufficient) *
          (
            apply(d_field * array(rep(d_potential, length(p)), dim = c(dim(d_potential), length(p))), 3, sum)
            + 0 # improper prior
          )/2
        proposed_U =
          (
            0  # improper prior
            + .5* sum((new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
            - sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
          )
        current_K = sum(state$momenta$range_log_scale_sufficient^2)
        proposed_K = sum(state$momenta$p^2)
        
        if(!is.nan(current_U-proposed_U+current_K- proposed_K))
        {
          if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
          {
            state$momenta$range_log_scale_sufficient = p
            state$params$range_log_scale = new_range_log_scale
            state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            state$params$range_field = new_range_field
            state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
            acceptance_records$range_log_scale_sufficient[iter - 50*(iter %/% 50) ] = acceptance_records$range_log_scale_sufficient[iter - 50*(iter %/% 50) ] + 1
          }
        }
      }
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_log_scale_sufficient)>.91)state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_log_scale_sufficient)<.61)state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient - rnorm(1, .4, .05)
          acceptance_records$range_log_scale_sufficient =  0*acceptance_records$range_log_scale_sufficient
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_log_scale_sufficient)<.61)state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient - rnorm(1, .2, .05)
          acceptance_records$range_log_scale_sufficient =  0*acceptance_records$range_log_scale_sufficient
        }
      }
      
      # sufficient - sufficient
      whitened_range_field_crossprod = crossprod(as.matrix(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol %*% state$params$range_field))
      for(i in seq(5))
      {
        new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, .05)
        if(
          (min(eigen(Bidart::expmat(new_range_log_scale))$values)>exp(-8))&
          (max(eigen(Bidart::expmat(new_range_log_scale))$values)<exp(3))
        )
        {
          if(
            -.5 * vecchia_approx$n_locs * (determinant(Bidart::expmat(new_range_log_scale))$mod - determinant(Bidart::expmat(state$params$range_log_scale))$mod)
            -.5 * sum(
              diag(
                (t(solve(chol(Bidart::expmat(new_range_log_scale)))) %*% whitened_range_field_crossprod %*% solve(chol(Bidart::expmat(new_range_log_scale)))) -
                (t(solve(chol(Bidart::expmat(state$params$range_log_scale)))) %*% whitened_range_field_crossprod %*% solve(chol(Bidart::expmat(state$params$range_log_scale))))
              )
            )
            > log(runif(1))
          )
          {
            state$params$range_log_scale = new_range_log_scale
          }
        }
      }
      
      #      # ancillary - ancillary
      #      new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, exp(.5 * state$transition_kernels$range_log_scale_ancillary))
      #      new_range_field = state$params$range_field %*% solve(chol(Bidart::expmat(state$params$range_log_scale))) %*% chol(Bidart::expmat(new_range_log_scale))
      #      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
      #                                                                        range_beta = state$params$range_beta, NNarray = vecchia_approx$NNarray, 
      #                                                                        locs = data$locs, 
      #                                                                        range_field = new_range_field, 
      #                                                                        range_X = data$covariates$range_X$X_locs, nu = hierarchical_model$nu)
      #      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      #      new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
      #      if(
      #        - .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
      #        + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
      #        > log(runif(1))
      #      )
      #      {
      #        state$params$range_log_scale = new_range_log_scale
      #        state$params$field = new_field
      #        state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
      #        state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
      #        state$params$range_field = new_range_field
      #        state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
      #        acceptance_records$range_log_scale_ancillary[iter - 50*(iter %/% 50) ] = acceptance_records$range_log_scale_ancillary[iter - 50*(iter %/% 50) ] + 1
      #      }
      #      
      #      # updating MALA kernel
      #      if(iter_start + iter < 1000)
      #      {
      #        if(iter %/% 50 ==iter / 50)
      #        {
      #          if((mean(acceptance_records$range_log_scale_ancillary)>.41) & state$transition_kernels$range_log_scale_ancillary< (0))state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary + rnorm(1, .4, .05)
      #          if(mean(acceptance_records$range_log_scale_ancillary)<.11)state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary - rnorm(1, .4, .05)
      #          acceptance_records$range_log_scale_ancillary =  0*acceptance_records$range_log_scale_ancillary
      #        }
      #      }
      # ancillary - ancillary
      new_range_log_scale = state$params$range_log_scale
      new_range_field = state$params$range_field
      new_field = state$params$field
      current_U =
        (
          0 # improper prior 
          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        )
      # HMC
      state$momenta$range_log_scale_ancillary = sqrt(.6) * state$momenta$range_log_scale_ancillary + sqrt(.4)*rnorm(length(new_range_log_scale))
      p = state$momenta$range_log_scale_ancillary
      # Make a half step for momentum at the beginning
      # derivative of range field wrt log scale
      d_field = Bidart::derivative_field_wrt_scale(as.matrix(new_range_field), new_range_log_scale)
      # derivative of potential  wrt range field
      d_potential = 
        Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                      left_vector = as.vector(
                                        Matrix::solve(
                                          Matrix::t(new_sparse_chol), 
                                          - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                        ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                          * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                        )), 
                                      right_vector = new_field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                      NNarray = vecchia_approx$NNarray  
        )
      # chain rule
      p = p - exp(state$transition_kernels$range_log_scale_ancillary) *
        (
          apply(d_field * array(rep(d_potential, length(p)), dim = c(dim(d_potential), length(p))), 3, sum)
          + 0 # improper prior
        )/2
      # updating the position
      new_range_log_scale = new_range_log_scale + exp(state$transition_kernels$range_log_scale_ancillary) * p
      if(
        (min(eigen(Bidart::expmat(new_range_log_scale))$values)>exp(-8))&
        (max(eigen(Bidart::expmat(new_range_log_scale))$values)<exp(3))
      )
      {
        new_range_field = state$params$range_field %*% solve(chol(Bidart::expmat(state$params$range_log_scale))) %*% chol(Bidart::expmat(new_range_log_scale))
        new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                          range_beta = state$params$range_beta, NNarray = vecchia_approx$NNarray, 
                                                                          locs = data$locs, 
                                                                          range_field = new_range_field, 
                                                                          range_X = data$covariates$range_X$X_locs, 
                                                                          nu = hierarchical_model$nu)
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
        # Make a half step for momentum at the end
        # derivative of range field wrt log scale
        d_field = Bidart::derivative_field_wrt_scale(as.matrix(new_range_field), new_range_log_scale)
        # derivative of potential  wrt range field
        d_potential = 
          Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                        left_vector = as.vector(
                                          Matrix::solve(
                                            Matrix::t(new_sparse_chol), 
                                            - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                          ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                            * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                          )), 
                                        right_vector = new_field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                        NNarray = vecchia_approx$NNarray  
          )
        # chain rule
        p = p - exp(state$transition_kernels$range_log_scale_ancillary) *
          (
            apply(d_field * array(rep(d_potential, length(p)), dim = c(dim(d_potential), length(p))), 3, sum)
            + 0 # improper prior
          )/2
        
        current_K = sum (state$momenta$range_log_scale_ancillary ^2) / 2
        proposed_U =
          (
            0 # improper prior 
            + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
          )
        proposed_K = sum(p^2) / 2
        if(!is.nan(current_U-proposed_U+current_K- proposed_K))
        {
          if (log(runif(1)) < current_U-proposed_U + current_K- proposed_K)
          {
            state$momenta$range_log_scale_ancillary = p
            state$params$range_log_scale = new_range_log_scale
            state$params$field = new_field
            state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
            state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
            state$params$range_field = new_range_field
            state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
            acceptance_records$range_log_scale_ancillary[iter - 50*(iter %/% 50) ] = acceptance_records$range_log_scale_ancillary[iter - 50*(iter %/% 50) ] + 1
          }
        }
      }
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if((mean(acceptance_records$range_log_scale_ancillary)>.91) )state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_log_scale_ancillary)<.61)state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary - rnorm(1, .4, .05)
          acceptance_records$range_log_scale_ancillary =  0*acceptance_records$range_log_scale_ancillary
        }
      }
      
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))      
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_log_scale_ancillary)<.61)state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary - rnorm(1, .2, .05)
          acceptance_records$range_log_scale_ancillary =  0*acceptance_records$range_log_scale_ancillary
        }
      }
      
      # sufficient - sufficient
      whitened_range_field_crossprod = crossprod(as.matrix(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol %*% state$params$range_field))
      for(i in seq(5))
      {
        new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, .05)
        if(
          (min(eigen(Bidart::expmat(new_range_log_scale))$values)>exp(-8))&
          (max(eigen(Bidart::expmat(new_range_log_scale))$values)<exp(3))
        )
        {
          if(
            -.5 * vecchia_approx$n_locs * (determinant(Bidart::expmat(new_range_log_scale))$mod - determinant(Bidart::expmat(state$params$range_log_scale))$mod)
            -.5 * sum(
              diag(
                (t(solve(chol(Bidart::expmat(new_range_log_scale)))) %*% whitened_range_field_crossprod %*% solve(chol(Bidart::expmat(new_range_log_scale)))) -
                (t(solve(chol(Bidart::expmat(state$params$range_log_scale)))) %*% whitened_range_field_crossprod %*% solve(chol(Bidart::expmat(state$params$range_log_scale))))
              )
            )
            > log(runif(1))
          )
          {
            state$params$range_log_scale = new_range_log_scale
          }
        }
      }
    }
    #########
    # Noise #
    #########
    # recomputation in order to avoid errors
    state$sparse_chol_and_stuff$noise = Bidart::variance_field(beta = state$params$noise_beta, X = data$covariates$noise_X$X, field = state$params$noise_field[vecchia_approx$locs_match])
    # VEWY IMPOWTANT don't remove or comment
    squared_residuals = as.vector(state$sparse_chol_and_stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
    ##############
    # Noise beta #
    ##############
    # HMC update
    q = t(solve(data$covariates$noise_X$chol_solve_crossprod_X)) %*% state$params$noise_beta
    current_U =
      (
        +.5*t(q - hierarchical_model$beta_priors$noise_beta$mean_whitened) %*% hierarchical_model$beta_priors$noise_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$noise_beta$mean_whitened) # normal prior 
        +.5* sum(log(state$sparse_chol_and_stuff$noise)) # det
        +.5*sum(squared_residuals/state$sparse_chol_and_stuff$noise) # observations
      )
    # HMC whitened
    state$momenta$noise_beta = sqrt(.9) * state$momenta$noise_beta + sqrt(.1)*rnorm(length(state$momenta$noise_beta))
    p = state$momenta$noise_beta
    
    # Make a half step for momentum at the beginning
    p = p - exp(state$transition_kernels$noise_beta_mala) *
      (
        + 0.5* hierarchical_model$beta_priors$noise_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$noise_beta$mean_whitened) 
        + solve(solve(data$covariates$noise_X$chol_solve_crossprod_X), # solving by prior sparse chol because of whitening
              t(data$covariates$noise_X$X) %*%
                (
                  .5 # determinant part of normal likelihood
                  - (squared_residuals/state$sparse_chol_and_stuff$noise)/2 # exponential part of normal likelihood
                )
      ))/ 2
    # checking gradient with finite differences
    ###noise_beta_ = state$params$noise_beta 
    ###noise_beta_[1] = noise_beta_[1] + 0.0001
    ###noise_ = Bidart::variance_field(beta = noise_beta_, field = state$params$noise_field[vecchia_approx$locs_match], data$covariates$noise_X$X)
    ###U_ =
    ###  (
    ###    0 # improper prior 
    ###    +.5* sum(log(noise_)) # det
    ###    +.5*sum(squared_residuals/noise_) # observations
    ###  )
    ###print(10000*(U_- current_U))
    ###print(
    ###  (t(data$covariates$noise_X$X) %*%
    ###    (
    ###      .5 # determinant part of normal likelihood
    ###      - (squared_residuals/state$sparse_chol_and_stuff$noise)/2 # exponential part of normal likelihood
    ###    ))[1]
    ###)
    nsteps = 4 + rbinom(1, 1, .5)
    for(i in seq(nsteps))
    {
      # Make a full step for the position
      q = q + exp(state$transition_kernels$noise_beta_mala) * p
      new_noise_beta = t(data$covariates$noise_X$chol_solve_crossprod_X) %*% q
      new_noise = Bidart::variance_field(new_noise_beta, state$params$noise_field[vecchia_approx$locs_match], data$covariates$noise_X$X)
      # Make a full step for momentum at the end
      if(i!= nsteps)
      {
        p = p - exp(state$transition_kernels$noise_beta_mala) *
          (
            + 0.5* hierarchical_model$beta_priors$noise_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$noise_beta$mean_whitened) 
            + solve(solve(data$covariates$noise_X$chol_solve_crossprod_X), # solving by prior sparse chol because of whitening
                  t(data$covariates$noise_X$X) %*%
                    (
                      .5 # determinant part of normal likelihood
                      - (squared_residuals/new_noise)/2 # exponential part of normal likelihood
                    )
          ))
      }
    }
    # Make a half step for momentum at the end
    p = p - exp(state$transition_kernels$noise_beta_mala) *
      (
        + 0.5* hierarchical_model$beta_priors$noise_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$noise_beta$mean_whitened) 
        + solve(solve(data$covariates$noise_X$chol_solve_crossprod_X), # solving by prior sparse chol because of whitening
              t(data$covariates$noise_X$X) %*%
                (
                  .5 # determinant part of normal likelihood
                  - (squared_residuals/new_noise)/2 # exponential part of normal likelihood
                )
      ))/ 2
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (state$momenta$noise_beta ^2) / 2
    proposed_U = 
      (
        +.5*t(q - hierarchical_model$beta_priors$noise_beta$mean_whitened) %*% hierarchical_model$beta_priors$noise_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$noise_beta$mean_whitened) # normal prior 
        +.5* sum(log(new_noise)) # det
        +.5*sum(squared_residuals/new_noise) # observations
      )
    proposed_K = sum(p^2) / 2
    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
    {
      if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
      {
        state$momenta$noise_beta = p
        state$params$noise_beta[] = new_noise_beta
        state$sparse_chol_and_stuff$noise = new_noise
        acceptance_records$noise_beta_mala[iter - 50*(iter %/% 50) ] = acceptance_records$noise_beta_mala[iter - 50*(iter %/% 50) ] + 1
      }
    }
    # updating MALA kernel
    if(iter_start + iter < 1000)
    {
      if(iter %/% 50 ==iter / 50)
      {
        if(mean(acceptance_records$noise_beta_mala)>.91)state$transition_kernels$noise_beta_mala = state$transition_kernels$noise_beta_mala +rnorm(1, .4, .05)
        if(mean(acceptance_records$noise_beta_mala)<.61)state$transition_kernels$noise_beta_mala = state$transition_kernels$noise_beta_mala -rnorm(1, .4, .05)
        acceptance_records$noise_beta_mala =  0*acceptance_records$noise_beta_mala
      }
    }
    if((iter_start + iter > 1000)&(iter_start + iter < 2000))    
    {
      if(iter %/% 50 ==iter / 50)
      {
        if(mean(acceptance_records$noise_beta_mala)<.61)state$transition_kernels$noise_beta_mala = state$transition_kernels$noise_beta_mala -rnorm(1, .2, .05)
        acceptance_records$noise_beta_mala =  0*acceptance_records$noise_beta_mala
      }
    }
    
    ########################
    # Noise beta and field #
    ########################
    if(!is.null(hierarchical_model$hyperprior_covariance$noise_NNGP_prior))
    {
      #######################################################
      # Gaussian update of centered regression coefficients #
      #######################################################
      #centered_field = as.vector(state$params$noise_field + data$covariates$noise_X$X_locs %*% matrix(state$params$noise_beta[data$covariates$noise_X$which_locs], ncol = 1))
      #sparse_chol_X = data$covariates$noise_X$sparse_chol_X_locs
      #beta_covmat = data$covariates$noise_X$solve_crossprod_sparse_chol_X_locs
      #beta_mean =  c(as.vector(hierarchical_model$hyperprior_covariance$noise_NNGP_prior$sparse_chol %*% centered_field)  %*% sparse_chol_X %*% beta_covmat)
      #state$params$noise_beta [data$covariates$noise_X$which_locs] = as.vector(beta_mean + exp(.5 * state$params$noise_log_scale) * t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      #state$params$noise_field = centered_field - as.vector(data$covariates$noise_X$X_locs %*% matrix(state$params$noise_beta[data$covariates$noise_X$which_locs], ncol = 1))
      new_beta_and_field = Bidart::update_centered_beta(
        beta = state$params$noise_beta, 
        field = state$params$noise_field,
        X = data$covariates$noise_X,
        beta_prior = hierarchical_model$beta_priors$noise_beta, 
        NNGP_prior = hierarchical_model$hyperprior_covariance$noise_NNGP_prior, 
        log_scale = state$params$noise_log_scale
      )
      state$params$noise_beta[data$noise_X$which_locs] = new_beta_and_field$beta
      state$params$noise_field[] = new_beta_and_field$field
      
      
      ######################################################################
      # HMC update of the latent field with whitening & momentum recycling #
      ######################################################################
      if(iter_start+iter>50) # in some cases waiting for the fixed effects to adjust avoids HMC crashing
      {
        weighted_squared_residuals = as.vector(vecchia_approx$locs_match_matrix%*%(squared_residuals/state$sparse_chol_and_stuff$noise))
        state$momenta$noise_field = sqrt(.9) * state$momenta$noise_field + sqrt(.1)*rnorm(vecchia_approx$n_locs)
        p = state$momenta$noise_field
        q = exp(-.5*state$params$noise_log_scale)*as.vector(hierarchical_model$hyperprior_covariance$noise_NNGP_prior$sparse_chol %*% state$params$noise_field)
        current_U =
          (
            .5 * sum(q^2) # whitenend prior 
            + .5* sum(log(state$sparse_chol_and_stuff$noise)) + .5*sum(weighted_squared_residuals)
          )
        # Make a half step for momentum at the beginning
        p = p - exp(state$transition_kernels$noise_field_mala) *
          (q  # white noise prior derivative
           + exp(.5*state$params$noise_log_scale)*as.vector(
             Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$noise_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                           (.5 * vecchia_approx$obs_per_loc # determinant part of normal likelihood
                            - #state$params$noise_field  * 
                              weighted_squared_residuals/2 # exponential part of normal likelihood
                           ))) # natural derivative
          )/ 2
        for(i in seq(5))
        {
          # Make a full step for the position
          q = q + exp(state$transition_kernels$noise_field_mala) * p
          new_noise_field = exp(.5*state$params$noise_log_scale)*as.vector(Matrix::solve(hierarchical_model$hyperprior_covariance$noise_NNGP_prior$sparse_chol, q)) 
          new_noise = state$sparse_chol_and_stuff$noise * (exp(new_noise_field - state$params$noise_field)[vecchia_approx$locs_match])
          new_weighted_squared_residuals = as.vector(vecchia_approx$locs_match_matrix%*%(squared_residuals/new_noise))
          # Make a full step for momentum
          if(i != 5)
          {
            p = p - exp(state$transition_kernels$noise_field_mala) *
              (q  # white noise prior derivative
               + exp(.5*state$params$noise_log_scale)*as.vector(Matrix::solve(
                 Matrix::t(hierarchical_model$hyperprior_covariance$noise_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                 (.5 * vecchia_approx$obs_per_loc - 
                    #new_noise_field  * 
                    new_weighted_squared_residuals/2) # natural derivative
               ))
              )
          }
        }
        # Make a half step for momentum at the end.
        p = p - exp(state$transition_kernels$noise_field_mala) *
          (q  # white noise prior derivative
           + exp(.5*state$params$noise_log_scale)*as.vector(Matrix::solve(
             Matrix::t(hierarchical_model$hyperprior_covariance$noise_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
             (.5 * vecchia_approx$obs_per_loc - 
                #new_noise_field  * 
                new_weighted_squared_residuals/2) # natural derivative
           ))
          )/ 2
        # Evaluate potential and kinetic energies at start and end of trajectory
        current_K = sum (state$momenta$noise_field ^2) / 2
        proposed_U = 
          (
            .5 * sum(q^2) 
            + .5* sum(log(new_noise)) + .5*sum(new_weighted_squared_residuals)
          )
        proposed_K = sum(p^2) / 2
        if(!is.nan(current_U-proposed_U+current_K- proposed_K))
        {
          if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
          {
            state$params$noise_field = new_noise_field
            acceptance_records$noise_field_mala[iter - 50*(iter %/% 50) ] = acceptance_records$noise_field_mala[iter - 50*(iter %/% 50) ] + 1
            state$sparse_chol_and_stuff$noise = new_noise
            weighted_squared_residuals = new_weighted_squared_residuals
            state$momenta$noise_field = p
          }
        }
        # updating MALA kernel
        if(iter_start + iter < 1000)
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$noise_field_mala)>.91)state$transition_kernels$noise_field_mala= state$transition_kernels$noise_field_mala + rnorm(1, .4, .05)
            if(mean(acceptance_records$noise_field_mala)<.61)state$transition_kernels$noise_field_mala = state$transition_kernels$noise_field_mala - rnorm(1, .4, .05)
            acceptance_records$noise_field_mala =  0*acceptance_records$noise_field_mala
          }
        }
        if((iter_start + iter > 1000)&(iter_start + iter < 2000))
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$noise_field_mala)<.61)state$transition_kernels$noise_field_mala = state$transition_kernels$noise_field_mala - rnorm(1, .2, .05)
            acceptance_records$noise_field_mala =  0*acceptance_records$noise_field_mala
          }
        }
        
      }
      ####################################################
      # Variance of the nonstationary noise latent field #
      ####################################################
      # ancillary -- sufficient
      new_noise_log_scale = state$params$noise_log_scale + rnorm(1, 0, exp(state$transition_kernels$noise_log_scale))
      new_noise_field = state$params$noise_field * exp((new_noise_log_scale - state$params$noise_log_scale)/2)
      new_noise = state$sparse_chol_and_stuff$noise * exp(new_noise_field[vecchia_approx$locs_match] - state$params$noise_field[vecchia_approx$locs_match])
      if(
        (
          -.5* sum(log(new_noise)) 
          -.5*sum(squared_residuals/new_noise)
          +.5* sum(log(state$sparse_chol_and_stuff$noise)) 
          +.5*sum(squared_residuals/state$sparse_chol_and_stuff$noise)
          > log(runif(1)
          )
          & (new_noise_log_scale > -8)
          & (new_noise_log_scale < 3)
        )
      )
      {
        state$params$noise_log_scale = new_noise_log_scale
        state$params$noise_field = new_noise_field 
        state$sparse_chol_and_stuff$noise = new_noise
        acceptance_records$noise_log_scale[iter - 50*(iter %/% 50) ] = acceptance_records$noise_log_scale[iter - 50*(iter %/% 50) ]+1
      }
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$noise_log_scale)>.41)state$transition_kernels$noise_log_scale = state$transition_kernels$noise_log_scale + rnorm(1, .4, .05)
          if(mean(acceptance_records$noise_log_scale)<.11)state$transition_kernels$noise_log_scale = state$transition_kernels$noise_log_scale - rnorm(1, .4, .05)
          acceptance_records$noise_log_scale =  0*acceptance_records$noise_log_scale
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$noise_log_scale)<.11)state$transition_kernels$noise_log_scale = state$transition_kernels$noise_log_scale - rnorm(1, .2, .05)
          acceptance_records$noise_log_scale =  0*acceptance_records$noise_log_scale
        }
      }
      # sufficient -- sufficient
      unscaled_thingy = sum((hierarchical_model$hyperprior_covariance$noise_NNGP_prior$sparse_chol %*% state$params$noise_field)^2)
      
      
      for(i in seq(4))
      {
        new_noise_log_scale = state$params$noise_log_scale + rnorm(1, 0, .1)
        if(
          (
            -.5*vecchia_approx$n_locs*new_noise_log_scale -.5* exp(-new_noise_log_scale)*unscaled_thingy
            +.5*vecchia_approx$n_locs*state$params$noise_log_scale -.5* exp(-state$params$noise_log_scale)*unscaled_thingy
          ) > log(runif(1))
          & (new_noise_log_scale > -8)
          & (new_noise_log_scale < 3)        ){
          state$params$noise_log_scale = new_noise_log_scale
        }
      }
    }
    
    
    #########
    # Scale #
    #########
    ##############
    # Scale beta #
    ##############
      ##############################
      # ancillary - sufficient HMC #
      ##############################
    sparse_chol_diag_field = state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = state$params$field)
    q = t(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs)) %*% state$params$scale_beta
    current_U =
      (
        +0.5*t(q - hierarchical_model$beta_priors$scale_beta$mean_whitened) %*% hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
        +0.5*sum(log(state$sparse_chol_and_stuff$scale))# determinant part
        +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/state$sparse_chol_and_stuff$scale))^2)# covmat product part
      )
    # HMC whitened
    state$momenta$scale_beta_sufficient = sqrt(.9) * state$momenta$scale_beta_sufficient + sqrt(.1)*rnorm(ncol(data$covariates$scale_X$X_locs))
    p = state$momenta$scale_beta_sufficient
    # Make a half step for momentum at the beginning
    p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
      (
        + 0.5* hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
        + solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
              t(data$covariates$scale_X$X_locs) %*%
                (
                  .5  # determinant part 
                  -.5 * sqrt(1/state$sparse_chol_and_stuff$scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/state$sparse_chol_and_stuff$scale)))# natural derivative
                )
      ))/ 2
    # testing gradient vs finite difference
    ###scale_beta_ = state$params$scale_beta 
    ###scale_beta_[1] = scale_beta_[1] + .001
    ###scale_ = Bidart::variance_field(beta = scale_beta_, X = data$covariates$scale_X$X_locs, field = state$params$scale_field)
    ###finite_diff_grad = 
    ###  (
    ###    +0.5*sum(log(scale_))# determinant part
    ###    +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/scale_))^2)# covmat product part
    ###    -0.5*sum(log(state$sparse_chol_and_stuff$scale))# determinant part
    ###    -0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/state$sparse_chol_and_stuff$scale))^2)# covmat product part
    ###  )/.001
    ###print(finite_diff_grad)
    ###print(
    ###  (t(data$covariates$scale_X$X_locs) %*%
    ###    (
    ###      .5  # determinant part 
    ###      -.5 * sqrt(1/state$sparse_chol_and_stuff$scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/state$sparse_chol_and_stuff$scale)))# natural derivative
    ###      ))[1]
    ###  )
    n_hmc_steps = 4 + rbinom(1, 1, .5)
    for(i in seq(n_hmc_steps))
    {
      # Make a full step for the position
      q = q + exp(state$transition_kernels$scale_beta_sufficient_mala) * p
      new_scale_beta = t(data$covariates$scale_X$chol_solve_crossprod_X_locs) %*% q
      new_scale = Bidart::variance_field(beta = new_scale_beta, X = data$covariates$scale_X$X_locs, field = state$params$scale_field)
      if(i != n_hmc_steps)
      {
        p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
          (
            + 0.5* hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
            + solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  t(data$covariates$scale_X$X_locs) %*%
                    (
                      .5  # determinant part 
                      -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
                    )
          ))
      }
    }
    # Make a half step for momentum at the end.
    p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
      (
        + 0.5* hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
        + solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
              t(data$covariates$scale_X$X_locs) %*%
                (
                  .5  # determinant part 
                  -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
                )
      ))/ 2
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (state$momenta$scale_beta_sufficient ^2) / 2
    proposed_U = 
      (
        +0.5*t(q - hierarchical_model$beta_priors$scale_beta$mean_whitened) %*% hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
        +0.5*sum(log(new_scale))# determinant part
        +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/new_scale))^2)# covmat product part
      )
    proposed_K = sum(p^2) / 2
    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
    {
      if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
      {
        state$momenta$scale_beta_sufficient = p
        state$params$scale_beta[] = new_scale_beta
        acceptance_records$scale_beta_sufficient_mala[iter - 50*(iter %/% 50) ] = acceptance_records$scale_beta_sufficient_mala[iter - 50*(iter %/% 50) ] + 1
        state$sparse_chol_and_stuff$scale = new_scale
      }
    }
    # updating MALA kernel
    if(iter_start + iter < 1000)
    {
      if(iter %/% 50 ==iter / 50)
      {
        if(mean(acceptance_records$scale_beta_sufficient_mala)>.91)state$transition_kernels$scale_beta_sufficient_mala = state$transition_kernels$scale_beta_sufficient_mala +rnorm(1, .4, .05)
        if(mean(acceptance_records$scale_beta_sufficient_mala)<.61)state$transition_kernels$scale_beta_sufficient_mala = state$transition_kernels$scale_beta_sufficient_mala -rnorm(1, .4, .05)
        acceptance_records$scale_beta_sufficient_mala =  0*acceptance_records$scale_beta_sufficient_mala
      }
    }
    if((iter_start + iter > 1000)&(iter_start + iter < 2000))    
    {
      if(iter %/% 50 ==iter / 50)
      {
        if(mean(acceptance_records$scale_beta_sufficient_mala)<.61)state$transition_kernels$scale_beta_sufficient_mala = state$transition_kernels$scale_beta_sufficient_mala -rnorm(1, .2, .05)
        acceptance_records$scale_beta_sufficient_mala =  0*acceptance_records$scale_beta_sufficient_mala
      }
    }
      ######################################
      # sufficient - sufficient analytical #
      ######################################
    if(!is.null(hierarchical_model$hyperprior_covariance$scale_NNGP_prior))
    {
      new_beta_and_field = Bidart::update_centered_beta(
        beta = state$params$scale_beta, 
        field = state$params$scale_field,
        X = data$covariates$scale_X,
        beta_prior = hierarchical_model$beta_priors$scale_beta, 
        NNGP_prior = hierarchical_model$hyperprior_covariance$scale_NNGP_prior, 
        log_scale = state$params$scale_log_scale
      )
      state$params$scale_beta[] = new_beta_and_field$beta
      state$params$scale_field[] = new_beta_and_field$field
    }
    
      #############################
      # ancillary - ancillary HMC #
      #############################
    q = t(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs)) %*% state$params$scale_beta
    current_U =
      (
        +0.5*t(q - hierarchical_model$beta_priors$scale_beta$mean_whitened) %*% hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
        +.5 * sum((data$observed_field - state$sparse_chol_and_stuff$lm_fit - state$params$field[vecchia_approx$locs_match] )^2/state$sparse_chol_and_stuff$noise)
      )
    # MALA whitened
    state$momenta$scale_beta_ancillary = sqrt(.9) * state$momenta$scale_beta_ancillary   + sqrt(.1)*rnorm(ncol(data$covariates$scale_X$X_locs))
    p = state$momenta$scale_beta_ancillary
    # Make a. half step for momentum at the beginning
    p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
      (
        + 0.5* hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
        + solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
              t(data$covariates$scale_X$X_locs) %*%
                (.5 * state$params$field * 
                   as.vector(vecchia_approx$locs_match_matrix %*% 
                               ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                  state$sparse_chol_and_stuff$noise)
                   ))
      ))/ 2
    n_hmc_steps = 4 + rbinom(1, 1, .5)
    for(i in seq(n_hmc_steps))
    {
      # Make a full step for the position
      q = q + exp(state$transition_kernels$scale_beta_ancillary_mala) * p
      new_scale_beta = t(data$covariates$scale_X$chol_solve_crossprod_X_locs) %*% q
      new_scale = Bidart::variance_field(beta = new_scale_beta, X = data$covariates$scale_X$X_locs, field = state$params$scale_field)
      new_field  = state$params$field*sqrt(new_scale)/sqrt(state$sparse_chol_and_stuff$scale)
      # Make a full step for momentum.
      if(i != n_hmc_steps)
      {
        p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
          (
            + 0.5* hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
            + solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                  t(data$covariates$scale_X$X_locs) %*%
                    (.5 * new_field * 
                       as.vector(vecchia_approx$locs_match_matrix %*% 
                                   ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                      state$sparse_chol_and_stuff$noise)
                       )))
          )
      }
    }
    # Make a half step for momentum at the end.
    p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
      (
        + 0.5* hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
        + solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
              t(data$covariates$scale_X$X_locs) %*%
                (.5 * new_field * 
                   as.vector(vecchia_approx$locs_match_matrix %*% 
                               ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                  state$sparse_chol_and_stuff$noise)
                   ))
      ))/ 2
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (state$momenta$scale_beta_ancillary ^2) / 2
    proposed_U = 
      (
        +0.5*t(q - hierarchical_model$beta_priors$scale_beta$mean_whitened) %*% hierarchical_model$beta_priors$scale_beta$precision_whitened %*% (q - hierarchical_model$beta_priors$scale_beta$mean_whitened) # prior 
        +.5 * sum((data$observed_field - state$sparse_chol_and_stuff$lm_fit - new_field[vecchia_approx$locs_match] )^2/state$sparse_chol_and_stuff$noise)
      )
    proposed_K = sum(p^2) / 2
    if(!is.nan(current_U-proposed_U+current_K- proposed_K))
    {
      if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
      {
        state$momenta$scale_beta_ancillary = p
        state$params$scale_beta[] = new_scale_beta
        state$params$field = new_field
        acceptance_records$scale_beta_ancillary_mala[iter - 50*(iter %/% 50) ] = acceptance_records$scale_beta_ancillary_mala[iter - 50*(iter %/% 50) ] + 1
        state$sparse_chol_and_stuff$scale = new_scale
      }
    }
    # updating MALA kernel
    if(iter_start + iter < 1000)
    {
      if(iter %/% 50 ==iter / 50)
      {
        if(mean(acceptance_records$scale_beta_ancillary_mala)>.91)state$transition_kernels$scale_beta_ancillary_mala = state$transition_kernels$scale_beta_ancillary_mala +rnorm(1, .4, .05)
        if(mean(acceptance_records$scale_beta_ancillary_mala)<.61)state$transition_kernels$scale_beta_ancillary_mala = state$transition_kernels$scale_beta_ancillary_mala -rnorm(1, .4, .05)
        acceptance_records$scale_beta_ancillary_mala =  0*acceptance_records$scale_beta_ancillary_mala
      }
    }
    if((iter_start + iter > 1000)&(iter_start + iter < 2000))
    {
      if(iter %/% 50 ==iter / 50)
      {
        if(mean(acceptance_records$scale_beta_ancillary_mala)<.61)state$transition_kernels$scale_beta_ancillary_mala = state$transition_kernels$scale_beta_ancillary_mala -rnorm(1, .2, .05)
        acceptance_records$scale_beta_ancillary_mala =  0*acceptance_records$scale_beta_ancillary_mala
      }
    }
    
      #####################################
      # sufficient - ancillary analytical #
      #####################################
    if(!is.null(hierarchical_model$hyperprior_covariance$scale_NNGP_prior))
    {
      new_beta_and_field = Bidart::update_centered_beta(
        beta = state$params$scale_beta, 
        field = state$params$scale_field,
        X = data$covariates$scale_X,
        beta_prior = hierarchical_model$beta_priors$scale_beta, 
        NNGP_prior = hierarchical_model$hyperprior_covariance$scale_NNGP_prior, 
        log_scale = state$params$scale_log_scale
      )
      state$params$scale_beta[] = new_beta_and_field$beta
      state$params$scale_field[] = new_beta_and_field$field
    }
    
    #############################################
    # Nonstationary with log-scale latent field #
    #############################################
    if(!is.null(hierarchical_model$hyperprior_covariance$scale_NNGP_prior))
    {
      
      # recomputation to avoid errors
      state$sparse_chol_and_stuff$scale = Bidart::variance_field(beta = state$params$scale_beta, X = data$covariates$scale_X$X_locs, field = state$params$scale_field)

      ######################
      # scale latent field #
      ######################
      if(iter+iter_start >50) # in some cases waiting for the fixed effects to adjust avoids HMC crashing
      {
        sparse_chol_diag_field = state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = state$params$field)
        ###############################################
        # MALA update with sufficient parametrization #
        ###############################################
        q = exp(-.5*state$params$scale_log_scale)*as.vector(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol %*% state$params$scale_field)
        current_U =
          (
            .5 * sum(q^2) # whitenend prior 
            +0.5*sum(log(state$sparse_chol_and_stuff$scale))# determinant part
            +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/state$sparse_chol_and_stuff$scale))^2)# covmat product part
          )
        # MALA whitened
        state$momenta$scale_field_sufficient = sqrt(.9) * state$momenta$scale_field_sufficient + sqrt(.1)*rnorm(vecchia_approx$n_locs)
        p = state$momenta$scale_field_sufficient
        # Make a. half step for momentum at the beginning
        p = p - exp(state$transition_kernels$scale_field_sufficient_mala) *
          (q  # white noise prior derivative
           + exp(.5*state$params$scale_log_scale)*as.vector(
             Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                           .5  # determinant part 
                           -.5 * sqrt(1/state$sparse_chol_and_stuff$scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/state$sparse_chol_and_stuff$scale)))# natural derivative
             )))/ 2
        
        for(i in seq(5))
        {
          # Make a full step for the position
          q = q + exp(state$transition_kernels$scale_field_sufficient_mala) * p
          new_scale_field = exp(.5*state$params$scale_log_scale)*as.vector(Matrix::solve(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol, q)) 
          new_scale = state$sparse_chol_and_stuff$scale * exp(new_scale_field - state$params$scale_field)
          if(i != 5)
          {
            # Make a full step for momentum
            p = p - exp(state$transition_kernels$scale_field_sufficient_mala) *
              (q  # white noise prior derivative
               + exp(.5*state$params$scale_log_scale)*as.vector(
                 Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                               .5 -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
                 )))
          }
        }
        # Make a half step for momentum at the end.
        p = p - exp(state$transition_kernels$scale_field_sufficient_mala) *
          (q  # white noise prior derivative
           + exp(.5*state$params$scale_log_scale)*as.vector(
             Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                           .5 -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
             )))/ 2
        # Evaluate potential and kinetic energies at start and end of trajectory
        current_K = sum (state$momenta$scale_field_sufficient ^2) / 2
        proposed_U = 
          (
            .5 * sum(q^2) 
            +0.5*sum(log(new_scale))# determinant part
            +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/new_scale))^2)# covmat product part
          )
        proposed_K = sum(p^2) / 2
        
        if(!is.nan(current_U-proposed_U+current_K- proposed_K))
        {
          if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
          {
            state$momenta$scale_field_sufficient = p
            state$params$scale_field = new_scale_field
            acceptance_records$scale_field_sufficient_mala[iter - 50*(iter %/% 50) ] = acceptance_records$scale_field_sufficient_mala[iter - 50*(iter %/% 50) ] + 1
            state$sparse_chol_and_stuff$scale = new_scale
          }
        }
        # updating MALA kernel
        if(iter_start + iter < 1000)
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$scale_field_sufficient_mala)>.91)state$transition_kernels$scale_field_sufficient_mala = state$transition_kernels$scale_field_sufficient_mala +rnorm(1, .4, .05)
            if(mean(acceptance_records$scale_field_sufficient_mala)<.61)state$transition_kernels$scale_field_sufficient_mala = state$transition_kernels$scale_field_sufficient_mala -rnorm(1, .4, .05)
            acceptance_records$scale_field_sufficient_mala =  0*acceptance_records$scale_field_sufficient_mala
          }
        }
        if((iter_start + iter > 1000)&(iter_start + iter < 2000))        
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$scale_field_sufficient_mala)<.61)state$transition_kernels$scale_field_sufficient_mala = state$transition_kernels$scale_field_sufficient_mala -rnorm(1, .2, .05)
            acceptance_records$scale_field_sufficient_mala =  0*acceptance_records$scale_field_sufficient_mala
          }
        }
        ##############################################
        # MALA update with ancillary parametrization #
        ##############################################
        q = exp(-.5*state$params$scale_log_scale)*as.vector(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol %*% state$params$scale_field)
        current_U =
          (
            .5 * sum(q^2) # whitenend prior 
            +.5 * sum((data$observed_field - state$sparse_chol_and_stuff$lm_fit - state$params$field[vecchia_approx$locs_match] )^2/state$sparse_chol_and_stuff$noise)
          )
        # MALA whitened
        state$momenta$scale_field_ancillary = sqrt(.9) * state$momenta$scale_field_ancillary + sqrt(.1)*rnorm(vecchia_approx$n_locs)
        p = state$momenta$scale_field_ancillary
        # Make a. half step for momentum at the beginning
        p = p - exp(state$transition_kernels$scale_field_ancillary_mala) *
          (q  # white noise prior derivative
           + exp(.5*state$params$scale_log_scale)*as.vector(Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol), # solving by prior matrix
                                                                          .5 * state$params$field * 
                                                                            as.vector(vecchia_approx$locs_match_matrix %*% 
                                                                                        ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                                                                           state$sparse_chol_and_stuff$noise)
                                                                            )))
          )/ 2
        for(i in seq(5))
        {
          # Make a full step for the position
          q = q + exp(state$transition_kernels$scale_field_ancillary_mala) * p
          new_scale_field = exp(.5*state$params$scale_log_scale)*as.vector(Matrix::solve(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol, q)) 
          new_scale = state$sparse_chol_and_stuff$scale * exp(new_scale_field - state$params$scale_field)
          new_field  = state$params$field*sqrt(new_scale)/sqrt(state$sparse_chol_and_stuff$scale)
          # Make a full step for momentum.
          if(i !=5)
          {
            p = p - exp(state$transition_kernels$scale_field_ancillary_mala) *
              (q  # white noise prior derivative
               + exp(.5*state$params$scale_log_scale)*as.vector(Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol),
                                                                              .5 * new_field * 
                                                                                as.vector(vecchia_approx$locs_match_matrix %*% 
                                                                                            ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                                                                               state$sparse_chol_and_stuff$noise)
                                                                                )))
              )
          }
        }
        # Make a half step for momentum at the end.
        p = p - exp(state$transition_kernels$scale_field_ancillary_mala) *
          (q  # white noise prior derivative
           + exp(.5*state$params$scale_log_scale)*as.vector(Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol),
                                                                          .5 * new_field * 
                                                                            as.vector(vecchia_approx$locs_match_matrix %*% 
                                                                                        ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                                                                           state$sparse_chol_and_stuff$noise)
                                                                            )))
          )/ 2
        # Evaluate potential and kinetic energies at start and end of trajectory
        current_K = sum (state$momenta$scale_field_ancillary ^2) / 2
        proposed_U = 
          (
            .5 * sum(q^2) # whitenend prior 
            +.5 * sum((data$observed_field - state$sparse_chol_and_stuff$lm_fit - new_field[vecchia_approx$locs_match] )^2/state$sparse_chol_and_stuff$noise)
          )
        proposed_K = sum(p^2) / 2
        
        if(!is.nan(current_U-proposed_U+current_K- proposed_K))
        {
          if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
          {
            state$momenta$scale_field_ancillary = p
            state$params$scale_field = new_scale_field
            state$params$field = new_field
            acceptance_records$scale_field_ancillary_mala[iter - 50*(iter %/% 50) ] = acceptance_records$scale_field_ancillary_mala[iter - 50*(iter %/% 50) ] + 1
            state$sparse_chol_and_stuff$scale = new_scale
          }
        }
        # updating MALA kernel
        if(iter_start + iter < 1000)
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$scale_field_ancillary_mala)>.91)state$transition_kernels$scale_field_ancillary_mala = state$transition_kernels$scale_field_ancillary_mala +rnorm(1, .4, .05)
            if(mean(acceptance_records$scale_field_ancillary_mala)<.61)state$transition_kernels$scale_field_ancillary_mala = state$transition_kernels$scale_field_ancillary_mala -rnorm(1, .4, .05)
            acceptance_records$scale_field_ancillary_mala =  0*acceptance_records$scale_field_ancillary_mala
          }
        }
        if((iter_start + iter > 1000)&(iter_start + iter < 2000))
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$scale_field_ancillary_mala)<.61)state$transition_kernels$scale_field_ancillary_mala = state$transition_kernels$scale_field_ancillary_mala -rnorm(1, .2, .05)
            acceptance_records$scale_field_ancillary_mala =  0*acceptance_records$scale_field_ancillary_mala
          }
        }
      }
      ####################################################
      # Variance of the nonstationary scale latent field #
      ####################################################
      # ancillary -- sufficient 
      # a change in hyperprior scale changes (rescales) the scale latent field, which is then compared with the latent field
      new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, exp(state$transition_kernels$scale_sufficient_log_scale)) 
      new_scale_field = state$params$scale_field * exp((new_scale_log_scale - state$params$scale_log_scale)/2)
      new_scale = state$sparse_chol_and_stuff$scale * exp(new_scale_field - state$params$scale_field)
      if(
        (
          +.5* sum(log(state$sparse_chol_and_stuff$scale)) 
          -.5* sum(log(new_scale)) # log determinant
          +.5*sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          -.5*sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(new_scale)))^2) # Gaussian density of the latent field
          > log(runif(1))
        )
        & (new_scale_log_scale > -8)
        & (new_scale_log_scale < 3)
      )
      {
        state$params$scale_log_scale = new_scale_log_scale
        state$params$scale_field = new_scale_field 
        state$sparse_chol_and_stuff$scale = new_scale
        acceptance_records$scale_sufficient_log_scale[iter - 50*(iter %/% 50)] = acceptance_records$scale_sufficient_log_scale[iter - 50*(iter %/% 50)] +1
      }
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_sufficient_log_scale)>.51)state$transition_kernels$scale_sufficient_log_scale = state$transition_kernels$scale_sufficient_log_scale +rnorm(1, .4, .05)
          if(mean(acceptance_records$scale_sufficient_log_scale)<.11)state$transition_kernels$scale_sufficient_log_scale = state$transition_kernels$scale_sufficient_log_scale -rnorm(1, .4, .05)
          acceptance_records$scale_sufficient_log_scale =  0*acceptance_records$scale_sufficient_log_scale
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))      
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_sufficient_log_scale)<.11)state$transition_kernels$scale_sufficient_log_scale = state$transition_kernels$scale_sufficient_log_scale -rnorm(1, .2, .05)
          acceptance_records$scale_sufficient_log_scale =  0*acceptance_records$scale_sufficient_log_scale
        }
      }
      
      
      # sufficient -- sufficient
      unscaled_thingy = sum((hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol %*% state$params$scale_field)^2)
      for(i in seq(4))
      {
        new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, .1)
        if(
          (
            +.5*vecchia_approx$n_locs*state$params$scale_log_scale  +.5* exp(-state$params$scale_log_scale)*unscaled_thingy
            -.5*vecchia_approx$n_locs*new_scale_log_scale  -.5* exp(-new_scale_log_scale)*unscaled_thingy
            > log(runif(1))
          )
          & (new_scale_log_scale > -8)
          & (new_scale_log_scale < 3)
        )
        {
          state$params$scale_log_scale = new_scale_log_scale
        }
      }
      # ancillary -- ancillary
      new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, exp(state$transition_kernels$scale_ancillary_log_scale))
      new_scale_field = state$params$scale_field * exp((new_scale_log_scale - state$params$scale_log_scale)/2)
      new_scale = state$sparse_chol_and_stuff$scale * exp(new_scale_field - state$params$scale_field)
      new_field = state$params$field * sqrt(new_scale)/sqrt(state$sparse_chol_and_stuff$scale)
      if(
        (
          -.5* sum((state$sparse_chol_and_stuff$lm_residuals -          new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) 
          +.5* sum((state$sparse_chol_and_stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) 
          > log(runif(1))
        )
        & (new_scale_log_scale > -8)
        & (new_scale_log_scale < 3)
      )
      {
        state$params$field = new_field
        state$params$scale_log_scale = new_scale_log_scale
        state$params$scale_field = new_scale_field 
        state$sparse_chol_and_stuff$scale = new_scale
        acceptance_records$scale_ancillary_log_scale[iter - 50*(iter %/% 50)] = acceptance_records$scale_ancillary_log_scale[iter - 50*(iter %/% 50)] +1
      }
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_ancillary_log_scale)>.51)state$transition_kernels$scale_ancillary_log_scale = state$transition_kernels$scale_ancillary_log_scale +rnorm(1, .4, .05)
          if(mean(acceptance_records$scale_ancillary_log_scale)<.11)state$transition_kernels$scale_ancillary_log_scale = state$transition_kernels$scale_ancillary_log_scale -rnorm(1, .4, .05)
          acceptance_records$scale_ancillary_log_scale =  0*acceptance_records$scale_ancillary_log_scale
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_ancillary_log_scale)<.11)state$transition_kernels$scale_ancillary_log_scale = state$transition_kernels$scale_ancillary_log_scale -rnorm(1, .2, .05)
          acceptance_records$scale_ancillary_log_scale =  0*acceptance_records$scale_ancillary_log_scale
        }
      }
      
      # sufficient -- sufficient
      unscaled_thingy = sum((hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol %*% state$params$scale_field)^2)
      for(i in seq(4))
      {
        new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, .1)
        if(
          (
            +.5*vecchia_approx$n_locs*state$params$scale_log_scale  +.5* exp(-state$params$scale_log_scale)*unscaled_thingy
            -.5*vecchia_approx$n_locs*new_scale_log_scale  -.5* exp(-new_scale_log_scale)*unscaled_thingy
            > log(runif(1))
          )
          & (new_scale_log_scale > -8)
          & (new_scale_log_scale < 3)
        )
        {
          state$params$scale_log_scale = new_scale_log_scale
        }
      }
    }
    ###########################
    # Regression coefficients #
    ###########################
    t1 = Sys.time()
    # sufficient parametrization of latent field
    beta_covmat = solve(as.matrix(t(data$covariates$X$X) %*% Matrix::Diagonal(x = 1/state$sparse_chol_and_stuff$noise) %*% data$covariates$X$X))
    if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat)))
    {
      if(all(eigen(beta_covmat)$d >0))
      {
        beta_mean = c((((data$observed_field-state$params$field[vecchia_approx$locs_match]) / state$sparse_chol_and_stuff$noise) %*% data$covariates$X$X) %*% beta_covmat)
        state$params$beta[]   = c(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      }}
    # interweaving centered sampling in case of location-wise data to improve beta sampling
    centered_field = as.vector(state$params$field + data$covariates$X$X_locs%*%matrix(state$params$beta[data$covariates$X$which_locs], ncol = 1))
    sparse_chol_X = as.matrix(state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = 1/sqrt(state$sparse_chol_and_stuff$scale)) %*% data$covariates$X$X_locs)
    beta_precision = crossprod(sparse_chol_X)
    beta_covmat = solve(beta_precision, tol = min(rcond(beta_precision),.Machine$double.eps))
    if(all(!is.infinite(beta_covmat) & !is.nan(beta_covmat)))
    {
      if(all(eigen(beta_covmat)$d >0))
      {
        beta_mean =  c(as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (centered_field/sqrt(state$sparse_chol_and_stuff$scale)))  %*% sparse_chol_X %*% beta_covmat)
        state$params$beta[data$covariates$X$which_locs]   = as.vector(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
        state$params$field = centered_field - as.vector(data$covariates$X$X_locs %*% matrix(state$params$beta[data$covariates$X$which_locs], ncol = 1))
      }}
    # updating stuff 
    state$sparse_chol_and_stuff$lm_fit       = as.vector(data$covariates$X$X%*%     matrix(state$params$beta, ncol = 1))
    state$sparse_chol_and_stuff$lm_fit_locs  = as.vector(data$covariates$X$X_locs%*%matrix(state$params$beta[data$covariates$X$which_locs], ncol = 1))
    state$sparse_chol_and_stuff$lm_residuals = as.vector(data$observed_field-              state$sparse_chol_and_stuff$lm_fit)
    #print(Sys.time()-t1)
    ################
    # Latent field #
    ################    
    locs_partition = vecchia_approx$locs_partition[,runif(1, 1, ncol(vecchia_approx$locs_partition))]
    additional_precision = vecchia_approx$locs_match_matrix %*% (1/state$sparse_chol_and_stuff$noise)
    additional_mean = vecchia_approx$locs_match_matrix %*% ((data$observed_field - state$sparse_chol_and_stuff$lm_fit)/state$sparse_chol_and_stuff$noise)
    chol_sub_list = lapply(
      seq(length(unique(locs_partition))), function(idx)
      {
        selected_locs = which(locs_partition==idx)
        state$sparse_chol_and_stuff$sparse_chol[,selected_locs] %*% Matrix::Diagonal(x = 1/sqrt(state$sparse_chol_and_stuff$scale[selected_locs]))
      }
    )
    chol_Q_AA_sub_list = 
      lapply(
        seq(length(unique(locs_partition))), function(idx)
        {
          selected_locs = which(locs_partition==idx)
          Q_AA_sub =  Matrix::crossprod(chol_sub_list[[idx]])
          Matrix::diag(Q_AA_sub) = Matrix::diag(Q_AA_sub) + additional_precision[selected_locs]
          expanded = Matrix::expand(Matrix::Cholesky(Q_AA_sub, perm = T))
        }
      )
    Sys.time() - t1
    
    for(idx in seq(length(unique(locs_partition))))
    {
      selected_locs = which(locs_partition==idx)
      cond_mean = - 
        Matrix::t(chol_Q_AA_sub_list[[idx]]$P) %*%
        Matrix::solve(Matrix::t(chol_Q_AA_sub_list[[idx]]$L),
                      Matrix::solve(chol_Q_AA_sub_list[[idx]]$L, # inverse of precision matrix...
                                    chol_Q_AA_sub_list[[idx]]$P %*%
                                      (as.vector(Matrix::crossprod(chol_sub_list[[idx]], state$sparse_chol_and_stuff$sparse_chol%*%((state$params$field/sqrt(state$sparse_chol_and_stuff$scale))*(locs_partition!=idx)))) # rest of the latent field
                                       - additional_mean[selected_locs]))) # Gaussian observations 
      state$params$field [selected_locs] = as.vector(as.vector(cond_mean) +  Matrix::t(chol_Q_AA_sub_list[[idx]]$P) %*% Matrix::solve(Matrix::t(chol_Q_AA_sub_list[[idx]]$L), rnorm(length(selected_locs))) )
    }
    
    #######################
    # Storing the samples #
    #######################
    if(thinning * iter == floor(thinning*iter))
    {
      for(name in names(state$params))
      {
        params_records[[name]][,,thinning*iter] = state$params[[name]]
      }
    }
  }
  return(list("state" = state, "params_records" = params_records))
}
