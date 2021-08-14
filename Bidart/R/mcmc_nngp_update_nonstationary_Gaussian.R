
#' @export
mcmc_nngp_update_Gaussian = function(data,
                                     hierarchical_model, vecchia_approx, # model architecture
                                     state, # model state
                                     n_iterations_update = 400, thinning = .1, iter_start = 0, seed = 1,# practical settings
                                     field_n_chromatic = 1,  # number of chromatic steps for the latent field at each iteration
                                     field_n_mala = 1  # number of mala steps for the latent field at each iteration
                                     )
{
  #################
  # Sanity checks #
  #################
  # iterations and thinning
  if((floor(n_iterations_update)!=n_iterations_update) |  n_iterations_update<1)stop("n_iterations_update must be a positive round number")
  if((thinning<0)|(thinning>1))stop("thinning is a proportion and must be between 0 and 1")
  if((floor(iter_start)!=iter_start)|  iter_start<0)stop("iter_start must be a positive round number")
  # field
  if((floor(field_n_chromatic)!=field_n_chromatic)|  field_n_chromatic<0)stop("field_n_chromatic must be a positive round number")
  if((floor(field_n_mala)!=field_n_mala)|  field_n_mala<0)stop("field_n_mala must be a positive round number")
  if((field_n_chromatic==0) & (field_n_mala==0)) stop("Either field_n_chromatic or field_n_mala must be different from 0")
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
    #########
    # Range #
    #########
      #####################################
      # Stationary range beta (ancillary) #
      #####################################
    # Metropolis step
    if(is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior)&(ncol(data$covariates$range_X$X_locs)==1))
    {
      new_range_beta_0  = state$params$range_beta + exp(.5*state$transition_kernels$range_beta_ancillary) * rnorm(length(state$params$range_beta))
      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                range_beta = new_range_beta_0, NNarray = vecchia_approx$NNarray, 
                                                                locs = data$locs, 
                                                                range_field = NULL, 
                                                                range_X = NULL
      )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
      if(
          log(runif(1))<
          -.5 * sum((new_field         [vecchia_approx$locs_match]- state$sparse_chol_and_stuff$lm_residuals)^2/state$sparse_chol_and_stuff$noise)
          +.5 * sum((state$params$field[vecchia_approx$locs_match]- state$sparse_chol_and_stuff$lm_residuals)^2/state$sparse_chol_and_stuff$noise)
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
    }
      ######################################
      # Stationary range beta (sufficient) #
      ######################################
    # Metropolis step
    if(is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior)&(ncol(data$covariates$range_X$X_locs)==1))
    {
      new_range_beta_0  = state$params$range_beta + exp(.5*state$transition_kernels$range_beta_sufficient) * rnorm(length(state$params$range_beta))
      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                                                range_beta = new_range_beta_0, NNarray = vecchia_approx$NNarray, 
                                                                locs = data$locs, 
                                                                range_field = NULL, 
                                                                range_X = NULL
      )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      
      if(log(runif(1))< 
          sum(log(new_compressed_sparse_chol_and_grad[[1]][,1])) - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
          -.5 * sum((new_sparse_chol                         %*% (state$params$field /sqrt(state$sparse_chol_and_stuff$scale)))^2)
          +.5 * sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field /sqrt(state$sparse_chol_and_stuff$scale)))^2)
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
          if(mean(acceptance_records$range_beta_sufficient)<.21)state$transition_kernels$range_beta_sufficient = state$transition_kernels$range_beta_sufficient - rnorm(1, .4, .05)
          acceptance_records$range_beta_sufficient =  0*acceptance_records$range_beta_sufficient
        }
      }
    }
      ###############
      ########################################
      # Nonstationary range beta (ancillary) #
      ########################################
    if((!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior)|(ncol(data$covariates$range_X$X_locs)>1)))
    {
      q = t(solve(data$covariates$range_X$chol_solve_crossprod_X_locs)) %*% state$params$range_beta # whitening wrt covariates of the range
      current_U =
        (
          0 # improper prior 
          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        )
      # MALA whitened
      state$momenta$range_beta_ancillary = sqrt(.9) * state$momenta$range_beta_ancillary + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
      p = state$momenta$range_beta_ancillary
      # Make a half step for momentum at the beginning
      p = p - exp(state$transition_kernels$range_beta_ancillary) *
        (0  # improper prior derivative
         + as.matrix(
           solve(solve(data$covariates$range_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                 t(data$covariates$range_X$X_locs) %*%   # Jacobian of range field wrt range_beta
                   (# natural gradient of obs likelihood wrt range field
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
                                                                range_X = data$covariates$range_X$X_locs)
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
      # Make a half step for momentum at the end.
      p = p - exp(state$transition_kernels$range_beta_ancillary) *
        (0  # improper prior derivative
         + as.matrix(
           solve(solve(data$covariates$range_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                 t(data$covariates$range_X$X_locs) %*%   # Jacobian of range field wrt range_beta
                   (# natural gradient of obs likelihood wrt range field
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
        (
          0 # improper prior 
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
    }
      #######################################
      # Nonstationary range beta (centered) #
      #######################################
      if(!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior))
      {
        # rotating field and beta in order to get 3 decorrelated components (each component has spatial correlation tho)
        rotated_field = state$params$range_field  %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
        rotated_beta  = matrix(state$params$range_beta, nrow = ncol(data$covariates$range_X$X_locs)) %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
        # center
        centered_field = as.matrix(rotated_field + data$covariates$range_X$X_locs %*% rotated_beta)
        sparse_chol_X = data$covariates$range_X$sparse_chol_X_locs
        beta_covmat = data$covariates$range_X$solve_crossprod_sparse_chol_X_locs
        # beta mean =  w^t RtRX (XtRtRX)^{-1}(((Xt))) 
        beta_mean =  t(
          t(as.matrix(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol %*% centered_field))  %*% sparse_chol_X #vbw^t RtRX 
          %*% beta_covmat # (XtRtRX)^{-1} 
        )
        # sampling new beta
        rotated_beta = beta_mean +  t(chol(beta_covmat)) %*% matrix(rnorm(length(beta_mean)), nrow = nrow(beta_mean))
        # de-centering and de-rotating
        rotated_field = centered_field - as.matrix(data$covariates$range_X$X_locs %*% rotated_beta)
        state$params$range_field = rotated_field %*% chol(Bidart::expmat(state$params$range_log_scale))
        state$params$range_beta[] = rotated_beta %*% chol(Bidart::expmat(state$params$range_log_scale))
      }
      #########################################
      # Nonstationary range beta (sufficient) #
      #########################################
      if((!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior)|(ncol(data$covariates$range_X$X_locs)>1)))
      {
        q = t(solve(data$covariates$range_X$chol_solve_crossprod_X_locs)) %*% state$params$range_beta # whitening wrt covariates of the range
        current_U =
          (
            0 # improper prior 
            + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
            - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
          )
        # MALA whitened
        state$momenta$range_beta_sufficient = sqrt(.9) * state$momenta$range_beta_sufficient + sqrt(.1)*matrix(rnorm(length(q)), nrow = nrow(q))
        p = state$momenta$range_beta_sufficient
        # Make a half step for momentum at the beginning
        p = p - exp(state$transition_kernels$range_beta_sufficient) *
          (0  # improper prior derivative
           + as.matrix(
             solve(solve(data$covariates$range_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
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
                                                                          range_X = data$covariates$range_X$X_locs)
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        # Make a half step for momentum at the end.
        p = p - exp(state$transition_kernels$range_beta_sufficient) *
          (0  # improper prior derivative
           + as.vector(
             solve(solve(data$covariates$range_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
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
        proposed_U =
          (
            0 # improper prior 
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
      }
      #######################################
      # Nonstationary range beta (centered) #
      #######################################
    if(!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior))
    {
      # rotating field and beta in order to get 3 decorrelated components (each component has spatial correlation tho)
      rotated_field = state$params$range_field  %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
      rotated_beta  = matrix(state$params$range_beta, nrow = ncol(data$covariates$range_X$X_locs)) %*% solve(chol(Bidart::expmat(state$params$range_log_scale)))
      # center
      centered_field = as.matrix(rotated_field + data$covariates$range_X$X_locs %*% rotated_beta)
      sparse_chol_X = data$covariates$range_X$sparse_chol_X_locs
      beta_covmat = data$covariates$range_X$solve_crossprod_sparse_chol_X_locs
      # beta mean =  w^t RtRX (XtRtRX)^{-1}(((Xt))) 
      beta_mean =  t(
        t(as.matrix(hierarchical_model$hyperprior_covariance$range_NNGP_prior$sparse_chol %*% centered_field))  %*% sparse_chol_X #vbw^t RtRX 
        %*% beta_covmat # (XtRtRX)^{-1} 
      )
      # sampling new beta
      rotated_beta = beta_mean +  t(chol(beta_covmat)) %*% matrix(rnorm(length(beta_mean)), nrow = nrow(beta_mean))
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
                                                                range_X = data$covariates$range_X$X_locs)
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
          if(mean(acceptance_records$range_field_ancillary_mala)>.81)state$transition_kernels$range_field_ancillary_mala = state$transition_kernels$range_field_ancillary_mala + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_field_ancillary_mala)<.51)state$transition_kernels$range_field_ancillary_mala = state$transition_kernels$range_field_ancillary_mala - rnorm(1, .4, .05)
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
                                                                range_X = data$covariates$range_X$X_locs)
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
          if(mean(acceptance_records$range_field_sufficient_mala)>.8)state$transition_kernels$range_field_sufficient_mala = state$transition_kernels$range_field_sufficient_mala + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_field_sufficient_mala)<.51)state$transition_kernels$range_field_sufficient_mala = state$transition_kernels$range_field_sufficient_mala - rnorm(1, .4, .05)
          acceptance_records$range_field_sufficient_mala =  0*acceptance_records$range_field_sufficient_mala
        }
      }
    }
      ####################################################
      # Variance of the nonstationary range latent field #
      ####################################################
    if(!is.null(hierarchical_model$hyperprior_covariance$range_NNGP_prior)){

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
      state$momenta$range_log_scale_sufficient = sqrt(.6) * state$momenta$range_log_scale_sufficient + sqrt(.4)*rnorm(length(new_range_log_scale))
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
                                                                  range_X = data$covariates$range_X$X_locs)
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
          if(mean(acceptance_records$range_log_scale_sufficient)>.81)state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_log_scale_sufficient)<.51)state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient - rnorm(1, .4, .05)
          acceptance_records$range_log_scale_sufficient =  0*acceptance_records$range_log_scale_sufficient
        }
      }
      if((iter_start + iter > 1000) & (iter_start + iter < 1500))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_log_scale_sufficient)<.51)state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient - rnorm(1, .4, .05)
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
                                                                          range_X = data$covariates$range_X$X_locs)
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
          if((mean(acceptance_records$range_log_scale_ancillary)>.81) )state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary + rnorm(1, .4, .05)
          if(mean(acceptance_records$range_log_scale_ancillary)<.51)state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary - rnorm(1, .4, .05)
          acceptance_records$range_log_scale_ancillary =  0*acceptance_records$range_log_scale_ancillary
        }
      }
      
      if((iter_start + iter > 1000)&(iter_start + iter < 1500))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$range_log_scale_ancillary)<.51)state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary - rnorm(1, .4, .05)
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
    # MH sampling within Gibbs sweep over all regression coefficients (which are only an intercept in the stationary case)
    whitened_noise_beta = t(solve(data$covariates$noise_X$chol_solve_crossprod_X)) %*% state$params$noise_beta
    for(i in seq(ncol(data$covariates$noise_X$X_white)))
    {
      innovation = rnorm(1, 0, exp(state$transition_kernels$noise_beta[i]))
      new_noise = state$sparse_chol_and_stuff$noise * exp(data$covariates$noise_X$X_white[,i])^innovation
      if(
        -.5* sum(log(new_noise))-.5*sum(squared_residuals/new_noise)
        +.5* sum(log(state$sparse_chol_and_stuff$noise))+.5*sum(squared_residuals/state$sparse_chol_and_stuff$noise)
        > log(runif(1))
      )
      {
        state$sparse_chol_and_stuff$noise = new_noise 
        whitened_noise_beta[i] = whitened_noise_beta[i] + innovation
        acceptance_records$noise_beta[iter - 50*(iter %/% 50), i] = acceptance_records$noise_beta[iter - 50*(iter %/% 50), i] + 1
      }
    }
    state$params$noise_beta[] = t(data$covariates$noise_X$chol_solve_crossprod_X) %*% whitened_noise_beta
    # updating MALA kernel
    if(iter_start + iter <1000)
    {
      if(iter %/% 50 ==iter / 50)
      {
        acceptance_means = apply(acceptance_records$noise_beta, 2, mean)
        state$transition_kernels$noise_beta[acceptance_means>.41] = state$transition_kernels$noise_beta[acceptance_means>.41] + rnorm(sum(acceptance_means>.41), .4, .05)
        state$transition_kernels$noise_beta[acceptance_means<.11] = state$transition_kernels$noise_beta[acceptance_means<.11] - rnorm(sum(acceptance_means<.11), .4, .05)
        acceptance_records$noise_beta =  0*acceptance_records$noise_beta
      }
    }
    if(iter_start + iter <1500 & iter_start + iter >1000)
    {
      if(iter %/% 50 ==iter / 50)
      {
        acceptance_means = apply(acceptance_records$noise_beta, 2, mean)
        state$transition_kernels$noise_beta[acceptance_means<.11] = state$transition_kernels$noise_beta[acceptance_means<.11] - rnorm(sum(acceptance_means<.11), .4, .05)
        acceptance_records$noise_beta =  0*acceptance_records$noise_beta
      }
    }


    if(!is.null(hierarchical_model$hyperprior_covariance$noise_NNGP_prior))
    {
      #######################################################
      # Gaussian update of centered regression coefficients #
      #######################################################
      centered_field = as.vector(state$params$noise_field + data$covariates$noise_X$X_locs %*% matrix(state$params$noise_beta[data$covariates$noise_X$which_locs], ncol = 1))
      sparse_chol_X = data$covariates$noise_X$sparse_chol_X_locs
      beta_covmat = data$covariates$noise_X$solve_crossprod_sparse_chol_X_locs
      beta_mean =  c(as.vector(hierarchical_model$hyperprior_covariance$noise_NNGP_prior$sparse_chol %*% centered_field)  %*% sparse_chol_X %*% beta_covmat)
      state$params$noise_beta [data$covariates$noise_X$which_locs] = as.vector(beta_mean + exp(.5 * state$params$noise_log_scale) * t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      state$params$noise_field = centered_field - as.vector(data$covariates$noise_X$X_locs %*% matrix(state$params$noise_beta[data$covariates$noise_X$which_locs], ncol = 1))
 
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
        if(iter_start + iter <1000)
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$noise_field_mala)>.81)state$transition_kernels$noise_field_mala= state$transition_kernels$noise_field_mala + rnorm(1, .4, .05)
            if(mean(acceptance_records$noise_field_mala)<.51)state$transition_kernels$noise_field_mala = state$transition_kernels$noise_field_mala - rnorm(1, .4, .05)
            acceptance_records$noise_field_mala =  0*acceptance_records$noise_field_mala
          }
        }
        if(iter_start + iter <1500 & iter_start + iter >1000)
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$noise_field_mala)<.51)state$transition_kernels$noise_field_mala = state$transition_kernels$noise_field_mala - rnorm(1, .4, .05)
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
      if(iter_start + iter <1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$noise_log_scale)>.41)state$transition_kernels$noise_log_scale = state$transition_kernels$noise_log_scale + rnorm(1, .4, .05)
          if(mean(acceptance_records$noise_log_scale)<.11)state$transition_kernels$noise_log_scale = state$transition_kernels$noise_log_scale - rnorm(1, .4, .05)
          acceptance_records$noise_log_scale =  0*acceptance_records$noise_log_scale
        }
      }
      if(iter_start + iter <1500 & iter_start + iter >1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$noise_log_scale)<.11)state$transition_kernels$noise_log_scale = state$transition_kernels$noise_log_scale - rnorm(1, .4, .05)
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
          )
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
      
    # sufficient HMC
   sparse_chol_diag_field = state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = state$params$field)
    q = t(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs)) %*% state$params$scale_beta
    current_U =
      (
        0 # improper prior 
        +0.5*sum(log(state$sparse_chol_and_stuff$scale))# determinant part
        +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/state$sparse_chol_and_stuff$scale))^2)# covmat product part
      )
    # HMC whitened
    state$momenta$scale_beta_sufficient = sqrt(.9) * state$momenta$scale_beta_sufficient + sqrt(.1)*rnorm(ncol(data$covariates$scale_X$X_locs))
    p = state$momenta$scale_beta_sufficient
    # Make a. half step for momentum at the beginning
    p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
      solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
            t(data$covariates$scale_X$X_locs) %*%
              (
                .5  # determinant part 
                -.5 * sqrt(1/state$sparse_chol_and_stuff$scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/state$sparse_chol_and_stuff$scale)))# natural derivative
                )
            )/ 2
    
    nsteps = rpois(1, 3)
    for(i in seq(nsteps))
    {
      # Make a full step for the position
      q = q + exp(state$transition_kernels$scale_beta_sufficient_mala) * p
      new_scale_beta = t(data$covariates$scale_X$chol_solve_crossprod_X_locs) %*% q
      new_scale = Bidart::variance_field(beta = new_scale_beta, X = data$covariates$scale_X$X_locs, field = state$params$scale_field)
      if(i != nsteps)
      {
        p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
          solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                t(data$covariates$scale_X$X_locs) %*%
                  (
                    .5  # determinant part 
                    -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
                  )
          )
      }
    }
    # Make a half step for momentum at the end.
    p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
      solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
            t(data$covariates$scale_X$X_locs) %*%
              (
                .5  # determinant part 
                -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
              )
      )/ 2
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (state$momenta$scale_beta_sufficient ^2) / 2
    proposed_U = 
      (
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
    if(iter_start + iter <1000)
    {
      if(iter %/% 50 ==iter / 50)
      {
        if(mean(acceptance_records$scale_beta_sufficient_mala)>.81)state$transition_kernels$scale_beta_sufficient_mala = state$transition_kernels$scale_beta_sufficient_mala +rnorm(1, .4, .05)
        if(mean(acceptance_records$scale_beta_sufficient_mala)<.51)state$transition_kernels$scale_beta_sufficient_mala = state$transition_kernels$scale_beta_sufficient_mala -rnorm(1, .4, .05)
        acceptance_records$scale_beta_sufficient_mala =  0*acceptance_records$scale_beta_sufficient_mala
      }
    }
    
    

    if(!is.null(hierarchical_model$hyperprior_covariance$scale_NNGP_prior))
    {
      centered_field = as.vector(state$params$scale_field + data$covariates$scale_X$X_locs %*% matrix(state$params$scale_beta, ncol = 1))
      sparse_chol_X = data$covariates$scale_X$sparse_chol_X_locs
      beta_covmat = data$covariates$scale_X$solve_crossprod_sparse_chol_X_locs
      beta_mean =  c(as.vector(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol %*% centered_field)  %*% sparse_chol_X %*% beta_covmat)
      state$params$scale_beta [data$covariates$scale_X$which_locs] = as.vector(beta_mean + exp(.5 * state$params$scale_log_scale) * t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      state$params$scale_field = centered_field - as.vector(data$covariates$scale_X$X_locs %*% matrix(state$params$scale_beta, ncol = 1))
    }
    
    # ancillary 
    q = t(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs)) %*% state$params$scale_beta
    current_U =
      (
        0 # improper prior 
        +.5 * sum((data$observed_field - state$sparse_chol_and_stuff$lm_fit - state$params$field[vecchia_approx$locs_match] )^2/state$sparse_chol_and_stuff$noise)
      )
    # MALA whitened
    state$momenta$scale_beta_ancillary = sqrt(.9) * state$momenta$scale_beta_ancillary   + sqrt(.1)*rnorm(ncol(data$covariates$scale_X$X_locs))
    p = state$momenta$scale_beta_ancillary
    # Make a. half step for momentum at the beginning
    p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
      solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
            t(data$covariates$scale_X$X_locs) %*%
              (.5 * state$params$field * 
                 as.vector(vecchia_approx$locs_match_matrix %*% 
                             ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                state$sparse_chol_and_stuff$noise)
                           )))/ 2
    nsteps = rpois(1, 3)
    for(i in seq(nsteps))
    {
      # Make a full step for the position
      q = q + exp(state$transition_kernels$scale_beta_ancillary_mala) * p
      new_scale_beta = t(data$covariates$scale_X$chol_solve_crossprod_X_locs) %*% q
      new_scale = Bidart::variance_field(beta = new_scale_beta, X = data$covariates$scale_X$X_locs, field = state$params$scale_field)
      new_field  = state$params$field*sqrt(new_scale)/sqrt(state$sparse_chol_and_stuff$scale)
      # Make a full step for momentum.
      if(i !=nsteps)
      {
        p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
          solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
                t(data$covariates$scale_X$X_locs) %*%
                  (.5 * new_field * 
                     as.vector(vecchia_approx$locs_match_matrix %*% 
                                 ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                    state$sparse_chol_and_stuff$noise)
                     )))
      }
    }
    # Make a half step for momentum at the end.
    p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
      solve(solve(data$covariates$scale_X$chol_solve_crossprod_X_locs), # solving by prior sparse chol because of whitening
            t(data$covariates$scale_X$X_locs) %*%
              (.5 * new_field * 
                 as.vector(vecchia_approx$locs_match_matrix %*% 
                             ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                state$sparse_chol_and_stuff$noise)
                 )))/ 2
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (state$momenta$scale_beta_ancillary ^2) / 2
    proposed_U = 
      (
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
    if(iter_start + iter <1000)
    {
      if(iter %/% 50 ==iter / 50)
      {
        if(mean(acceptance_records$scale_beta_ancillary_mala)>.81)state$transition_kernels$scale_beta_ancillary_mala = state$transition_kernels$scale_beta_ancillary_mala +rnorm(1, .4, .05)
        if(mean(acceptance_records$scale_beta_ancillary_mala)<.51)state$transition_kernels$scale_beta_ancillary_mala = state$transition_kernels$scale_beta_ancillary_mala -rnorm(1, .4, .05)
        acceptance_records$scale_beta_ancillary_mala =  0*acceptance_records$scale_beta_ancillary_mala
      }
    }
    
    if(!is.null(hierarchical_model$hyperprior_covariance$scale_NNGP_prior))
    {
      centered_field = as.vector(state$params$scale_field + data$covariates$scale_X$X_locs %*% matrix(state$params$scale_beta, ncol = 1))
      sparse_chol_X = data$covariates$scale_X$sparse_chol_X_locs
      beta_covmat = data$covariates$scale_X$solve_crossprod_sparse_chol_X_locs
      beta_mean =  c(as.vector(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol %*% centered_field)  %*% sparse_chol_X %*% beta_covmat)
      state$params$scale_beta [data$covariates$scale_X$which_locs] = as.vector(beta_mean + exp(.5 * state$params$scale_log_scale) * t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
      state$params$scale_field = centered_field - as.vector(data$covariates$scale_X$X_locs %*% matrix(state$params$scale_beta, ncol = 1))
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
        state$momenta$scale_field_sufficient = sqrt(.95) * state$momenta$scale_field_sufficient + sqrt(.05)*rnorm(vecchia_approx$n_locs)
        p = state$momenta$scale_field_sufficient
        # Make a. half step for momentum at the beginning
        p = p - exp(state$transition_kernels$scale_field_sufficient_mala) *
          (q  # white noise prior derivative
           + exp(.5*state$params$scale_log_scale)*as.vector(
             Matrix::solve(Matrix::t(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol), # solving by prior sparse chol because of whitening
                           .5  # determinant part 
                           -.5 * sqrt(1/state$sparse_chol_and_stuff$scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/state$sparse_chol_and_stuff$scale)))# natural derivative
             )))/ 2
        nsteps = rpois(1, 3)
        for(i in seq(nsteps))
        {
          # Make a full step for the position
          q = q + exp(state$transition_kernels$scale_field_sufficient_mala) * p
          new_scale_field = exp(.5*state$params$scale_log_scale)*as.vector(Matrix::solve(hierarchical_model$hyperprior_covariance$scale_NNGP_prior$sparse_chol, q)) 
          new_scale = state$sparse_chol_and_stuff$scale * exp(new_scale_field - state$params$scale_field)
          if(i != nsteps)
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
        if(iter_start + iter <1000)
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$scale_field_sufficient_mala)>.81)state$transition_kernels$scale_field_sufficient_mala = state$transition_kernels$scale_field_sufficient_mala +rnorm(1, .4, .05)
            if(mean(acceptance_records$scale_field_sufficient_mala)<.51)state$transition_kernels$scale_field_sufficient_mala = state$transition_kernels$scale_field_sufficient_mala -rnorm(1, .4, .05)
            acceptance_records$scale_field_sufficient_mala =  0*acceptance_records$scale_field_sufficient_mala
          }
        }
        if(iter_start + iter <1500 & iter_start + iter >1000)
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$scale_field_sufficient_mala)<.51)state$transition_kernels$scale_field_sufficient_mala = state$transition_kernels$scale_field_sufficient_mala -rnorm(1, .4, .05)
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
        state$momenta$scale_field_ancillary = sqrt(.95) * state$momenta$scale_field_ancillary + sqrt(.05)*rnorm(vecchia_approx$n_locs)
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
        if(iter_start + iter <1000)
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$scale_field_ancillary_mala)>.81)state$transition_kernels$scale_field_ancillary_mala = state$transition_kernels$scale_field_ancillary_mala +rnorm(1, .4, .05)
            if(mean(acceptance_records$scale_field_ancillary_mala)<.51)state$transition_kernels$scale_field_ancillary_mala = state$transition_kernels$scale_field_ancillary_mala -rnorm(1, .4, .05)
            acceptance_records$scale_field_ancillary_mala =  0*acceptance_records$scale_field_ancillary_mala
          }
        }
        if(iter_start + iter <1500 & iter_start + iter >1000)
        {
          if(iter %/% 50 ==iter / 50)
          {
            if(mean(acceptance_records$scale_field_ancillary_mala)<.51)state$transition_kernels$scale_field_ancillary_mala = state$transition_kernels$scale_field_ancillary_mala -rnorm(1, .4, .05)
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
      if(iter_start + iter <1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_sufficient_log_scale)>.51)state$transition_kernels$scale_sufficient_log_scale = state$transition_kernels$scale_sufficient_log_scale +rnorm(1, .4, .05)
          if(mean(acceptance_records$scale_sufficient_log_scale)<.11)state$transition_kernels$scale_sufficient_log_scale = state$transition_kernels$scale_sufficient_log_scale -rnorm(1, .4, .05)
          acceptance_records$scale_sufficient_log_scale =  0*acceptance_records$scale_sufficient_log_scale
        }
      }
      if(iter_start + iter <1500 & iter_start + iter >1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_sufficient_log_scale)<.11)state$transition_kernels$scale_sufficient_log_scale = state$transition_kernels$scale_sufficient_log_scale -rnorm(1, .4, .05)
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
      if(iter_start + iter <1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_ancillary_log_scale)>.51)state$transition_kernels$scale_ancillary_log_scale = state$transition_kernels$scale_ancillary_log_scale +rnorm(1, .4, .05)
          if(mean(acceptance_records$scale_ancillary_log_scale)<.11)state$transition_kernels$scale_ancillary_log_scale = state$transition_kernels$scale_ancillary_log_scale -rnorm(1, .4, .05)
          acceptance_records$scale_ancillary_log_scale =  0*acceptance_records$scale_ancillary_log_scale
        }
      }
      if(iter_start + iter <1500 & iter_start + iter >1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_ancillary_log_scale)<.11)state$transition_kernels$scale_ancillary_log_scale = state$transition_kernels$scale_ancillary_log_scale -rnorm(1, .4, .05)
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
    # sufficient parametrization of latent field
    beta_covmat = solve(as.matrix(t(data$covariates$X$X) %*% Matrix::Diagonal(x = 1/state$sparse_chol_and_stuff$noise) %*% data$covariates$X$X))
    beta_mean = c((((data$observed_field-state$params$field[vecchia_approx$locs_match]) / state$sparse_chol_and_stuff$noise) %*% data$covariates$X$X) %*% beta_covmat)
    state$params$beta[]   = c(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
    # interweaving centered sampling in case of location-wise data to improve beta sampling
    centered_field = as.vector(state$params$field + data$covariates$X$X_locs%*%matrix(state$params$beta[data$covariates$X$which_locs], ncol = 1))
    sparse_chol_X = as.matrix(state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = 1/sqrt(state$sparse_chol_and_stuff$scale)) %*% data$covariates$X$X_locs)
    beta_precision = crossprod(sparse_chol_X)
    beta_covmat = solve(beta_precision, tol = min(rcond(beta_precision),.Machine$double.eps))
    beta_mean =  c(as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (centered_field/sqrt(state$sparse_chol_and_stuff$scale)))  %*% sparse_chol_X %*% beta_covmat)
    state$params$beta[data$covariates$X$which_locs]   = as.vector(beta_mean + t(chol(beta_covmat)) %*% rnorm(length(beta_mean)))
    state$params$field = centered_field - as.vector(data$covariates$X$X_locs %*% matrix(state$params$beta[data$covariates$X$which_locs], ncol = 1))
    # updating stuff 
    state$sparse_chol_and_stuff$lm_fit       = as.vector(data$covariates$X$X%*%     matrix(state$params$beta, ncol = 1))
    state$sparse_chol_and_stuff$lm_fit_locs  = as.vector(data$covariates$X$X_locs%*%matrix(state$params$beta[data$covariates$X$which_locs], ncol = 1))
    state$sparse_chol_and_stuff$lm_residuals = as.vector(data$observed_field-              state$sparse_chol_and_stuff$lm_fit)
    ################
    # Latent field #
    ################
      #####################
      # Chromatic sampler #
      #####################
    # spatial-location-wise sum of Linear regression residuals, weighted by noise precision
    weighted_residuals_sum = as.vector(vecchia_approx$locs_match_matrix%*%(as.vector(data$observed_field-  data$covariates$X$X%*%matrix(state$params$beta, ncol = 1))/state$sparse_chol_and_stuff$noise)) 
    # Diagonal of posterior precision, obtained by summing the diagonal of the prior and some extra precision brought by Gaussian observations
    posterior_precision_diag = as.vector(
      state$sparse_chol_and_stuff$precision_diag/state$sparse_chol_and_stuff$scale # diag from prior precision
      + vecchia_approx$locs_match_matrix %*%(1/state$sparse_chol_and_stuff$noise) # diag from observations precision
    )
    
    if(field_n_chromatic!=0)
    {
      for(i in seq(field_n_chromatic)){
        for(color_idx in unique(vecchia_approx$coloring))
        {
          selected_locs = which(vecchia_approx$coloring==color_idx)
          #conditional mean
          cond_mean = # mu A/B = inv(Q_AA)Q_AB (X_B - mu_B)
            - (1/posterior_precision_diag[selected_locs]) * # inverse of Q_AA
            ( 
              as.vector(Matrix::crossprod(
                state$sparse_chol_and_stuff$sparse_chol[,selected_locs] %*% Matrix::Diagonal(x = 1/sqrt(state$sparse_chol_and_stuff$scale[selected_locs])), 
                state$sparse_chol_and_stuff$sparse_chol %*%(state$params$field * (vecchia_approx$coloring!=color_idx)/sqrt(state$sparse_chol_and_stuff$scale)))) # rest of the latent field
              - weighted_residuals_sum[selected_locs]# Gaussian observations 
            )
          # field sampling
          state$params$field [selected_locs] = as.vector(cond_mean +rnorm(length(selected_locs))/sqrt(posterior_precision_diag[selected_locs]))
        }
      }
    }
      ########
      # MALA #
      ########
    if(field_n_mala!=0)
    {
      for(i in seq(field_n_mala)){
        state$momenta$field = sqrt(.9) * state$momenta$field + sqrt(.1)*rnorm(vecchia_approx$n_locs)
        p = state$momenta$field
        # Whitened MALA
        q = as.vector(state$sparse_chol_and_stuff$sparse_chol  %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))
        current_U =
          .5 * sum(q^2) + # whitenend prior 
          .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        # Make a. half step for momentum at the beginning
        p = p - exp(state$transition_kernels$latent_field_mala) *
          (q  # white noise prior derivative
           + as.vector( Matrix::solve(Matrix::t(state$sparse_chol_and_stuff$sparse_chol), # solving prior for whitening
                                      Matrix::Diagonal(x = sqrt(state$sparse_chol_and_stuff$scale)) %*%
                                        (
                                          as.vector(vecchia_approx$locs_match_matrix %*% 
                                                      ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                                         state$sparse_chol_and_stuff$noise)
                                          )
                                        )
           )
           )
          )/2
        # Make a full step for the position
        q = q + exp(state$transition_kernels$latent_field_mala) * p
        new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(state$sparse_chol_and_stuff$sparse_chol, q)) 
        # Make a half step for momentum at the end.
        p = p - exp(state$transition_kernels$latent_field_mala) *
          (q  # white noise prior derivative
           + as.vector( Matrix::solve(Matrix::t(state$sparse_chol_and_stuff$sparse_chol), # solving prior for whitening
                                      Matrix::Diagonal(x = sqrt(state$sparse_chol_and_stuff$scale)) %*%
                                        (
                                          as.vector(vecchia_approx$locs_match_matrix %*% 
                                                      ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                                         state$sparse_chol_and_stuff$noise)
                                          )
                                        )
           )
           )
          )/2
        p = - p
        # Evaluate potential and kinetic energies at start and end of trajectory
        current_K = sum (state$momenta$field ^2) / 2
        proposed_U = 
          .5 * sum(q^2) +
          .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise)
        proposed_K = sum(p^2) / 2
        if(!is.nan(current_U-proposed_U+current_K- proposed_K))
        {
          if (log(runif(1)) < current_U-proposed_U+current_K- proposed_K)
          {
            current_U = proposed_U
            state$momenta$field = p
            state$params$field = new_field
            acceptance_records$latent_field_mala[iter - 50*(iter %/% 50) ] = acceptance_records$latent_field_mala[iter - 50*(iter %/% 50) ] + 1
          }
        }
      }
      # updating MALA kernel
      if(iter_start + iter <1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$latent_field_mala)>.81)state$transition_kernels$latent_field_mala = state$transition_kernels$latent_field_mala + rnorm(1, .4, .05)
          if(mean(acceptance_records$latent_field_mala)<.51)state$transition_kernels$latent_field_mala = state$transition_kernels$latent_field_mala - rnorm(1, .4, .05)
          acceptance_records$latent_field_mala =  0*acceptance_records$latent_field_mala
        }
      }
      if(iter_start + iter <1500 & iter_start + iter >1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$latent_field_mala / field_n_mala)<.51)state$transition_kernels$latent_field_mala = state$transition_kernels$latent_field_mala - rnorm(1, .4, .05)
          acceptance_records$latent_field_mala =  0*acceptance_records$latent_field_mala
        }
      }
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
