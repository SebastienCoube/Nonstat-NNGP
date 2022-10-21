
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
    print(iter)
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
      new_compressed_sparse_chol_and_grad = 
        Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                    range_beta = new_range_beta_0, NNarray = vecchia_approx$NNarray, 
                                    locs = data$locs, 
                                    range_X = NULL, 
                                    compute_derivative = F, 
                                    nu = hierarchical_model$nu, 
                                    KL = NULL, use_KL = F, locs_idx = NULL
        )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
      if(
        log(runif(1))<
        -.5 * sum((new_field         [vecchia_approx$locs_match]- state$sparse_chol_and_stuff$lm_residuals)^2/state$sparse_chol_and_stuff$noise)
        +.5 * sum((state$params$field[vecchia_approx$locs_match]- state$sparse_chol_and_stuff$lm_residuals)^2/state$sparse_chol_and_stuff$noise)
        + Bidart::beta_prior_ll(beta = new_range_beta_0, n_KL = 0, 
                                beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                log_scale = NULL) # normal prior
        - Bidart::beta_prior_ll(beta = state$params$range_beta, n_KL = 0, 
                                beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                log_scale = NULL) # normal prior
        # prior
        
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
      new_compressed_sparse_chol_and_grad = 
        Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                    range_beta = new_range_beta_0, NNarray = vecchia_approx$NNarray, 
                                    locs = data$locs, 
                                    range_X = NULL, 
                                    compute_derivative = F, 
                                    nu = hierarchical_model$nu, 
                                    KL = NULL, use_KL = F, locs_idx = NULL
        )
      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
      
      if(log(runif(1))< 
         sum(log(new_compressed_sparse_chol_and_grad[[1]][,1])) - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
         -.5 * sum((new_sparse_chol                         %*% (state$params$field /sqrt(state$sparse_chol_and_stuff$scale)))^2)
         +.5 * sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field /sqrt(state$sparse_chol_and_stuff$scale)))^2)
         # prior
         + Bidart::beta_prior_ll(beta = new_range_beta_0, n_KL = 0, 
                                 beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                 beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                 log_scale = NULL) # normal prior
         - Bidart::beta_prior_ll(beta = state$params$range_beta, n_KL = 0, 
                                 beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                 beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                 log_scale = NULL) # normal prior
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

    ########################################
    # Nonstationary range beta (ancillary) #
    ########################################
    if(length(grep("nonstat", hierarchical_model$covfun))==1)
    {
      q = data$covariates$range$chol_crossprod_X_locs %*% state$params$range_beta # whitening wrt covariates of the range
      current_U =
        (
          - Bidart::beta_prior_ll(beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                                  beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                  beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                  log_scale = state$params$range_log_scale) # normal prior
          # normal prior 
          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
        )
      # MALA whitened
      state$momenta$range_beta_ancillary = sqrt(.95) * state$momenta$range_beta_ancillary + sqrt(.05)*matrix(rnorm(length(q)), nrow = nrow(q))
      p = state$momenta$range_beta_ancillary
      # Make a half step for momentum at the beginning
      p = p - exp(state$transition_kernels$range_beta_ancillary) *
        as.matrix(
          solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                - Bidart::beta_prior_ll_derivative(beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                  beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                  beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                  log_scale = state$params$range_log_scale) # normal prior
                 # normal prior derivative                
                + Bidart::X_KL_crossprod(
                  X = data$covariates$range_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, locs_idx = vecchia_approx$hctam_scol_1,
                  Y = # Jacobian of range field wrt range_beta
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
                      )))))/ 2

      
#source("Bidart/R/Useful_stuff.R")
#state$params$range_beta[] = rnorm(length(state$params$range_beta[]))
#state$params$range_log_scale[] = 0
#d1 =         - Bidart::beta_prior_ll_derivative(beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
#                                           beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
#                                           beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
#                                           log_scale = state$params$range_log_scale) # normal prior
#beta_ = state$params$range_beta 
#beta_[2] = beta_[2] + 1/10000
#d2 = 10000*(
#  beta_prior_ll(beta = beta_, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
#                        beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
#                        beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
#                        log_scale = state$params$range_log_scale) # normal prior
#  -
#  beta_prior_ll(beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
#                        beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
#                        beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
#                        log_scale = state$params$range_log_scale) # normal prior
#)
#d1/d2


#source("Bidart/R/Useful_stuff.R")
## recomputing current sparse chol
#state$params$range_beta[] = rnorm(length(state$params$range_beta[]))
#state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = 
#  compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
#                      range_beta = state$params$range_beta, NNarray = vecchia_approx$NNarray, 
#                      locs = data$locs, 
#                      range_X = data$covariates$range_X$X_locs, 
#                      nu = hierarchical_model$nu, 
#                      KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL,
#                      compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
#  )
#state$sparse_chol_and_stuff$sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
#
## compute gradient using derivative of sparse chol
#d1 = X_KL_crossprod(
#  X = data$covariates$range_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, locs_idx = vecchia_approx$hctam_scol_1,
#  Y = # Jacobian of range field wrt range_beta
#    (
#      # natural gradient of obs likelihood wrt range field
#      derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
#                                    left_vector = as.vector(
#                                      Matrix::solve(
#                                        Matrix::t(state$sparse_chol_and_stuff$sparse_chol), 
#                                        - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
#                                                      ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
#                                        * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
#                                      )), 
#                                    right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
#                                    NNarray = vecchia_approx$NNarray  
#      )))[1,1]
## compute gradient using finite diff
#new_range_beta = state$params$range_beta
#new_range_beta[1, 1] = new_range_beta[1, 1] + .0001
#new_compressed_sparse_chol_and_grad = 
#  compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
#                              range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
#                              locs = data$locs, 
#                              range_X = data$covariates$range_X$X_locs, 
#                              nu = hierarchical_model$nu, 
#                              KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL,
#                              compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
#  )
#new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
#new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
#d2 = 10000 * (
#  + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
#  - 
#    + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
#)
#d1/d2
      
      nsteps = 1
      for(i in seq(nsteps))
      {
        # Make a full step for the position
        q = q + exp(state$transition_kernels$range_beta_ancillary) * p
        new_range_beta = solve(data$covariates$range_X$chol_crossprod_X_locs, q )
        new_compressed_sparse_chol_and_grad = 
          Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                      range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
                                      locs = data$locs, 
                                      range_X = data$covariates$range_X$X_locs, 
                                      nu = hierarchical_model$nu, 
                                      KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL,
                                      compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
          )
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
        # Make a full step for momentum
        if(i!=nsteps)
        {
          
          p = p - exp(state$transition_kernels$range_beta_ancillary) *
            as.matrix(
              solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                    - Bidart::beta_prior_ll_derivative
                    (beta = new_range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                      beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                      beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                      log_scale = state$params$range_log_scale) # normal prior
                    #normal prior derivative                
                    + Bidart::X_KL_crossprod(
                      X = data$covariates$range_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, locs_idx = vecchia_approx$hctam_scol_1,
                      Y = # Jacobian of range field wrt range_beta
                        (
                          # natural gradient of obs likelihood wrt range field
                          Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                        left_vector = as.vector(
                                                          Matrix::solve(
                                                            Matrix::t(new_sparse_chol), 
                                                            - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                          ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                                            * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                          )), 
                                                        right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                        NNarray = vecchia_approx$NNarray  
                          )))))
        }
      }
      # Make a half step for momentum at the end.
      p = p - exp(state$transition_kernels$range_beta_ancillary) *
        as.matrix(
          solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                - Bidart::beta_prior_ll_derivative
                (beta = new_range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                  beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                  beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                  log_scale = state$params$range_log_scale) # normal prior
                 #normal prior derivative                
                + Bidart::X_KL_crossprod(
                  X = data$covariates$range_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, locs_idx = vecchia_approx$hctam_scol_1,
                  Y = # Jacobian of range field wrt range_beta
                    (
                      # natural gradient of obs likelihood wrt range field
                      Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                    left_vector = as.vector(
                                                      Matrix::solve(
                                                        Matrix::t(new_sparse_chol), 
                                                        - as.vector(vecchia_approx$locs_match_matrix %*%  # gradient of  Gaussian observations ll wrt latent field
                                                                      ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals) / state$sparse_chol_and_stuff$noise))
                                                        * sqrt(state$sparse_chol_and_stuff$scale) # part of sparse chol
                                                      )), 
                                                    right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                    NNarray = vecchia_approx$NNarray  
                      )))))/ 2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$range_beta_ancillary ^2) / 2
      proposed_U =
        (
          - Bidart::beta_prior_ll(beta = new_range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                                  beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                  beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                  log_scale = state$params$range_log_scale) # normal prior
          # normal prior 
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
    #########################################
    # Nonstationary range beta (sufficient) #
    #########################################
    if(length(grep("nonstat", hierarchical_model$covfun))==1)
    {
      q = data$covariates$range$chol_crossprod_X_locs %*% state$params$range_beta # whitening wrt covariates of the range
      current_U =
        (
          - Bidart::beta_prior_ll(beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                                  beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                  beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                  log_scale = state$params$range_log_scale) # normal prior
          # normal prior 
          + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
        )
      
      # MALA whitened
      state$momenta$range_beta_sufficient = sqrt(.95) * state$momenta$range_beta_sufficient + sqrt(.05)*matrix(rnorm(length(q)), nrow = nrow(q))
      p = state$momenta$range_beta_sufficient
      # Make a half step for momentum at the beginning
      p = p - exp(state$transition_kernels$range_beta_sufficient) *
        as.matrix(
          solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                - Bidart::beta_prior_ll_derivative
                (beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                  beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                  beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                  log_scale = state$params$range_log_scale) # normal prior
                # normal prior derivative                
                + Bidart::X_KL_crossprod(
                  X = data$covariates$range_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, locs_idx = vecchia_approx$hctam_scol_1,
                  Y = # Jacobian of range field wrt range_beta
                    (# natural gradient of obs likelihood wrt range field
                      Bidart::derivative_sandwiches(derivatives = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                    left_vector = as.vector(state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                    right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                    NNarray = vecchia_approx$NNarray  
                      )
                      - Bidart::log_determinant_derivatives(sparse_chol_and_grad = state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
                    ))))/ 2
      
      nsteps = 1
      for(i in seq(nsteps))
      {
        # Make a full step for the position
        q = q + exp(state$transition_kernels$range_beta_sufficient) * p
        new_range_beta = solve(data$covariates$range_X$chol_crossprod_X_locs, q )
        new_compressed_sparse_chol_and_grad =
          Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                      range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
                                      locs = data$locs, 
                                      range_X = data$covariates$range_X$X_locs, 
                                      nu = hierarchical_model$nu, 
                                      KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL,
                                      compute_derivative = T, locs_idx = vecchia_approx$hctam_scol_1
          )
        new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
        # Make a full step for momentum
        if(i!= nsteps)
        {
          p = p - exp(state$transition_kernels$range_beta_sufficient) *
            as.matrix(
              solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                    - Bidart::beta_prior_ll_derivative
                    (beta = new_range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                      beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                      beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                      log_scale = state$params$range_log_scale) # normal prior
                    # normal prior derivative                
                    + Bidart::X_KL_crossprod(
                      X = data$covariates$range_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, locs_idx = vecchia_approx$hctam_scol_1,
                      Y = # Jacobian of range field wrt range_beta
                        (# natural gradient of obs likelihood wrt range field
                          Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                        left_vector = as.vector(new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                        right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                        NNarray = vecchia_approx$NNarray  
                          )
                          - Bidart::log_determinant_derivatives(sparse_chol_and_grad = new_compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
                        ))))
        }
      }
      # Make a half step for momentum at the end.
      p = p - exp(state$transition_kernels$range_beta_sufficient) *
        as.matrix(
          solve(t(data$covariates$range$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                - Bidart::beta_prior_ll_derivative
                (beta = new_range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                  beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                  beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                  log_scale = state$params$range_log_scale) # normal prior
                # normal prior derivative                
                + Bidart::X_KL_crossprod(
                  X = data$covariates$range_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, locs_idx = vecchia_approx$hctam_scol_1,
                  Y = # Jacobian of range field wrt range_beta
                    (# natural gradient of obs likelihood wrt range field
                      Bidart::derivative_sandwiches(derivatives = new_compressed_sparse_chol_and_grad[[2]], # derivative of the (unscaled) NNGP factor
                                                    left_vector = as.vector(new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))), # left vector = whitened latent field
                                                    right_vector = state$params$field/sqrt(state$sparse_chol_and_stuff$scale), # scaled latent field, the scaling actually belongs to the derivative since the derivative must be scaled
                                                    NNarray = vecchia_approx$NNarray  
                      )
                      - Bidart::log_determinant_derivatives(sparse_chol_and_grad = new_compressed_sparse_chol_and_grad, NNarray = vecchia_approx$NNarray)# derivative of determinant
                    ))))/ 2
      # Evaluate potential and kinetic energies at start and end of trajectory
      current_K = sum (state$momenta$range_beta_sufficient ^2) / 2
      proposed_U =
        (
           - Bidart::beta_prior_ll(beta = new_range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                                   beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                   beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                   log_scale = state$params$range_log_scale) # normal prior
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
    
    ####################################################
    # Variance of the nonstationary range latent field #
    ####################################################
    if(hierarchical_model$range_KL){
##      # ancillary-sufficient ####
##      
##      new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, exp(.5 * state$transition_kernels$range_log_scale_sufficient))
##      old_eigen = eigen(Bidart::expmat(state$params$range_log_scale))$values
##      new_eigen = eigen(Bidart::expmat(new_range_log_scale))$values
##      new_range_beta = state$params$range_beta
##      new_range_beta[-seq(data$covariates$range_X$n_regressors),] = new_range_beta[-seq(data$covariates$range_X$n_regressors),] %*% 
##        solve(chol(Bidart::expmat(state$params$range_log_scale))) %*% chol(Bidart::expmat(new_range_log_scale))
##      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun,  
##                                                                        range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
##                                                                        locs = data$locs, range_X = data$covariates$range_X$X_locs, 
##                                                                        KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, 
##                                                                        nu = hierarchical_model$nu, locs_idx = vecchia_approx$hctam_scol_1, compute_derivative = T)
##      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
##      if(
##        (
##          + .5* sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
##          - sum(log(state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][,1]))
##          - .5* sum((new_sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
##          + sum(log(new_compressed_sparse_chol_and_grad[[1]][,1]))
##          > log(runif(1))
##        )
##      )
##      {
##        if(        
##          all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
##          all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2]))
##        )
##        {
##          state$params$range_log_scale = new_range_log_scale
##          state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
##          state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
##          state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
##        }
##        acceptance_records$range_log_scale_sufficient[iter - 50*(iter %/% 50) ] = acceptance_records$range_log_scale_sufficient[iter - 50*(iter %/% 50) ] + 1
##      }
##      # updating MALA kernel
##      if(iter_start + iter < 1000)
##      {
##        if(iter %/% 50 ==iter / 50)
##        {
##          if(mean(acceptance_records$range_log_scale_sufficient)>.41  & state$transition_kernels$range_log_scale_sufficient< (0))state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient + rnorm(1, .4, .05)
##          if(mean(acceptance_records$range_log_scale_sufficient)<.11)state$transition_kernels$range_log_scale_sufficient = state$transition_kernels$range_log_scale_sufficient - rnorm(1, .4, .05)
##          acceptance_records$range_log_scale_sufficient =  0*acceptance_records$range_log_scale_sufficient
##        }
##      }
      # sufficient - sufficient ####
      for(i in seq(5))
      {
        new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, .05)
        old_eigen = eigen(Bidart::expmat(state$params$range_log_scale))$values
        new_eigen = eigen(Bidart::expmat(new_range_log_scale))$values
        if(
          all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
          all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2]))&
          (
            + Bidart::beta_prior_ll(beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                                    beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                    beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                    log_scale = new_range_log_scale)
            - Bidart::beta_prior_ll(beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
                                    beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
                                    beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
                                    log_scale = state$params$range_log_scale) 
            > log(runif(1))
          )
          
        )
        {
          state$params$range_log_scale = new_range_log_scale
        }
      }
###      # ancillary - ancillary ####
###      new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, exp(.5 * state$transition_kernels$range_log_scale_ancillary))
###      old_eigen = eigen(Bidart::expmat(state$params$range_log_scale))$values
###      new_eigen = eigen(Bidart::expmat(new_range_log_scale))$values
###      new_range_beta = state$params$range_beta
###      new_range_beta[-seq(data$covariates$range_X$n_regressors),] = new_range_beta[-seq(data$covariates$range_X$n_regressors),] %*% 
###        solve(chol(Bidart::expmat(state$params$range_log_scale))) %*% chol(Bidart::expmat(new_range_log_scale))
###      new_compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun,  
###                                                                        range_beta = new_range_beta, NNarray = vecchia_approx$NNarray, 
###                                                                        locs = data$locs, range_X = data$covariates$range_X$X_locs, 
###                                                                        KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, 
###                                                                        nu = hierarchical_model$nu, locs_idx = vecchia_approx$hctam_scol_1, compute_derivative = T)
###      new_sparse_chol = Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = new_compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], triangular = T)
###      new_field = sqrt(state$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(new_sparse_chol, state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale))))
###      if(
###        ( 
###          - .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
###          + .5 * sum((state$sparse_chol_and_stuff$lm_residuals -  state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) # observation ll
###          > log(runif(1))
###        )
###      )
###      {
###        if(
###          all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
###          all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2]))
###        )
###        {
###          state$params$range_log_scale = new_range_log_scale
###          state$params$field = new_field
###          state$sparse_chol_and_stuff$sparse_chol= new_sparse_chol
###          state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = new_compressed_sparse_chol_and_grad
###          state$sparse_chol_and_stuff$precision_diag = as.vector((state$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
###        }
###        acceptance_records$range_log_scale_ancillary[iter - 50*(iter %/% 50) ] = acceptance_records$range_log_scale_ancillary[iter - 50*(iter %/% 50) ] + 1
###      }
###      # updating MALA kernel
###      if(iter_start + iter < 1000)
###      {
###        if(iter %/% 50 ==iter / 50)
###        {
###          if((mean(acceptance_records$range_log_scale_ancillary)>.41) & state$transition_kernels$range_log_scale_ancillary< (0))state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary + rnorm(1, .4, .05)
###          if(mean(acceptance_records$range_log_scale_ancillary)<.11)state$transition_kernels$range_log_scale_ancillary = state$transition_kernels$range_log_scale_ancillary - rnorm(1, .4, .05)
###          acceptance_records$range_log_scale_ancillary =  0*acceptance_records$range_log_scale_ancillary
###        }
###      }
###      
###      # sufficient - sufficient ######
###      for(i in seq(5))
###      {
###        new_range_log_scale = state$params$range_log_scale + rnorm(length(state$params$range_log_scale), 0, .05)
###        old_eigen = eigen(Bidart::expmat(state$params$range_log_scale))$values
###        new_eigen = eigen(Bidart::expmat(new_range_log_scale))$values
###        if(
###          all(min(new_eigen)>exp(hierarchical_model$range_log_scale_prior[1]))&
###          all(max(new_eigen)<exp(hierarchical_model$range_log_scale_prior[2]))&
###          (
###            + Bidart::beta_prior_ll(beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
###                                    beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
###                                    beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
###                                    log_scale = new_range_log_scale)
###            - Bidart::beta_prior_ll(beta = state$params$range_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$range_KL, 
###                                    beta_mean = hierarchical_model$beta_priors$range_beta_mean, 
###                                    beta_precision =  hierarchical_model$beta_priors$range_beta_precision, 
###                                    log_scale = state$params$range_log_scale) 
###            > log(runif(1))
###          )
###          
###        )
###        {
###          state$params$range_log_scale = new_range_log_scale
###        }
###      }
###      
    }
    #########
    # Noise #
    #########
      ##############
      # Noise beta #
      ##############
    # recomputation in order to avoid errors
    state$sparse_chol_and_stuff$noise = Bidart::variance_field(
      beta = state$params$noise_beta, X = data$covariates$noise_X$X, 
      KL = hierarchical_model$KL, use_KL = hierarchical_model$noise_KL)
    # VEWY IMPOWTANT don't remove or comment
    squared_residuals = as.vector(state$sparse_chol_and_stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2
    # HMC update
    q = data$covariates$noise_X$chol_crossprod_X %*% state$params$noise_beta
    current_U =
      (
        - Bidart::beta_prior_ll(beta = state$params$noise_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$noise_KL, 
                      beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                      beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
                      log_scale = state$params$noise_log_scale) # normal prior 
        +.5* sum(log(state$sparse_chol_and_stuff$noise)) # det
        +.5*sum(squared_residuals/state$sparse_chol_and_stuff$noise) # observations
      )
    # HMC whitened
    state$momenta$noise_beta = sqrt(.9) * state$momenta$noise_beta + sqrt(.1)*rnorm(length(state$momenta$noise_beta))
    p = state$momenta$noise_beta
    
    # Make a half step for momentum at the beginning
    p = p - exp(state$transition_kernels$noise_beta_mala) *
      (
        + solve(t(data$covariates$noise_X$chol_crossprod_X), # solving by prior chol because of whitening
              - Bidart::beta_prior_ll_derivative(beta = state$params$noise_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$noise_KL, 
                                       beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                                       beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
                                       log_scale = state$params$noise_log_scale) # normal prior
              + Bidart::X_KL_crossprod(X = data$covariates$noise_X$X, KL = hierarchical_model$KL, use_KL = hierarchical_model$noise_KL, 
                               Y = 
                                 (
                                    + .5 # determinant part of normal likelihood
                                    - (squared_residuals/state$sparse_chol_and_stuff$noise)/2 # exponential part of normal likelihood
                                 ))
      ))/ 2
    # checking gradient with finite differences, to update
    ###noise_beta_ = state$params$noise_beta 
    ###noise_beta_[1] = noise_beta_[1] + 0.0001
    ###noise_ =Bidart::variance_field(beta = noise_beta_, field = state$params$noise_field[vecchia_approx$locs_match], data$covariates$noise_X$X)
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
    nsteps = 8 + rbinom(1, 1, .5)
    for(i in seq(nsteps))
    {
      # Make a full step for the position
      q = q + exp(state$transition_kernels$noise_beta_mala) * p
      new_noise_beta = solve(data$covariates$noise_X$chol_crossprod_X, q)
      new_noise = Bidart::variance_field(beta = new_noise_beta, KL = hierarchical_model$KL, use_KL = hierarchical_model$noise_KL, 
                                 X = data$covariates$noise_X$X)
      # Make a full step for momentum at the end
      if(i!= nsteps)
      {
        p = p - exp(state$transition_kernels$noise_beta_mala) *
          (
            + solve(t(data$covariates$noise_X$chol_crossprod_X), # solving by prior sparse chol because of whitening
                  - Bidart::beta_prior_ll_derivative(beta = new_noise_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$noise_KL, 
                                           beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                                           beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
                                           log_scale = state$params$noise_log_scale) # normal prior
                  + Bidart::X_KL_crossprod(X = data$covariates$noise_X$X, KL = hierarchical_model$KL, use_KL = hierarchical_model$noise_KL, 
                                       (
                                         + .5 # determinant part of normal likelihood
                                         - (squared_residuals/new_noise)/2 # exponential part of normal likelihood
                                       ))
          ))
      }
    }
    # Make a half step for momentum at the end
    p = p - exp(state$transition_kernels$noise_beta_mala) *
      (
        + solve(t(data$covariates$noise_X$chol_crossprod_X), # solving by prior sparse chol because of whitening
              - Bidart::beta_prior_ll_derivative(beta = new_noise_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$noise_KL, 
                            beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                            beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
                            log_scale = state$params$noise_log_scale) # normal prior  
              +
                + Bidart::X_KL_crossprod(X = data$covariates$noise_X$X, KL = hierarchical_model$KL, use_KL = hierarchical_model$noise_KL, 
                                   (
                                     + .5 # determinant part of normal likelihood
                                     - (squared_residuals/new_noise)/2 # exponential part of normal likelihood
                                   ))
      ))/ 2
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (state$momenta$noise_beta ^2) / 2
    proposed_U = 
      (
        - Bidart::beta_prior_ll(beta = new_noise_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$noise_KL, 
                      beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                      beta_precision =  hierarchical_model$beta_priors$noise_beta_precision, 
                      log_scale = state$params$noise_log_scale) # normal prior        
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
    
      ###################
      # Noise log scale # 
      ###################
    if(hierarchical_model$noise_KL)
    {
      # ancillary -- sufficient ####
      new_noise_log_scale = state$params$noise_log_scale + rnorm(1, 0, exp(state$transition_kernels$noise_log_scale))
      new_noise_beta = state$params$noise_beta
      new_noise_beta[-seq(data$covariates$noise_X$n_regressors)] = new_noise_beta[-seq(data$covariates$noise_X$n_regressors)] *
        exp((new_noise_log_scale - state$params$noise_log_scale)/2)
      new_noise =Bidart::variance_field(beta = new_noise_beta, KL = hierarchical_model$KL, use_KL = hierarchical_model$noise_KL, X = data$covariates$noise_X$X)
      if(
        (
          -.5* sum(log(new_noise)) 
          -.5*sum(squared_residuals/new_noise)
          +.5* sum(log(state$sparse_chol_and_stuff$noise)) 
          +.5*sum(squared_residuals/state$sparse_chol_and_stuff$noise)
          > log(runif(1)
          )
        )
      )
      {
        if(
          (new_noise_log_scale > hierarchical_model$noise_log_scale_prior[1])&
          (new_noise_log_scale < hierarchical_model$noise_log_scale_prior[2])
        )
        {
          state$params$noise_log_scale = new_noise_log_scale
          state$params$noise_beta = new_noise_beta 
          state$sparse_chol_and_stuff$noise = new_noise
        }
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
      
      # sufficient -- sufficient ####
      for(i in seq(4))
      {
        new_noise_log_scale = state$params$noise_log_scale + rnorm(1, 0, .1)
        if(
          (
            + Bidart::beta_prior_ll(beta = state$params$noise_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$noise_KL, 
                          beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                          beta_precision = hierarchical_model$beta_priors$noise_beta_precision, 
                          log_scale = new_noise_log_scale) - 
            Bidart::beta_prior_ll(beta = state$params$noise_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$noise_KL, 
                          beta_mean = hierarchical_model$beta_priors$noise_beta_mean, 
                          beta_precision = hierarchical_model$beta_priors$noise_beta_precision, 
                          log_scale = state$params$noise_log_scale)
          ) > log(runif(1))
          &(new_noise_log_scale > hierarchical_model$noise_log_scale_prior[1])
          &(new_noise_log_scale < hierarchical_model$noise_log_scale_prior[2])
          )
        {
          state$params$noise_log_scale = new_noise_log_scale
        }
      }
    }
    
    
    #########
    # Scale #
    #########
    #########################################
    # scale beta ancillary - sufficient HMC #
    #########################################
    sparse_chol_diag_field = state$sparse_chol_and_stuff$sparse_chol %*% Matrix::Diagonal(x = state$params$field)
    q = data$covariates$scale_X$chol_crossprod_X_locs %*% state$params$scale_beta
    current_U =
      (
- Bidart::beta_prior_ll(beta = state$params$scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                log_scale = state$params$scale_log_scale) # normal prior 
        +0.5*sum(log(state$sparse_chol_and_stuff$scale))# determinant part
        +0.5*sum(as.vector(sparse_chol_diag_field%*%sqrt(1/state$sparse_chol_and_stuff$scale))^2)# covmat product part
      )
    # HMC whitened
    state$momenta$scale_beta_sufficient = sqrt(.9) * state$momenta$scale_beta_sufficient + sqrt(.1)*rnorm(length(state$momenta$scale_beta_sufficient))
    p = state$momenta$scale_beta_sufficient
    # Make a half step for momentum at the beginning
    p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
      (
        + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by t chol because of whitening
                - Bidart::beta_prior_ll_derivative(beta = state$params$scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                                           beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                           beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                           log_scale = state$params$scale_log_scale) # normal prior
                +  Bidart::X_KL_crossprod(X = data$covariates$scale_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, locs_idx = vecchia_approx$hctam_scol_1,
                  Y = (
                    .5  # determinant part 
                    -.5 * sqrt(1/state$sparse_chol_and_stuff$scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/state$sparse_chol_and_stuff$scale)))# natural derivative
                  ))
        ))/ 2
    n_hmc_steps = 8 + rbinom(1, 1, .5)
    for(i in seq(n_hmc_steps))
    {
      # Make a full step for the position
      q = q + exp(state$transition_kernels$scale_beta_sufficient_mala) * p
      new_scale_beta = solve(data$covariates$scale_X$chol_crossprod_X_locs, q)
      new_scale =Bidart::variance_field(beta = new_scale_beta, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, 
                                X = data$covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
      if(i != n_hmc_steps)
      {
        p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
          (
            + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                    - Bidart::beta_prior_ll_derivative(beta = new_scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                                               beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                               beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                               log_scale = state$params$scale_log_scale) # normal prior  
                    + Bidart::X_KL_crossprod(X = data$covariates$scale_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, locs_idx = vecchia_approx$hctam_scol_1,
                      (
                        .5  # determinant part 
                        -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
                      ))
            ))
      }
    }
    # Make a half step for momentum at the end.
    p = p - exp(state$transition_kernels$scale_beta_sufficient_mala) *
      (
        + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                - Bidart::beta_prior_ll_derivative(beta = new_scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                                           beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                           beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                           log_scale = state$params$scale_log_scale) # normal prior  
                + Bidart::X_KL_crossprod(X = data$covariates$scale_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, locs_idx = vecchia_approx$hctam_scol_1,
                  (
                    .5  # determinant part 
                    -.5 * sqrt(1/new_scale) * as.vector(Matrix::crossprod(sparse_chol_diag_field, sparse_chol_diag_field %*% sqrt(1/new_scale)))# natural derivative
                  ))
        ))/ 2
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (state$momenta$scale_beta_sufficient ^2) / 2
    proposed_U = 
      (
- Bidart::beta_prior_ll(beta = new_scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                log_scale = state$params$scale_log_scale) # normal prior 
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
    # updating HMC kernel
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
    
    ########################################
    # scale beta ancillary - ancillary HMC #
    ########################################
    q = data$covariates$scale_X$chol_crossprod_X_locs %*% state$params$scale_beta
    current_U =
      (
        - Bidart::beta_prior_ll(beta = state$params$scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                        beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                        beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                        log_scale = state$params$scale_log_scale) # normal prior 
        +.5 * sum((data$observed_field - state$sparse_chol_and_stuff$lm_fit - state$params$field[vecchia_approx$locs_match] )^2/state$sparse_chol_and_stuff$noise)
      )
    # MALA whitened
    state$momenta$scale_beta_ancillary = sqrt(.9) * state$momenta$scale_beta_ancillary   + sqrt(.1)*rnorm(length(state$momenta$scale_beta_ancillary))
    p = state$momenta$scale_beta_ancillary
    # Make a. half step for momentum at the beginning
    p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
      (
        + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior chol because of whitening
                - Bidart::beta_prior_ll_derivative(beta = state$params$scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                                           beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                           beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                           log_scale = state$params$scale_log_scale) # normal prior 
                + Bidart::X_KL_crossprod(X = data$covariates$scale_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, locs_idx = vecchia_approx$hctam_scol_1,
                  (.5 * state$params$field * 
                     as.vector(vecchia_approx$locs_match_matrix %*% 
                                 ((state$params$field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                    state$sparse_chol_and_stuff$noise)
                     )))
        ))/ 2
    n_hmc_steps = 8 + rbinom(1, 1, .5)
    for(i in seq(n_hmc_steps))
    {
      # Make a full step for the position
      q = q + exp(state$transition_kernels$scale_beta_ancillary_mala) * p
      new_scale_beta = solve(data$covariates$scale_X$chol_crossprod_X_locs, q)
      new_scale =Bidart::variance_field(beta = new_scale_beta, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, 
                                X = data$covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
      new_field  = state$params$field*sqrt(new_scale)/sqrt(state$sparse_chol_and_stuff$scale)
      # Make a full step for momentum.
      if(i != n_hmc_steps)
      {
        p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
          (
            + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                    - Bidart::beta_prior_ll_derivative(beta = new_scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                                               beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                               beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                               log_scale = state$params$scale_log_scale) # normal prior  
                    +
                      Bidart::X_KL_crossprod(X = data$covariates$scale_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, locs_idx = vecchia_approx$hctam_scol_1,
                      (.5 * new_field * 
                         as.vector(vecchia_approx$locs_match_matrix %*% 
                                     ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                        state$sparse_chol_and_stuff$noise)
                         ))))
          )
      }
    }
    # Make a half step for momentum at the end.
    p = p - exp(state$transition_kernels$scale_beta_ancillary_mala) *
      (
        + solve(t(data$covariates$scale_X$chol_crossprod_X_locs), # solving by prior sparse chol because of whitening
                - Bidart::beta_prior_ll_derivative(beta = new_scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                                           beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                           beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                           log_scale = state$params$scale_log_scale) # normal prior  
                + Bidart::X_KL_crossprod(X = data$covariates$scale_X$X_locs, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, locs_idx = vecchia_approx$hctam_scol_1,
                  (.5 * new_field * 
                     as.vector(vecchia_approx$locs_match_matrix %*% 
                                 ((new_field[vecchia_approx$locs_match] - state$sparse_chol_and_stuff$lm_residuals)/
                                    state$sparse_chol_and_stuff$noise)
                     )))
        ))/ 2
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_K = sum (state$momenta$scale_beta_ancillary ^2) / 2
    proposed_U = 
      (
        - Bidart::beta_prior_ll(beta = new_scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                        beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                        beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                        log_scale = state$params$scale_log_scale) # normal prior  
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
    
    
    ###################
    # scale log scale #
    ###################
    if(hierarchical_model$scale_KL){
      # ancillary -- sufficient  ####
      # a change in hyperprior scale changes (rescales) the scale, which is then compared with the latent field w
      new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, exp(state$transition_kernels$scale_log_scale_sufficient)) 
      new_scale_beta = state$params$scale_beta
      new_scale_beta[-seq(data$covariates$scale_X$n_regressors)] = new_scale_beta[-seq(data$covariates$scale_X$n_regressors)] *
        exp((new_scale_log_scale - state$params$scale_log_scale)/2)
      new_scale =Bidart::variance_field(beta = new_scale_beta, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, X = data$covariates$scale_X$X_locs)
      if(
        (
          +.5* sum(log(state$sparse_chol_and_stuff$scale)) 
          -.5* sum(log(new_scale)) # log determinant
          +.5*sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(state$sparse_chol_and_stuff$scale)))^2)
          -.5*sum((state$sparse_chol_and_stuff$sparse_chol %*% (state$params$field/sqrt(new_scale)))^2) # Gaussian density of the latent field
          > log(runif(1))
        )
      )
      {
        if(
          (new_scale_log_scale > hierarchical_model$scale_log_scale_prior[1])&
          (new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
        )
        {
          state$params$scale_log_scale = new_scale_log_scale
          state$params$scale_beta = new_scale_beta 
          state$sparse_chol_and_stuff$scale = new_scale
        }
        acceptance_records$scale_log_scale_sufficient[iter - 50*(iter %/% 50)] = acceptance_records$scale_log_scale_sufficient[iter - 50*(iter %/% 50)] +1
      }
      # updating kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_log_scale_sufficient)>.51)state$transition_kernels$scale_log_scale_sufficient = state$transition_kernels$scale_log_scale_sufficient +rnorm(1, .4, .05)
          if(mean(acceptance_records$scale_log_scale_sufficient)<.11)state$transition_kernels$scale_log_scale_sufficient = state$transition_kernels$scale_log_scale_sufficient -rnorm(1, .4, .05)
          acceptance_records$scale_log_scale_sufficient =  0*acceptance_records$scale_log_scale_sufficient
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))      
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_log_scale_sufficient)<.11)state$transition_kernels$scale_log_scale_sufficient = state$transition_kernels$scale_log_scale_sufficient -rnorm(1, .2, .05)
          acceptance_records$scale_log_scale_sufficient =  0*acceptance_records$scale_log_scale_sufficient
        }
      }
      
      
      # sufficient -- sufficient ####
      for(i in seq(4))
      {
        new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, .1)
        if(
          (
            + Bidart::beta_prior_ll(beta = state$params$scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                          beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                          beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                          log_scale = new_scale_log_scale)     
            - Bidart::beta_prior_ll(beta = state$params$scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                            beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                            beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                            log_scale = state$params$scale_log_scale)     
            > log(runif(1))
          )
          &(new_scale_log_scale > hierarchical_model$scale_log_scale_prior[1])
          &(new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
        )
        {
          state$params$scale_log_scale = new_scale_log_scale
        }
      }
      # ancillary -- ancillary ####
      new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, exp(state$transition_kernels$scale_log_scale_sufficient)) 
      new_scale_beta = state$params$scale_beta
      new_scale_beta[-seq(data$covariates$scale_X$n_regressors)] = new_scale_beta[-seq(data$covariates$scale_X$n_regressors)] *
        exp((new_scale_log_scale - state$params$scale_log_scale)/2)
      new_scale =Bidart::variance_field(new_scale_beta, KL = hierarchical_model$KL, use_KL = hierarchical_model$scale_KL, 
                                X = data$covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
      new_field = state$params$field * sqrt(new_scale)/sqrt(state$sparse_chol_and_stuff$scale)
      if(
        (
          -.5* sum((state$sparse_chol_and_stuff$lm_residuals -          new_field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) 
          +.5* sum((state$sparse_chol_and_stuff$lm_residuals - state$params$field[vecchia_approx$locs_match])^2/state$sparse_chol_and_stuff$noise) 
          > log(runif(1))
        )
      )
      {
        if(
          (new_scale_log_scale >  hierarchical_model$scale_log_scale_prior[1])
          &(new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
        )
        {  
          state$params$field = new_field
          state$params$scale_log_scale = new_scale_log_scale
          state$params$scale_beta = new_scale_beta 
          state$sparse_chol_and_stuff$scale = new_scale
        }
        acceptance_records$scale_log_scale_ancillary[iter - 50*(iter %/% 50)] = acceptance_records$scale_log_scale_ancillary[iter - 50*(iter %/% 50)] +1
      }
      # updating MALA kernel
      if(iter_start + iter < 1000)
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_log_scale_ancillary)>.51)state$transition_kernels$scale_log_scale_ancillary = state$transition_kernels$scale_log_scale_ancillary +rnorm(1, .4, .05)
          if(mean(acceptance_records$scale_log_scale_ancillary)<.11)state$transition_kernels$scale_log_scale_ancillary = state$transition_kernels$scale_log_scale_ancillary -rnorm(1, .4, .05)
          acceptance_records$scale_log_scale_ancillary =  0*acceptance_records$scale_log_scale_ancillary
        }
      }
      if((iter_start + iter > 1000)&(iter_start + iter < 2000))
      {
        if(iter %/% 50 ==iter / 50)
        {
          if(mean(acceptance_records$scale_log_scale_ancillary)<.11)state$transition_kernels$scale_log_scale_ancillary = state$transition_kernels$scale_log_scale_ancillary -rnorm(1, .2, .05)
          acceptance_records$scale_log_scale_ancillary =  0*acceptance_records$scale_log_scale_ancillary
        }
      }
      # sufficient -- sufficient ####
      for(i in seq(4))
      {
        new_scale_log_scale = state$params$scale_log_scale + rnorm(1, 0, .1)
        if(
          (
            + Bidart::beta_prior_ll(beta = state$params$scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                                  beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                  beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                  log_scale = new_scale_log_scale)     
            - Bidart::beta_prior_ll(beta = state$params$scale_beta, n_KL = hierarchical_model$KL$n_KL*hierarchical_model$scale_KL, 
                                    beta_mean = hierarchical_model$beta_priors$scale_beta_mean, 
                                    beta_precision =  hierarchical_model$beta_priors$scale_beta_precision, 
                                    log_scale = state$params$scale_log_scale)     
            > log(runif(1))
          )
          &(new_scale_log_scale > hierarchical_model$scale_log_scale_prior[1])
          &(new_scale_log_scale < hierarchical_model$scale_log_scale_prior[2])
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
    beta_precision = crossprod(x = sparse_chol_X, y = sparse_chol_X)
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
