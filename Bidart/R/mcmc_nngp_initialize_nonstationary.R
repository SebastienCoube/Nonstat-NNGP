#' @export
process_covariates = function(X, observed_locs, vecchia_approx, NNGP_prior_info = NULL, locs_match_matrix = F)
{
  # covariates in the observed field #
  res = list()
  # creating model matrix
  # extracting a model matrix and storing the original argument
  if(!is.null(X))
  {
    res$arg = X
    res$X = model.matrix(~., X)
  }
  # extracting a model matrix with only intercept and storing a message about the lack of original argument if no X is provided
  if(is.null(X))
  {
    res$arg = "No covariates were provided"
    res$X = matrix(model.matrix(~., as.data.frame(rep(1, nrow(observed_locs))))[,-1], nrow(observed_locs))
    colnames(res$X) = "(Intercept)"
  }
  # scaling (minus the intercept)
  res$X_mean =  apply(res$X, 2, mean) 
  res$X_sd =  apply(res$X, 2, sd) 
  res$X =  scale(res$X, scale = F)
  res$X[,1] =  1
  # pre- computing XTX
  res$solve_crossprod_X = solve(crossprod(res$X))
  res$chol_solve_crossprod_X = chol(res$solve_crossprod_X)
  # identifying  which X do not vary within location
  res$which_locs = c()
  for(i in seq(ncol(res$X))) 
  {
    if(all(duplicated(cbind(observed_locs, res$X[,i])) == vecchia_approx$duplicated_locs)) res$which_locs = c(res$which_locs, i)
  }
  res$X_locs = matrix(res$X[vecchia_approx$hctam_scol_1,res$which_locs], ncol = length(res$which_locs))
  colnames(res$X_locs) = colnames(res$X)[res$which_locs]
  res$solve_crossprod_X_locs = solve(crossprod(res$X_locs))
  res$chol_solve_crossprod_X_locs = chol(res$solve_crossprod_X_locs)
  # if needed multiplying by log-NNGP prior sparse chol
  if(!is.null(NNGP_prior_info))
  {
    res$sparse_chol_X_locs = as.matrix(NNGP_prior_info$sparse_chol %*% res$X_locs)
    res$crossprod_sparse_chol_X_locs = crossprod(res$sparse_chol_X_locs)
    res$solve_crossprod_sparse_chol_X_locs = solve(res$crossprod_sparse_chol_X_locs, tol = min(rcond(res$crossprod_sparse_chol_X_locs),.Machine$double.eps))
    res$chol_solve_crossprod_sparse_chol_X_locs = chol(res$solve_crossprod_sparse_chol_X_locs)
  }
  res
}

#' @export
beta_prior_sanity_check = function(beta_mean, beta_precision, precision_parameter_name, mean_parameter_name = NULL, X)
{
  if(is.null(mean_parameter_name))mean_parameter_name = precision_parameter_name
  # input must be list or NULL
  if(!is.list(beta_precision) & !is.null(beta_precision))stop(paste(precision_parameter_name, "_beta_precision must be a list or NULL", sep = ""))
  if(!is.list(beta_mean) & !is.null(beta_mean))stop(paste(precision_parameter_name, "_beta_precision must be a list or NULL", sep = ""))
  # expanding NULL into list of NULLs
  if(is.null(beta_precision))
  {
    if(!is.null(X))beta_precision = lapply(seq(ncol(X)+1 ), function(x)NULL)
    if(is.null(X)) beta_precision = lapply(1, function(x)NULL)
  }
  if(is.null(beta_mean))
  {
    if(!is.null(X))beta_mean = lapply(seq(ncol(X)+1 ), function(x)NULL)
    if(is.null(X)) beta_mean = lapply(1, function(x)NULL)
  }
  # checking lengths
  if( is.null(X) &(length(beta_mean)!=1 | length(beta_precision)!=1))stop(
    paste(mean_parameter_name, "_X is NULL, then ", mean_parameter_name, "_beta_mean and ", precision_parameter_name, "_beta_precision must have length equal to 1")
  )
  if(!is.null(X))
  {
    if(length(unique(c(ncol(X)+1, length(beta_mean), length(beta_precision))))>1)stop(
      paste(mean_parameter_name, "_X has ", ncol(X)+1, "columns (couting the intercept), ", mean_parameter_name, "_beta_mean has length equal to", length(beta_mean), precision_parameter_name ,"_beta_precision has length equal to", length(beta_precision), "please make those all equal")
    )
  }
  # sanity checks
  if(!all(sapply(beta_mean, is.null) == sapply(beta_precision, is.null)))stop("The NULL elements in ", mean_parameter_name, "_beta_mean and ", precision_parameter_name, "_beta_precision should match")
  if(!all(sapply(beta_mean, is.numeric) + sapply(beta_mean, is.null)))stop("all elements of ", mean_parameter_name, "_beta_mean must be numeric or NULL")
  if(!all(sapply(beta_precision, is.numeric) + sapply(beta_precision, is.null)))stop("all elements of ", precision_parameter_name, "_beta_precision must be numeric or NULL")
  if(any(unlist(beta_precision)<0))stop("all elements of ", precision_parameter_name, "_beta_precision must be greater than or equal to 0 (or NULL)")
  return(list("mean" = beta_mean, "precision" = beta_precision))
}


#' @export
beta_prior_fill_missing = function(beta_prior,  parameter_name, target_beta_0, X)
{
  res = list(); res$arg = beta_prior
  #only one regressor
  if(ncol(X$X)==1)
  {
    # first coefficient, using target_beta_0 if base model is needed 
    if(is.null(beta_prior$mean[[1]])) 
    {
      beta_prior$mean = target_beta_0
      beta_prior$precision = 1/25
      message("values 1 (intercept) of ", parameter_name, "_beta_mean and ", parameter_name, "_beta_precision were automatially set to a waek prior")
    }
    beta_prior$mean = matrix(unlist(beta_prior$mean));beta_prior$precision = matrix(unlist(beta_prior$precision))
  }
  #more than one regressor
  if(ncol(X$X)>1)
  {
    # rest of the coefficients, putting weak centered prior if needed
    for(i in seq(2, length(beta_prior$mean)))
    {
      if(is.null(beta_prior$mean[[i]])) 
      {
        beta_prior$mean[[i]] = 0
        beta_prior$precision[[i]] = (1/25)*var(X$X[,i])
        message(paste("values ", i, " of ", parameter_name, "_beta_mean and ", parameter_name, "_beta_precision were automatially set to a weak centered prior", sep = ""))
      }
    }
    # first coefficient, using target_beta_0 if base model is needed 
    if(is.null(beta_prior$mean[[1]])) 
    {
      beta_prior$mean[[1]] = sum(unlist(beta_prior$mean[-1])*X$X_mean[-1]) + target_beta_0
      beta_prior$precision[[1]] = 1/25
      message("values 1 (intercept) of ", parameter_name, "_beta_mean and ", parameter_name, "_beta_precision were automatially set to a waek prior")
    }
    # linear transform to keep up with centering of the variables
    C = diag(rep(1, length(beta_prior$mean)))
    C[1, -1] = - X$X_mean[-1]
    beta_prior$mean = solve(C, unlist(beta_prior$mean))
    beta_prior$precision = (t(C) %*% diag(beta_prior$precision) %*% C)
  }
  res$mean = beta_prior$mean
  res$precision = beta_prior$precision
  if(parameter_name %in% c("scale", "range", "anisotropy"))
  {
    res$mean_whitened = t(solve(X$chol_solve_crossprod_X_locs)) %*% beta_prior$mean
    res$precision_whitened = X$chol_solve_crossprod_X_locs %*% beta_prior$precision  %*% t(X$chol_solve_crossprod_X_locs)
  }
  if(parameter_name %in% c("noise"))
  {
    res$mean_whitened = t(solve(X$chol_solve_crossprod_X)) %*% beta_prior$mean
    res$precision_whitened = X$chol_solve_crossprod_X %*% beta_prior$precision  %*% t(X$chol_solve_crossprod_X)
  }
  res
}

merge_beta_prior_range_precision = function(beta_prior_range, beta_prior_anisotropy)
{
  res = list()
  res$arg_range = beta_prior_range$arg
  res$arg_anisotropy = beta_prior_anisotropy$arg
  M = matrix(c(sqrt(2)/2, sqrt(2)/2, 0, sqrt(2)/2, -sqrt(2)/2, 0, 0,0,1),3) %x% diag(nrow = length(beta_prior_range$mean))
  res$mean = M %*% c(beta_prior_range$mean, beta_prior_anisotropy$mean, beta_prior_anisotropy$mean)
  res$precision = t(solve(M)) %*% as.matrix(Matrix::bdiag(beta_prior_range$precision, beta_prior_anisotropy$precision, beta_prior_anisotropy$precision)) %*% solve(M)
  res$mean_whitened = M %*% c(beta_prior_range$mean_whitened, beta_prior_anisotropy$mean_whitened, beta_prior_anisotropy$mean_whitened)
  res$precision_whitened = t(solve(M)) %*% Matrix::bdiag(beta_prior_range$precision_whitened, beta_prior_anisotropy$precision_whitened, beta_prior_anisotropy$precision_whitened) %*% solve(M)
  return(res)
}


#' @export
NNGP_prior = function(range_param, vecchia_approx, log_NNGP_nu, log_NNGP_matern_covfun, locs)
{
  if(is.null(range_param))return(NULL) 
  if(!is.null(range_param)) 
  {
    res = list()
    res$range = range_param
    res$smoothness = log_NNGP_nu
    compressed_sparse_chol = GpGp::vecchia_Linv(covparms = c(1, range_param, log_NNGP_nu, 0), covfun_name = log_NNGP_matern_covfun, locs = locs, NNarray = vecchia_approx$NNarray)
    res$sparse_chol = Matrix::sparseMatrix(
      i = vecchia_approx$sparse_chol_row_idx, 
      j = vecchia_approx$sparse_chol_column_idx, 
      x = compressed_sparse_chol[vecchia_approx$NNarray_non_NA], 
      triangular = T)
    res$precision_diag =    as.vector((compressed_sparse_chol[vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
    return(res)
  }
}



#' @export
mcmc_nngp_initialize_nonstationary = 
  function(observed_locs = NULL, #spatial locations
           observed_field = NULL, # Response variable
           X = NULL, # Covariates per observation
           m = 5, #number of Nearest Neighbors
           nu =1.5, #Matern smoothness
           reordering = "maxmin", #Reordering
           covfun = "exponential_isotropic", response_model = "Gaussian", # covariance model and response model
           noise_X = NULL, noise_beta_mean = NULL, noise_beta_precision = NULL, noise_range = NULL, # range for latent field of parameters, if NULL no latent field
           scale_X = NULL, scale_beta_mean = NULL, scale_beta_precision = NULL, scale_range = NULL, # range for latent field of parameters, if NULL no latent field
           range_X = NULL, range_beta_mean = NULL, range_beta_precision = NULL, range_range = NULL, anisotropy_beta_precision = NULL, # range for latent field of parameters, if NULL no latent field
           log_NNGP_matern_covfun = NULL, log_NNGP_nu = NULL, # covariance function for the hyperpriors
           n_chains = 3,  # number of MCMC chains
           seed = 1
  )
  {
    
    # time  
    t_begin = Sys.time()
    # seed
    set.seed(seed)
    # cleansing RAM
    gc()
    #################
    # Sanity checks #
    #################
    
    # allowed functions
    if(!reordering[1] %in% c("maxmin", "random", "coord", "dist_to_point", "middleout"))stop("reordering should be chosen among : maxmin, random, coord, dist_to_point, middleout")
    if(response_model!="Gaussian") stop("response_model should be chosen among : Gaussian")
    if(!covfun %in% c("exponential_isotropic", "exponential_sphere",  "exponential_spacetime", "exponential_spheretime", "exponential_anisotropic2D", 
                      "matern_isotropic",      "matern_sphere",       "matern_spacetime",      "matern_spheretime",      "matern_anisotropic2D", 
                      "nonstationary_exponential_isotropic", "nonstationary_exponential_isotropic_sphere", 
                      "nonstationary_matern_isotropic",      "nonstationary_matern_isotropic_sphere", 
                      "nonstationary_exponential_anisotropic", "nonstationary_exponential_anisotropic_sphere"))stop(
                        "covfun should be chosen among : 
                        exponential_isotropic,     exponential_sphere, exponential_spacetime,   exponential_spheretime,     exponential_anisotropic2D, 
                        matern_isotropic, matern_sphere, matern_spacetime, matern_spheretime, matern_anisotropic2D, 
                        nonstationary_exponential_isotropic, nonstationary_exponential_isotropic_sphere" , "nonstationary_matern_isotropic", 
                        "nonstationary_matern_isotropic_sphere", "nonstationary_exponential_anisotropic", "nonstationary_exponential_anisotropic_sphere")
    # format
    if(!is.matrix(observed_locs))stop("observed_locs should be a matrix")
    if(!is.vector(observed_field))stop("observed_field should be a vector")
    if(!is.data.frame(X            ) & !is.null(X))stop("X should be a data.frame or nothing")
    if(!is.data.frame(noise_X) & !is.null(noise_X))stop("noise_X should be a data.frame or nothing")
    if(!is.data.frame(scale_X) & !is.null(scale_X))stop("scale_X should be a data.frame or nothing")
    if(!is.data.frame(range_X) & !is.null(range_X))stop("range_X should be a data.frame or nothing")
    # positivity of ranges
    if(any(!noise_range>0)) stop("noise_range should be greater than 0")
    if(any(!scale_range>0)) stop("scale_range should be greater than 0")
    if(any(!range_range>0)) stop("range_range should be greater than 0")
    #length of observations
    if(
      length(unique(c(
        length(observed_field),
        nrow(observed_locs),
        nrow(X),
        nrow(scale_X),
        nrow(noise_X),
        nrow(range_X)))
      )!=1) stop(
        paste("Lengths are not matching : observed_field has", length(observed_field), "observations,",
              "observed_locs has", nrow(observed_locs), "rows,", 
              "X has", nrow(X), "rows,", 
              "scale_X has", nrow(scale_X), "rows,", 
              "noise_X has", nrow(noise_X), "rows,", 
              "range_X has", nrow(range_X), "rows."
        )
      )
    if ((nu != 1)&(nu != 2)&(nu != 1.5))stop("only nu = 1, 1.5, 2 for now (1.5 recommended)")
    if ((nu == 1)|(nu == 2))message("nu = 1.5 is much faster")
    #nonstationary range is provided field and/or covariates
    #if((length(grep("nonstationary", covfun))!=0) & is.null(range_X)&(is.null(range_range)))stop(paste("Nonstationary covariance function", covfun, "should be provided with at least :  1 or 2 range parameter (s) for the random field (range_range)  or/and : one set of covariates (range_X)"))
    #stationary range is not provided field or covariates
    if((length(grep("nonstationary", covfun))==0) & (!is.null(range_X)|(!is.null(range_range))))stop(paste("Stationary covariance function", covfun, "should have no range parameter for the random field (range_range)  and no covariates (range_X)"))
    # useless log_NNGP_matern_covfun provided
    if(!is.null(log_NNGP_matern_covfun)&(is.null(scale_range)&is.null(noise_range)&is.null(range_range)&is.null(log_NNGP_nu))) stop("A covariance function was indicated for the log-NNGP prior. Indicate hyperprior range/smoothness arguments or remove the covariance function.")
    # no log_NNGP_matern_covfun provided
    if(is.null(log_NNGP_matern_covfun)&(!is.null(scale_range)|!is.null(noise_range)|!is.null(range_range)|!is.null(log_NNGP_nu))) stop("No covariance function was indicated for the log-NNGP prior. Indicate a covariance function or remove hyperprior range/smoothness arguments.")
    # bad log_NNGP_matern_covfun provided
    if(!is.null(log_NNGP_matern_covfun)){if(!log_NNGP_matern_covfun %in% c("matern_spheretime", "matern_isotropic", "matern_spacetime", "matern_sphere"))stop("log_NNGP_matern_covfun should be chosen among : matern_spheretime, matern_isotropic, matern_spacetime, matern_sphere")}
    # message about discrepancy
    if (length(grep("sphere", log_NNGP_matern_covfun)) != length(grep("sphere", covfun))) message("The hyperprior covariance is on the sphere while the nonstationary covariance is on the plane (or the other way around)")

    # prior length and NULL patterns
    noise_beta_prior = beta_prior_sanity_check(beta_mean = noise_beta_mean, beta_precision = noise_beta_precision, precision_parameter_name = "noise", X = noise_X)
    range_beta_prior = beta_prior_sanity_check(beta_mean = range_beta_mean, beta_precision = range_beta_precision, precision_parameter_name = "range", X = range_X)
    scale_beta_prior = beta_prior_sanity_check(beta_mean = scale_beta_mean, beta_precision = scale_beta_precision, precision_parameter_name = "scale", X = scale_X)
    if((length(grep("anisotropic", covfun))==1) & (length(grep("nonstationary", covfun))==1))anisotropy_beta_prior = beta_prior_sanity_check(beta_mean = anisotropy_beta_precision, beta_precision = anisotropy_beta_precision, mean_parameter_name = "range", precision_parameter_name = "anisotropy", X = range_X)
    
    ###############
    # Re-ordering #
    ###############
    # remove duplicated locations
    duplicated_locs = duplicated (observed_locs)
    locs = observed_locs[duplicated_locs==F,]
    if(reordering[1] == "maxmin")locs_reordering = GpGp::order_maxmin(locs, lonlat = !identical(grep("sphere", covfun) , integer(0)))
    if(reordering[1] == "random")locs_reordering = sample(seq(nrow(locs)))
    if(reordering[1] == "coord")locs_reordering = GpGp::order_coordinate(locs = locs, coordinate = as.numeric(reordering[2]))
    if(reordering[1] == "dist_to_point")locs_reordering = GpGp::order_dist_to_point(locs, loc0 = as.numeric(reordering[2]), lonlat = !identical(grep("sphere", covfun) , integer(0)))
    if(reordering[1] == "middleout")locs_reordering = GpGp::order_middleout(locs, lonlat = !identical(grep("sphere", covfun) , integer(0)))
    locs = locs[locs_reordering,]
    # extracting number of locations as shortcut
    n = nrow(locs)
    
    #########################
    # Vecchia approximation #
    #########################
    
    # This object gathers the NNarray table used by GpGp package and related objects
    
    vecchia_approx = list()
    # saving reordering method
    vecchia_approx$reordering = reordering
    # storing numbers
    vecchia_approx$n_locs = n
    vecchia_approx$n_obs = length(observed_field)
    # matching observed locations with reordered, unrepeated locations
    locs_match = match(split(observed_locs, row(observed_locs)), split(locs, row(locs)))
    vecchia_approx$locs_match = locs_match
    vecchia_approx$locs_match_matrix = Matrix::sparseMatrix(i = vecchia_approx$locs_match, j = seq(vecchia_approx$n_obs), x = 1)
    # doing reversed operation : for a given unrepeated location, tell which observations correspond
    vecchia_approx$hctam_scol = split(seq(vecchia_approx$n_obs), locs_match)
    vecchia_approx$hctam_scol_1 = sapply(vecchia_approx$hctam_scol, function(x)x[1])
    # count how many observations correspond to one location
    vecchia_approx$obs_per_loc = unlist(sapply(vecchia_approx$hctam_scol, length))
    #extracting NNarray =  nearest neighbours for Vecchia approximation
    vecchia_approx$NNarray = GpGp::find_ordered_nn(locs, m)
    
    #computations from vecchia_approx$NNarray in order to create sparse Cholesky using Matrix::sparseMatrix
    #non_NA indices from vecchia_approx$NNarray
    vecchia_approx$NNarray_non_NA = !is.na(vecchia_approx$NNarray)
    #column idx of the uncompressed sparse Cholesky factor
    vecchia_approx$sparse_chol_column_idx = vecchia_approx$NNarray[vecchia_approx$NNarray_non_NA]
    #row idx of the uncompressed sparse Cholesky factor
    vecchia_approx$sparse_chol_row_idx = row(vecchia_approx$NNarray)[vecchia_approx$NNarray_non_NA]
    # adjacency matrix of MRF
    vecchia_approx$MRF_adjacency_mat =  Matrix::crossprod(Matrix::sparseMatrix(i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, x = 1))
    # stupid trick to coerce adjacency matrix format...
    vecchia_approx$MRF_adjacency_mat[1, 2] = 0
    vecchia_approx$MRF_adjacency_mat[1, 2] = 1
    vecchia_approx$MRF_adjacency_mat@x = rep(1, length(vecchia_approx$MRF_adjacency_mat@x))
    #vecchia_approx$coloring = naive_greedy_coloring(vecchia_approx$MRF_adjacency_mat)
    # duplicated locs
    vecchia_approx$duplicated_locs = duplicated_locs
    # partition of locs for field update
    vecchia_approx$locs_partition = sapply(seq(n/10000+1, 2*(n/10000)+1 ), function(i)
      kmeans(locs, centers = i)$cluster)
    
    #################################
    # Info about hierarchical model #
    #################################
    
    hierarchical_model = list()
    hierarchical_model$response_model = response_model
    hierarchical_model$covfun = covfun
    hierarchical_model$nu = nu
    hierarchical_model$beta_priors = list()
      
    # log-NNGP hyperpriors
    hierarchical_model$hyperprior_covariance = list()
    hierarchical_model$hyperprior_covariance$log_NNGP_matern_covfun = log_NNGP_matern_covfun
    hierarchical_model$hyperprior_covariance$range_NNGP_prior = NNGP_prior(range_range, vecchia_approx, log_NNGP_nu, log_NNGP_matern_covfun, locs = locs)
    hierarchical_model$hyperprior_covariance$scale_NNGP_prior = NNGP_prior(scale_range, vecchia_approx, log_NNGP_nu, log_NNGP_matern_covfun, locs = locs)
    hierarchical_model$hyperprior_covariance$noise_NNGP_prior = NNGP_prior(noise_range, vecchia_approx, log_NNGP_nu, log_NNGP_matern_covfun, locs = locs)
    

    ##############
    # covariates #
    ##############
    covariates = list()
    
    covariates$X = process_covariates(X, observed_locs, vecchia_approx)  
    
    covariates$range_X = process_covariates(range_X, observed_locs, vecchia_approx, hierarchical_model$hyperprior_covariance$range_NNGP_prior)
    if(!identical(covariates$range_X$which_locs, seq(ncol(covariates$range_X$X_locs))))stop("The covariates range_X cannot vary within one spatial location of observed_locs")
    
    covariates$scale_X = process_covariates(scale_X, observed_locs, vecchia_approx, hierarchical_model$hyperprior_covariance$scale_NNGP_prior)
    if(!identical(covariates$scale_X$which_locs, seq(ncol(covariates$scale_X$X))))stop("The covariates scale_X cannot vary within one spatial location of observed_locs")
    
    covariates$noise_X = process_covariates(noise_X, observed_locs, vecchia_approx, hierarchical_model$hyperprior_covariance$noise_NNGP_prior, locs_match_matrix = T)
    
    
    ####################
    # range beta prior #
    ####################
    # prior 
    # using high-range model as base model
    hierarchical_model$beta_priors$range_beta = beta_prior_fill_missing(beta_prior = range_beta_prior, parameter_name = "range", target_beta_0 = log(max(dist(locs[sample(seq(vecchia_approx$n_locs), 1000),]))), X = covariates$range_X)
    if((length(grep("anisotropic", covfun))==1) & (length(grep("nonstationary", covfun))==1))
    {
      anisotropy_beta_prior$mean = lapply(anisotropy_beta_prior$mean, function(x)
        {
        if(is.numeric(x)) return(0)
        if(is.null(x))return(NULL)
      })
      anisotropy_beta_prior = beta_prior_fill_missing(beta_prior = anisotropy_beta_prior, parameter_name = "anisotropy", target_beta_0 = 0, X = covariates$range_X)
      hierarchical_model$beta_priors$range_beta = merge_beta_prior_range_precision(hierarchical_model$beta_priors$range_beta, anisotropy_beta_prior)
    }
    
    ################
    # Chain states #
    ################
    # for each chain, creating sub-lists in order to stock all the stuff that is related to one chain, including : 
    # transition_kernel_sd : a list that stocks the (current) automatically-tuned transition kernels standard deviations
    # params : a list that stocks the (current) parameters of the model, including covariance parameters, the value of the sampled field, etc
    # records : records of the MCMC iterations
    
    states = list()
    # creating sub-list : each sub-list is a chain
    for(i in seq(n_chains))
    {
      states[[paste("chain", i, sep = "_")]] = list()
      #parameters of interest to the model
      states[[i]]$params = list()
      #stuff that is useful for computations
      states[[i]]$sparse_chol_and_stuff = list()
      # HMC momenta
      states[[i]]$momenta = list()
      
      #########
      # Range #
      #########
      
      # mean log-range (only parameter in the case of stationary model)
      if(covfun %in% c("exponential_spacetime"                     , "matern_spacetime"                     )) states[[i]]$params$range_beta = matrix(c(sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),-ncol(locs)]))     )-log(seq(50, 500, 1)), 1         ), sample(log(max(dist(locs[sample(seq(vecchia_approx$n_locs), 100),ncol(locs)])))-log(seq(5, 10, 1)), 1))           , nrow = 1)
      if(covfun %in% c("exponential_isotropic"                     , "matern_isotropic"                     )) states[[i]]$params$range_beta = matrix(  sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), 1         )                                                                                                                    , nrow = 1)
      if(covfun %in% c("exponential_sphere"                        , "matern_sphere"                        )) states[[i]]$params$range_beta = matrix(  sample(log(max(fields::rdist.earth(locs[sample(seq(vecchia_approx$n_locs), 1000),           ]))/6000)-log(seq(50, 500, 1)), 1         )                                                                                                            , nrow = 1)
      if(covfun %in% c("exponential_spacetime_sphere"              , "matern_spacetime_sphere"              )) states[[i]]$params$range_beta = matrix(c(sample(log(max(fields::rdist.earth(locs[sample(seq(vecchia_approx$n_locs), 1000), c(1, 2)   ]))/6000)-log(seq(50, 500, 1)), 1         ), sample(log(max(dist(locs[sample(seq(vecchia_approx$n_locs), 100),ncol(locs)])))-log(seq(5, 10, 1)), 1))   , nrow = 1)
      if(covfun %in% c("nonstationary_exponential_isotropic"       , "nonstationary_matern_isotropic"       )) states[[i]]$params$range_beta = matrix(  sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), 1         )                                                                                                                    , nrow = 1)
      if(covfun %in% c("nonstationary_exponential_isotropic_sphere", "nonstationary_matern_isotropic_sphere")) states[[i]]$params$range_beta = matrix(  sample(log(max(fields::rdist.earth(locs[sample(seq(vecchia_approx$n_locs), 1000),           ]))/6000)-log(seq(50, 500, 1)), 1         )                                                                                                            , nrow = 1)
      if(covfun %in% c("exponential_anisotropic2D"                 , "matern_anisotropic2D"                 )) states[[i]]$params$range_beta = matrix(c(sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), ncol(locs)),                       rep(0, ncol(locs)*(ncol(locs)+1)/2-ncol(locs)))                                             , nrow = 1)
      if(covfun == "nonstationary_exponential_anisotropic")                                                    states[[i]]$params$range_beta = matrix(c(sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), ncol(locs)),                       rep(0, ncol(locs)*(ncol(locs)+1)/2-ncol(locs)))                                             , nrow = 1)
      if(covfun == "nonstationary_exponential_anisotropic_sphere")                                             states[[i]]$params$range_beta = matrix(c(sample(log(max(fields::rdist.earth(locs[sample(seq(vecchia_approx$n_locs), 1000), c(1, 2)   ]))/6000)-log(seq(50, 500, 1)), 2         ), rep(1, ncol(locs)-2), rep(0, ncol(locs)*(ncol(locs)+1)/2-ncol(locs)))                                     , nrow = 1)
      # initiate null regressors for other covariates                                                                                                                                          
      if(!is.null(covariates$range_X$X_locs))states[[i]]$params$range_beta = rbind(states[[i]]$params$range_beta, matrix(0#rnorm((ncol(covariates$range_X$X_locs)-1) * length(states[[i]]$params$range_beta))
                                                                                                                    , ncol(covariates$range_X$X_locs)-1, length(states[[i]]$params$range_beta)))
      row.names(states[[i]]$params$range_beta) = colnames(covariates$range_X$X_locs)
      states[[i]]$momenta$range_beta_ancillary = matrix(rnorm(length(states[[i]]$params$range_beta)), nrow(states[[i]]$params$range_beta))
      states[[i]]$momenta$range_beta_sufficient = matrix(rnorm(length(states[[i]]$params$range_beta)), nrow(states[[i]]$params$range_beta))
        
      # field 
      if(!is.null(hierarchical_model$hyperprior_covariance$range))
      {
        # scale is given by a PD matrix whose Cholesky coefficients are given here
        if(covfun == "nonstationary_exponential_isotropic" | covfun == "nonstationary_exponential_isotropic_sphere")
        {
          states[[i]]$params$range_field = matrix(rnorm(1 * vecchia_approx$n_locs), ncol = 1)
          states[[i]]$params[["range_log_scale"]] = rnorm(1, -4, 1)
          states[[i]]$momenta[["range_log_scale_ancillary"]] = rnorm(1)
          states[[i]]$momenta[["range_log_scale_sufficient"]] = rnorm(1)
        }
        if(covfun == "nonstationary_exponential_anisotropic" | covfun == "nonstationary_exponential_anisotropic_sphere") 
        {
          states[[i]]$params$range_field = matrix(rnorm(ncol(locs)*(ncol(locs)+1)/2 * vecchia_approx$n_locs), ncol = ncol(locs)*(ncol(locs)+1)/2)
          states[[i]]$params[["range_log_scale"]] =             c(-2, -2, -2, 0, 0, 0)
          states[[i]]$momenta[["range_log_scale_ancillary"]] =  rnorm(ncol(states[[i]]$params$range_field)*(ncol(states[[i]]$params$range_field)+1)/2)
          states[[i]]$momenta[["range_log_scale_sufficient"]] = rnorm(ncol(states[[i]]$params$range_field)*(ncol(states[[i]]$params$range_field)+1)/2)
        }
        # range field
        states[[i]]$params$range_field =  as.matrix(Matrix::solve(hierarchical_model$hyperprior_covariance$range$sparse_chol, states[[i]]$params$range_field)) %*% chol(Bidart::expmat(states[[i]]$params$range_log_scale))
        states[[i]]$momenta$range_field_ancillary  =  matrix(rnorm(length(states[[i]]$params$range_field)), nrow = nrow(states[[i]]$params$range_field))
        states[[i]]$momenta$range_field_sufficient =  matrix(rnorm(length(states[[i]]$params$range_field)), nrow = nrow(states[[i]]$params$range_field))
      }
      
      ######################
      # Transition kernels #  
      ######################
      
      # Starting points for transition kernels, will be adaptively tuned
      states[[i]]$transition_kernels = list()
      # Transition kernel state
      # transition kernel variance is given as the log
      # can be used in both stationary and nonstationary cases respectively as a random walk Metropolis or MALA step size 
      # have an ancillary and a sufficient version when applicable
      # range
      states[[i]]$transition_kernels$range_scale_joint = -4
      states[[i]]$transition_kernels$range_field_sufficient_mala = -2
      states[[i]]$transition_kernels$range_field_ancillary_mala  = -2
      states[[i]]$transition_kernels$range_log_scale_sufficient = -4
      states[[i]]$transition_kernels$range_log_scale_ancillary =  -4
      states[[i]]$transition_kernels$range_beta_sufficient = -4
      states[[i]]$transition_kernels$range_beta_ancillary  = -4
      # scale
      states[[i]]$transition_kernels$scale_field_sufficient_mala = -3
      states[[i]]$transition_kernels$scale_field_ancillary_mala = -3
      states[[i]]$transition_kernels$scale_beta_sufficient_mala = -2
      states[[i]]$transition_kernels$scale_beta_ancillary_mala  = -2
      states[[i]]$transition_kernels$scale_sufficient_beta = rep(0, ncol(covariates$scale_X$X_locs))
      states[[i]]$transition_kernels$scale_ancillary_beta = rep(0, ncol(covariates$scale_X$X_locs))
      states[[i]]$transition_kernels$scale_sufficient_log_scale = -1
      states[[i]]$transition_kernels$scale_ancillary_log_scale = -1
      # noise variance
      states[[i]]$transition_kernels$noise_field_mala = -3
      states[[i]]$transition_kernels$noise_beta = rep(0, ncol(covariates$noise_X$X))
      states[[i]]$transition_kernels$noise_beta_mala = -1
      states[[i]]$transition_kernels$noise_log_scale = -3
      # MALA for latent field 
      # states[[i]]$transition_kernels$latent_field_mala = -2
      
      ##############################################################################
      # Linear regression coefficients, noise variance, and gaussian process scale #
      ##############################################################################
      # Cases are treated following data model
      
      #################
      # Gaussian case #
      #################
      if(response_model == "Gaussian")
      {
        ###################################
        # Linear regression coefficients  #
        ###################################
        # Naive OLS to determine starting values for beta, the field, 
        naive_ols =  lm(observed_field~covariates$X$X-1)
        #starting points for regression coeffs
        perturb = t(chol(vcov(naive_ols)))%*%rnorm(length(naive_ols$coefficients))
        states[[i]]$params[["beta"]] = naive_ols$coefficients + perturb
        row.names(states[[i]]$params[["beta"]]) = colnames(covariates$X$X)
        # Residuals of the OLS model that have to be explained by the latent field and the noise
        states[[i]]$sparse_chol_and_stuff$lm_fit = as.vector(covariates$X$X%*%matrix(states[[i]]$params[["beta"]], ncol = 1))
        states[[i]]$sparse_chol_and_stuff$lm_fit_locs = as.vector(covariates$X$X_locs%*%matrix(states[[i]]$params[["beta"]][covariates$X$which_locs], ncol = 1))
        states[[i]]$sparse_chol_and_stuff$lm_residuals = as.vector(observed_field-  states[[i]]$sparse_chol_and_stuff$lm_fit)
        ##################
        # Noise variance #
        ##################
        # beta is just an intercept in stationary case
        states[[i]]$params$noise_beta    = matrix(rep(0, ncol(covariates$noise_X$X)), ncol = 1) #random starting values
        states[[i]]$params$noise_beta[1] = log(var(states[[i]]$sparse_chol_and_stuff$lm_residuals)) - log(2) +rnorm(1, 0, .5) # setting sensible value for the intercept
        row.names(states[[i]]$params$noise_beta) = colnames(covariates$noise_X$X)
        states[[i]]$momenta$noise_beta = rnorm(ncol(covariates$noise_X$X))
        # field  when required
        if(!is.null(hierarchical_model$hyperprior_covariance$noise))
        {
          states[[i]]$params$noise_log_scale = rnorm(1, -4, 1) 
          states[[i]]$params$noise_field = exp(.5*states[[i]]$params$noise_log_scale) * as.matrix(Matrix::solve(hierarchical_model$hyperprior_covariance$noise$sparse_chol,  rnorm(vecchia_approx$n_locs)))
          states[[i]]$momenta$noise_field = rnorm(vecchia_approx$n_locs)
        }
        # effective variance field, shall be used in density computations
        states[[i]]$sparse_chol_and_stuff$noise = Bidart::variance_field(states[[i]]$params$noise_beta, states[[i]]$params$noise_field[vecchia_approx$locs_match], covariates$noise_X$X)
        # prior 
        if(i==1) 
        {
          # using high-noise model as base model
          hierarchical_model$beta_priors$noise_beta = beta_prior_fill_missing(beta_prior = noise_beta_prior, parameter_name = "noise", target_beta_0 = log(var(states[[i]]$sparse_chol_and_stuff$lm_residuals)), X = covariates$noise_X)
        }
        
        #########
        # Scale #
        #########
        # beta is just an intercept in stationary case
        states[[i]]$params$scale_beta    = matrix(0, ncol(covariates$scale_X$X_locs), ncol = 1) #random starting values
        states[[i]]$momenta$scale_beta_ancillary = rnorm(ncol(covariates$scale_X$X_locs))
        states[[i]]$momenta$scale_beta_sufficient = rnorm(ncol(covariates$scale_X$X_locs))
        states[[i]]$params$scale_beta[1] = log(var(states[[i]]$sparse_chol_and_stuff$lm_residuals)) - log(2) + rnorm(1, 0, .5) # setting sensible value for the intercept
        row.names(states[[i]]$params$scale_beta) = colnames(covariates$scale_X$X_locs)
        
        # prior 
        if(i==1) 
        {
          hierarchical_model$beta_priors$scale_beta = beta_prior_fill_missing(beta_prior = scale_beta_prior, parameter_name = "scale", target_beta_0 = log(var(states[[i]]$sparse_chol_and_stuff$lm_residuals)), X = covariates$scale_X)
        }
        
        
        # field  when required
        if(!is.null(hierarchical_model$hyperprior_covariance$scale))
        {
          states[[i]]$params$scale_log_scale = rnorm(1, -4, 1)
          states[[i]]$transition_kernels$scale_log_scale = rnorm(1)
          states[[i]]$params$scale_field = exp(.5*states[[i]]$params$scale_log_scale) * as.matrix(Matrix::solve(hierarchical_model$hyperprior_covariance$scale$sparse_chol,  rnorm(vecchia_approx$n_locs)))
          states[[i]]$momenta$scale_field_ancillary =  rnorm(length(states[[i]]$params$scale_field))
          states[[i]]$momenta$scale_field_sufficient = rnorm(length(states[[i]]$params$scale_field))
        }
        # effective variance field, shall be used in density computations
        states[[i]]$sparse_chol_and_stuff$scale = Bidart::variance_field(states[[i]]$params$scale_beta, states[[i]]$params$scale_field, covariates$scale_X$X_locs)
        #######
      }
      ####################
      # NNGP sparse chol #
      ####################
      states[[i]]$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, range_X = covariates$range_X$X_locs, range_beta = states[[i]]$params$range_beta, range_field = states[[i]]$params$range_field, NNarray = vecchia_approx$NNarray, locs = locs, nu = nu)
      states[[i]]$sparse_chol_and_stuff$sparse_chol = Matrix::sparseMatrix(x =  states[[i]]$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, triangular = T)
      states[[i]]$sparse_chol_and_stuff$precision_diag = as.vector((states[[i]]$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
      ################
      # Latent field #
      ################
      states[[i]]$params$field = sqrt(states[[i]]$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(states[[i]]$sparse_chol_and_stuff$sparse_chol, rnorm(vecchia_approx$n_locs)))
      states[[i]]$momenta$field = rnorm(length(states[[i]]$params$field))
    }
    
    #######################
    # Chain records setup #
    #######################
    
    # records is a list that stocks the recorded parameters of the model, including covariance parameters, the value of the sampled field, etc. In terms of RAM, those are the biggest bit !
    # iteration is a 2-colums matrix that records the iteration at the end of each chains join and the associated CPU time
    records = list()
    for(i in seq(n_chains))
    {
      records[[paste("chain", i, sep = "_")]] = list()
      records[[i]] =  sapply(states[[i]]$params, function(x)NULL)
      
    }
    iterations = list()
    iterations$checkpoints =  matrix(c(0, as.numeric(Sys.time()-t_begin, unit = "mins")), ncol = 2)
    colnames(iterations$checkpoints) = c("iteration", "time")
    iterations$thinning = c()
    
    ##########
    
    message(paste("Setup done,", as.numeric(Sys.time()- t_begin, units = "secs"), "s elapsed" ))
    return(list("data" = list("locs" = locs, "observed_field" = observed_field, "observed_locs" = observed_locs, "covariates" = covariates), "hierarchical_model" = hierarchical_model, 
                "vecchia_approx" = vecchia_approx, "states" = states, "records" = records, "t_begin" = t_begin, "seed" = seed, "iterations" = iterations))
  }  