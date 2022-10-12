#' @export
process_covariates = function(X, observed_locs, vecchia_approx, KL = NULL, use_KL = F)
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
  X_ = res$X
  if(use_KL) X_ = cbind(res$X, KL$basis %*% diag(KL$KL_decomposition$d))
  # pre- computing XTX
  crossprod_X = crossprod(X_)
  res$chol_crossprod_X = chol(crossprod_X)
  res$n_regressors = ncol(res$X)
  # identifying  which X do not vary within location
  res$which_locs = c()
  for(i in seq(ncol(res$X))) 
  {
    if(all(duplicated(cbind(observed_locs, res$X[,i])) == vecchia_approx$duplicated_locs)) res$which_locs = c(res$which_locs, i)
  }
  res$X_locs = matrix(res$X[vecchia_approx$hctam_scol_1,res$which_locs], ncol = length(res$which_locs))
  colnames(res$X_locs) = colnames(res$X)[res$which_locs]
  X_locs_ = res$X_locs
  if(use_KL)X_locs_ = cbind(X_locs_, (KL$basis %*% diag(KL$KL_decomposition$d))[vecchia_approx$hctam_scol_1,])
  crossprod_X_locs = crossprod(X_locs_)
  res$chol_crossprod_X_locs = chol(crossprod_X_locs)
  res
}


#' Initializes a list containing the observations, the hierarchical model, and the MCMC chains.
#' 
#' @param observed_locs a matrix of spatial coordinates where observations are done
#' @param observed_field a vector of observations of the interest variable
#' @param X a data.frame of covariates explaining the interest variable through fixed linear effects
#' @param m number of nearest neighbors to do Vecchia's approximation
#' @param nu Matern smoothness, either 0.5 or 1.5
#' @param reordering indicates the ordering heuristic, either "random", "maxmin", "middleout", list("coord", which coordinate is used), list("dist_to_point", coordinate of the starting point)
#' @param covfun covariance function for the Gaussian latent field, either "exponential_isotropic", "exponential_sphere", "exponential_anisotropic2D", "matern_isotropic",      "matern_sphere",      "matern_anisotropic2D", "nonstationary_exponential_isotropic", "nonstationary_exponential_isotropic_sphere", "nonstationary_matern_isotropic",      "nonstationary_matern_isotropic_sphere", "nonstationary_exponential_anisotropic", "nonstationary_exponential_anisotropic_sphere"
#' @param noise_X a data.frame of covariates explaining the Gaussian noise variance through fixed linear effects
#' @param scale_X a data.frame of covariates explaining the Gaussian process marginal variance through fixed linear effects
#' @param range_X a data.frame of covariates explaining the Gaussian process range through fixed linear effects
#' @param noise_beta_mean vector indicating the prior mean for the regression coefficients explaining the Gaussian noise variance through fixed linear effects. Filled automatically if NULL. The number of rows is equal to the number of variables in X_noise after adding an intercept and expanding the factors, the number of columns is 1.
#' @param scale_beta_mean vector indicating the prior mean for the regression coefficients explaining the Gaussian process marginal variance through fixed linear effects. Filled automatically if NULL. The number of rows is equal to the number of variables in X_scale after adding an intercept, the number of columns is 1.
#' @param range_beta_mean vector (3-columns matrix in anisotropy case) indicating the prior mean for the regression coefficients explaining the Gaussian process range through fixed linear effects. Filled automatically if NULL.The number of rows is equal to the number of variables in X_range after adding an intercept, the number of columns is 1 if isotropic function or 3 if anisotropic function.
#' @param noise_beta_precision matrix for the prior precision for the regression coefficients explaining the Gaussian noise variance through fixed linear effects. Filled automatically if NULL.
#' @param scale_beta_precision matrix for the prior precision for the regression coefficients explaining the Gaussian process marginal variance through fixed linear effects. Filled automatically if NULL.
#' @param range_beta_precision matrix for the prior precision for the regression coefficients explaining the Gaussian process range through fixed linear effects. Filled automatically if NULL. In anisotropic case, the matrix has size 3*n_var, each 3-block indicating the precision for a determinant-direction-direction 3-uplet: intecrept-det, intercept-dir, intercept-dir, V1-det, V1-dir, V1-dir, etc...
 
mcmc_nngp_initialize_nonstationary = 
  function(observed_locs = NULL, #spatial locations
           observed_field = NULL, # Response variable
           X = NULL, # Covariates per observation
           m = 5, #number of Nearest Neighbors
           nu =1.5, #Matern smoothness
           reordering = "maxmin", #Reordering
           covfun = "exponential_isotropic",
           noise_X = NULL, noise_beta_mean = NULL, noise_beta_precision = NULL, noise_KL = F, noise_log_scale_prior = NULL, 
           scale_X = NULL, scale_beta_mean = NULL, scale_beta_precision = NULL, scale_KL = F, scale_log_scale_prior = NULL, 
           range_X = NULL, range_beta_mean = NULL, range_beta_precision = NULL, range_KL = F, range_log_scale_prior = NULL, 
           KL = NULL, 
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
    allowed_covfuns = c("exponential_isotropic", "exponential_sphere", 
                        "exponential_anisotropic2D", 
      "matern_isotropic",      "matern_sphere",      
      "matern_anisotropic2D", 
      "nonstationary_exponential_isotropic", "nonstationary_exponential_isotropic_sphere", 
      "nonstationary_matern_isotropic",      "nonstationary_matern_isotropic_sphere", 
      "nonstationary_exponential_anisotropic", "nonstationary_exponential_anisotropic_sphere", 
      "nonstationary_matern_anisotropic", "nonstationary_matern_anisotropic_sphere"
      )
    if(!covfun %in% allowed_covfuns)stop(do.call("paste", as.list(c("covfun should be chosen among :", allowed_covfuns))))
    # format
    if(!is.matrix(observed_locs))stop("observed_locs should be a matrix")
    if(!is.vector(observed_field))stop("observed_field should be a vector")
    
    if(!is.data.frame(X) & !is.null(X))stop("X should be a data.frame or nothing")
    
    if(!is.data.frame(noise_X) & !is.null(noise_X))stop("noise_X should be a data.frame or nothing")
    if(!is.data.frame(scale_X) & !is.null(scale_X))stop("scale_X should be a data.frame or nothing")
    if(!is.data.frame(range_X) & !is.null(range_X))stop("range_X should be a data.frame or nothing")
    
    if((is.null(KL)) & (noise_KL | range_KL | scale_KL))stop("either noise_KL, range_KL, or scale_KL is TRUE, while nothing was provided for KL")
    
    if(is.matrix(noise_beta_mean))stop("noise_beta_mean should be a matrix or nothing")
    if(is.matrix(scale_beta_mean))stop("scale_beta_mean should be a matrix or nothing")
    if(is.matrix(range_beta_mean))stop("range_beta_mean should be a matrix or nothing")
    
    if(is.matrix(noise_beta_precision))stop("noise_beta_precision should be a matrix or nothing")
    if(is.matrix(scale_beta_precision))stop("scale_beta_precision should be a matrix or nothing")
    if(is.matrix(range_beta_precision))stop("range_beta_precision should be a matrix or nothing")
    #length of observations
    if(
      !all(unique(c(
        length(observed_field),
        nrow(observed_locs),
        nrow(X),
        nrow(scale_X),
        nrow(noise_X),
        nrow(range_X),
        length(KL$idx)
        )) %in% c(0, length(observed_field))
      )) stop(
        paste("Lengths are not matching : observed_field has", length(observed_field), "observations,",
              "observed_locs has", nrow(observed_locs), "rows,", 
              "X has", nrow(X), "rows,", 
              "scale_X has", nrow(scale_X), "rows (can only be either 0 or the length of the observations),", 
              "noise_X has", nrow(noise_X), "rows (can only be either 0 or the length of the observations),", 
              "range_X has", nrow(range_X), "rows (can only be either 0 or the length of the observations),", 
              "KL has", length(KL$idx), "locations (can only be either 0 or the length of the observations)"
        )
      )
    
    # smoothness
    if (nu != 1.5)stop("only nu = 1.5")
    #stationary range is provided field or covariates
    if((length(grep("nonstationary", covfun))==0) & (!is.null(range_X)))stop(paste("Stationary covariance function", covfun, "should have no covariates (range_X)"))

  
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
    
    ##############
    # covariates #
    ##############
    

    covariates = list()
    
    covariates$X = process_covariates(X, observed_locs, vecchia_approx)  
    
    covariates$range_X = process_covariates(range_X, observed_locs, vecchia_approx, KL, range_KL)
    if(!identical(covariates$range_X$which_locs, seq(ncol(covariates$range_X$X_locs))))stop("The covariates range_X cannot vary within one spatial location of observed_locs")
    
    covariates$scale_X = process_covariates(scale_X, observed_locs, vecchia_approx, KL, scale_KL)
    if(!identical(covariates$scale_X$which_locs, seq(ncol(covariates$scale_X$X))))stop("The covariates scale_X cannot vary within one spatial location of observed_locs")
    
    covariates$noise_X = process_covariates(noise_X, observed_locs, vecchia_approx, KL, noise_KL)
    
    
    ############################################
    # More sanity checks for priors dimensions #
####    ############################################
####    
####      # priors for scale
####    if(
####      !all(unique(c(
####        ncol(noise_beta_precision),
####        nrow(noise_beta_precision),
####        nrow(noise_beta_mean),
####        ncol(covariates$noise_X$X)
####      )) %in% c(0, ncol(covariates$noise_X$X))
####      )) stop(
####        paste("Lengths are not matching : ", 
####              "noise_beta_precision has "  , ncol(noise_beta_precision)         , "columns",                                                      
####              "noise_beta_precision has "  , nrow(noise_beta_precision)         , "rows",                                                      
####              "noise_beta_mean has "       , nrow(noise_beta_mean)              , "rows",                                                 
####              "pre-processed noise_X has " , ncol(covariates$noise_X$X)         , "columns after adding an intercept (use model.matrix)"))     
####    # priors for noise
####    if(
####      all(unique(c(
####        ncol(scale_beta_precision),
####        nrow(scale_beta_precision),
####        nrow(scale_beta_mean),     
####        ncol(covariates$scale_X$X)       
####      )) %in% c(0, ncol(covariates$scale_X$X))
####      )) stop(
####        paste("Lengths are not matching : ", 
####              "scale_beta_precision has ", ncol(scale_beta_precision)          , "columns",                                                      
####              "scale_beta_precision has ", nrow(scale_beta_precision)          , "rows",                                                         
####              "scale_beta_mean has "     , nrow(scale_beta_mean)               , "rows",                                                    
####              "pre-processed scale_X has " , ncol(covariates$scale_X$X), "columns after adding an intercept (use model.matrix)"))         
####      # priors for range
####    if(
####      all(unique(c(
####        ncol(range_beta_precision)/(1+2*length(grep("anisotropic", covfun))),
####        nrow(range_beta_precision)/(1+2*length(grep("anisotropic", covfun))),
####        nrow(range_beta_mean),       
####        ncol(covariates$range_X$X)             
####      )) %in% c(0, ncol(covariates$range_X$X))
####      )) stop(
####        paste("Lengths are not matching : ", 
####              "range_beta_precision has ", ncol(range_beta_precision)           , "columns (should be equal to nrow(range_beta_mean) or 3*nrow(range_beta_mean) in case of nonstationary anisotropic function)",                                                      
####              "range_beta_precision has ", nrow(range_beta_precision)           , "rows    (should be equal to nrow(range_beta_mean) or 3*nrow(range_beta_mean) in case of nonstationary anisotropic function)",                                                      
####              "range_beta_mean has "     , nrow(range_beta_mean)                , "rows",                                                    
####              "pre-processed range_X has " , ncol(covariates$range_X$X), "columns after adding an intercept (use model.matrix)"))         
####      # priors for range mean wrt anisotropy
####    if(!is.null(range_beta_mean))
####    {
####      if(!ncol(range_beta_mean) %in%(c(1, 3)))stop("range_beta_mean must have 1 column (when isotropic function) or 3 columns (anisotropic function)")
####      if((ncol(range_beta_mean)!=3)&(length(grep("anisotropic", covfun))==1))stop("range_beta_mean must have 1 column (when isotropic function) or 3 columns (anisotropic function)")
####      if((ncol(range_beta_mean)!=1)&(length(grep("anisotropic", covfun))==0))stop("range_beta_mean must have 1 column (when isotropic function) or 3 columns (anisotropic function)")
####    }
    
    #################################
    # Info about hierarchical model #
    #################################
    
    hierarchical_model = list()
    hierarchical_model$covfun = covfun
    hierarchical_model$nu = nu
    hierarchical_model$beta_priors = list()
    hierarchical_model$KL = KL
    if(is.null(hierarchical_model$KL)) hierarchical_model$KL = list("n_KL" = 0)
    hierarchical_model$noise_KL = noise_KL
    hierarchical_model$scale_KL = scale_KL
    hierarchical_model$range_KL = range_KL
    
    if(is.null(noise_log_scale_prior)&noise_KL)
      {
      message("noise_log_scale_prior was automatically set to an uniform on (-8, 4)")
      noise_log_scale_prior = c(-8, 4)
      hierarchical_model$noise_log_scale_prior = matrix(noise_log_scale_prior)
      }
    if(is.null(scale_log_scale_prior)&scale_KL)
      {
      message("scale_log_scale_prior was automatically set to an uniform on (-8, 4)")
      scale_log_scale_prior = c(-8, 4)
      hierarchical_model$scale_log_scale_prior = matrix(scale_log_scale_prior)
      }
    if(is.null(range_log_scale_prior)&range_KL)
      {
      message("range_log_scale_prior was automatically set to an uniform on (-8, 4)")
      range_log_scale_prior = c(-8, 4)
      hierarchical_model$range_log_scale_prior = as.matrix(range_log_scale_prior) %x%rep(1, 1+2*length(grep("anisotropic", covfun)))
      }
    
    # OLS to get residual variance to make a guess 
    naive_ols =  lm(observed_field~covariates$X$X-1)
    lm_fit = as.vector(covariates$X$X%*%matrix(naive_ols$coefficients, ncol = 1))
    lm_residuals = as.vector(observed_field- lm_fit)
    
    hierarchical_model$beta_priors$noise_beta_mean = noise_beta_mean
    hierarchical_model$beta_priors$scale_beta_mean = scale_beta_mean
    hierarchical_model$beta_priors$range_beta_mean = range_beta_mean
    # Default mean prior computed from a reasonable case. 
    # The intercept is set to reasonable value and the rest is set to 0
    if(is.null(noise_beta_mean)) hierarchical_model$beta_priors$noise_beta_mean = matrix(c(log(var(lm_residuals)) - log(2),                     rep(0, covariates$noise_X$n_regressors-1)), ncol=1)
    if(is.null(scale_beta_mean)) hierarchical_model$beta_priors$scale_beta_mean = matrix(c(log(var(lm_residuals)) - log(2),                     rep(0, covariates$scale_X$n_regressors-1)), ncol=1)
    if(is.null(range_beta_mean)) hierarchical_model$beta_priors$range_beta_mean = matrix(0, covariates$range_X$n_regressors, 1 + 2*length(grep("anisotropic", covfun)))
    if(is.null(range_beta_mean)) hierarchical_model$beta_priors$range_beta_mean[1,1] = c(log(max(dist(locs[seq(1000), seq(2)])))-log(50))
    hierarchical_model$beta_priors$noise_beta_precision = noise_beta_precision
    hierarchical_model$beta_priors$scale_beta_precision = scale_beta_precision
    hierarchical_model$beta_priors$range_beta_precision = range_beta_precision
    if(is.null(noise_beta_precision)) hierarchical_model$beta_priors$noise_beta_precision =diag(.0001, covariates$noise_X$n_regressors, covariates$noise_X$n_regressors)
    if(is.null(scale_beta_precision)) hierarchical_model$beta_priors$scale_beta_precision =diag(.0001, covariates$scale_X$n_regressors, covariates$scale_X$n_regressors)
    if(is.null(range_beta_precision)) 
      {
      hierarchical_model$beta_priors$range_beta_precision = diag(.0001, covariates$range_X$n_regressors, covariates$range_X$n_regressors)
      hierarchical_model$beta_priors$range_beta_precision = hierarchical_model$beta_priors$range_beta_precision %x% diag(1, (1+2*length(grep("anisotropic", covfun))))
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
      if(covfun %in% c("exponential_isotropic"                     , "matern_isotropic"                     )) states[[i]]$params$range_beta = matrix(  sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), 1         )                                                                                                                    , nrow = 1)
      if(covfun %in% c("exponential_sphere"                        , "matern_sphere"                        )) states[[i]]$params$range_beta = matrix(  sample(log(max(fields::rdist.earth(locs[sample(seq(vecchia_approx$n_locs), 1000),           ]))/6000)-log(seq(50, 500, 1)), 1         )                                                                                                            , nrow = 1)
      if(covfun %in% c("nonstationary_exponential_isotropic"       , "nonstationary_matern_isotropic"       )) states[[i]]$params$range_beta = matrix(  sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), 1         )                                                                                                                    , nrow = 1)
      if(covfun %in% c("nonstationary_exponential_isotropic_sphere", "nonstationary_matern_isotropic_sphere")) states[[i]]$params$range_beta = matrix(  sample(log(max(fields::rdist.earth(locs[sample(seq(vecchia_approx$n_locs), 1000),           ]))/6000)-log(seq(50, 500, 1)), 1         )                                                                                                            , nrow = 1)
      if(covfun %in% c("exponential_anisotropic2D"                 , "matern_anisotropic2D"                 )) states[[i]]$params$range_beta = matrix(c(sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), ncol(locs)),                       rep(0, ncol(locs)*(ncol(locs)+1)/2-ncol(locs)))                                             , nrow = 1)
      if(covfun == "nonstationary_exponential_anisotropic")                                                    states[[i]]$params$range_beta = matrix(c(sample(log(max(       dist(locs[sample(seq(vecchia_approx$n_locs)        , 1000),           ]))     )-log(seq(50, 500, 1)), ncol(locs)),                       rep(0, ncol(locs)*(ncol(locs)+1)/2-ncol(locs)))                                             , nrow = 1)
      if(covfun == "nonstationary_exponential_anisotropic_sphere")                                             states[[i]]$params$range_beta = matrix(c(sample(log(max(fields::rdist.earth(locs[sample(seq(vecchia_approx$n_locs), 1000), c(1, 2)   ]))/6000)-log(seq(50, 500, 1)), 2         ), rep(1, ncol(locs)-2), rep(0, ncol(locs)*(ncol(locs)+1)/2-ncol(locs)))                                     , nrow = 1)
      # initiate null regressors for other covariates                                                                                                                                          
      states[[i]]$params$range_beta = rbind(states[[i]]$params$range_beta, matrix(0, ncol(covariates$range_X$X_locs)+hierarchical_model$KL$n_KL*range_KL-1, length(states[[i]]$params$range_beta)))
      if(range_KL)row.names(states[[i]]$params$range_beta) = c(colnames(covariates$range_X$X_locs), paste("KL", seq(hierarchical_model$KL$n_KL), sep = "_"))
      if(!range_KL)row.names(states[[i]]$params$range_beta) = c(colnames(covariates$range_X$X_locs))
      if(range_KL)states[[i]]$params$range_log_scale = runif(1 + 2*length(grep("aniso", covfun)), hierarchical_model$range_log_scale_prior[1,], hierarchical_model$range_log_scale_prior[2,])
      states[[i]]$momenta$range_beta_ancillary = matrix(rnorm(length(states[[i]]$params$range_beta)), nrow(states[[i]]$params$range_beta))
      states[[i]]$momenta$range_beta_sufficient = matrix(rnorm(length(states[[i]]$params$range_beta)), nrow(states[[i]]$params$range_beta))
        
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
      states[[i]]$transition_kernels$range_log_scale_sufficient = -4
      states[[i]]$transition_kernels$range_log_scale_ancillary =  -4
      states[[i]]$transition_kernels$range_beta_sufficient = -4
      states[[i]]$transition_kernels$range_beta_ancillary  = -4
      # scale
      states[[i]]$transition_kernels$scale_beta_sufficient_mala = -4
      states[[i]]$transition_kernels$scale_beta_ancillary_mala  = -4
      states[[i]]$transition_kernels$scale_log_scale_sufficient = -4
      states[[i]]$transition_kernels$scale_log_scale_ancillary =  -4
      # noise variance
      states[[i]]$transition_kernels$noise_beta_mala = -2
      states[[i]]$transition_kernels$noise_log_scale = -2
      
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
      states[[i]]$params$noise_beta    = matrix(rep(0, ncol(covariates$noise_X$X)+noise_KL*hierarchical_model$KL$n_KL), ncol = 1) #random starting values
      states[[i]]$params$noise_beta[1] = log(var(states[[i]]$sparse_chol_and_stuff$lm_residuals)) - log(2) +rnorm(1, 0, .5) # setting sensible value for the intercept
      if(!noise_KL) row.names(states[[i]]$params$noise_beta) = colnames(covariates$noise_X$X)
      if(noise_KL ) row.names(states[[i]]$params$noise_beta) = c(colnames(covariates$noise_X$X), paste("KL", seq(hierarchical_model$KL$n_KL), sep = "_"))
      states[[i]]$momenta$noise_beta = rnorm(ncol(covariates$noise_X$X)+noise_KL*hierarchical_model$KL$n_KL)
      # effective variance field, shall be used in density computations
      states[[i]]$sparse_chol_and_stuff$noise = Bidart::variance_field(beta = states[[i]]$params$noise_beta, KL = KL, use_KL = noise_KL, X = covariates$noise_X$X)
      # log scale if latent field
      if(noise_KL)states[[i]]$params$noise_log_scale = runif(1, hierarchical_model$noise_log_scale_prior[1], hierarchical_model$noise_log_scale_prior[2])
      
      #########
      # Scale #
      #########
      # beta is just an intercept in stationary case
      states[[i]]$params$scale_beta    = matrix(0, ncol(covariates$scale_X$X_locs) + scale_KL * hierarchical_model$KL$n_KL, ncol = 1) #random starting values
      states[[i]]$momenta$scale_beta_ancillary =  rnorm(ncol(covariates$scale_X$X_locs) + scale_KL * hierarchical_model$KL$n_KL)
      states[[i]]$momenta$scale_beta_sufficient = rnorm(ncol(covariates$scale_X$X_locs) + scale_KL * hierarchical_model$KL$n_KL)
      states[[i]]$params$scale_beta[1] = log(var(states[[i]]$sparse_chol_and_stuff$lm_residuals)) - log(2) + rnorm(1, 0, .5) # setting sensible value for the intercept
      if(!scale_KL)row.names(states[[i]]$params$scale_beta) = colnames(covariates$scale_X$X_locs)
      if(scale_KL)row.names(states[[i]]$params$scale_beta) = c(colnames(covariates$scale_X$X_locs),  paste("KL", seq(hierarchical_model$KL$n_KL), sep = "_"))
      if(scale_KL)states[[i]]$params$scale_log_scale = runif(1, hierarchical_model$scale_log_scale_prior[1], hierarchical_model$scale_log_scale_prior[2])
      # effective variance field, shall be used in density computations
      states[[i]]$sparse_chol_and_stuff$scale = Bidart::variance_field(beta = states[[i]]$params$scale_beta, KL = KL, use_KL = scale_KL, X = covariates$scale_X$X_locs, locs_idx = vecchia_approx$hctam_scol_1)
      #######
      
      ####################
      # NNGP sparse chol #
      ####################
      states[[i]]$sparse_chol_and_stuff$compressed_sparse_chol_and_grad = 
        Bidart::compute_sparse_chol(covfun_name = hierarchical_model$covfun, 
                                    range_X = covariates$range_X$X_locs, 
                                    range_beta = states[[i]]$params$range_beta, 
                                    KL = hierarchical_model$KL, use_KL = hierarchical_model$range_KL, 
                                    NNarray = vecchia_approx$NNarray, locs_idx = vecchia_approx$hctam_scol_1, 
                                    locs = locs, nu = nu)
      states[[i]]$sparse_chol_and_stuff$sparse_chol = Matrix::sparseMatrix(x =  states[[i]]$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA], i = vecchia_approx$sparse_chol_row_idx, j = vecchia_approx$sparse_chol_column_idx, triangular = T)
      states[[i]]$sparse_chol_and_stuff$precision_diag = as.vector((states[[i]]$sparse_chol_and_stuff$compressed_sparse_chol_and_grad[[1]][vecchia_approx$NNarray_non_NA]^2)%*%Matrix::sparseMatrix(i = seq(length(vecchia_approx$sparse_chol_column_idx)), j = vecchia_approx$sparse_chol_column_idx, x = rep(1, length(vecchia_approx$sparse_chol_row_idx))))
      ################
      # Latent field #
      ################
      states[[i]]$params$field = sqrt(states[[i]]$sparse_chol_and_stuff$scale) * as.vector(Matrix::solve(states[[i]]$sparse_chol_and_stuff$sparse_chol, rnorm(vecchia_approx$n_locs)))
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