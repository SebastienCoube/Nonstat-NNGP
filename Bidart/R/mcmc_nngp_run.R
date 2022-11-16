#' @export
mcmc_nngp_run_nonstationary = function(mcmc_nngp_list, 
                         burn_in = .5, seed = 1, # MCMC parameters
                         n_cores = NULL, thinning = .1,
                         plot_diags = T, 
                         big_range = F, 
                         plot_PSRF_fields = F, debug_outfile = NULL
                         )
{
  t_start = Sys.time()
  ######################################
  # Default parameters following model # 
  ######################################
  # number of cores
  n_cores = min(c(n_cores, length(mcmc_nngp_list$states), parallel::detectCores()-1))
  gc()
  #################
  # Sanity checks #
  #################
  # parallelization  
  if(!is.numeric(n_cores))stop("n_cores must be a positive round number")
  if((floor(n_cores)!=n_cores) |  n_cores<1)stop("n_cores must be a positive round number")
  # iterations and thinning
  if(!is.numeric(thinning))stop("thinning is a proportion and must be between 0 and 1")
  if((thinning<0)|(thinning>1))stop("thinning is a proportion and must be between 0 and 1")
  cl = parallel::makeCluster(n_cores, outfile = debug_outfile)
  mcmc_nngp_list_ = mcmc_nngp_list[-match(c("records", "states"), names(mcmc_nngp_list))]
  parallel::clusterExport(cl = cl, varlist = c("mcmc_nngp_list_", "mcmc_nngp_update_Gaussian", "seed"), envir = environment())
  rm(mcmc_nngp_list_);gc()
  #parallel::clusterEvalQ(cl = cl, expr = "library(Matrix);library(GpGp);library(Bidart);library(expm);library(MfUSampler)")
  #parallel::clusterEvalQ(cl = cl, expr = "library(irlba, lib.loc = '/home/user/s/scoube/R_packages/');library(Matrix, lib.loc = '/home/user/s/scoube/R_packages/');library(GpGp, lib.loc = '/home/user/s/scoube/R_packages/');library(Bidart, lib.loc = '/home/user/s/scoube/R_packages/');library(expm, lib.loc = '/home/user/s/scoube/R_packages/');library(MfUSampler, lib.loc = '/home/user/s/scoube/R_packages/')")
  iter_start = mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]
  # mcmc_sampling
  samples = 
  parallel::parLapply(cl = cl, X = seq(length(mcmc_nngp_list$states)), 
    fun 
     = function(i)Bidart::mcmc_nngp_update_Gaussian(
      data = mcmc_nngp_list_$data, 
      hierarchical_model = mcmc_nngp_list_$hierarchical_model, 
      vecchia_approx = mcmc_nngp_list_$vecchia_approx, 
      state = mcmc_nngp_list$states[[i]], 
      n_iterations_update = 100, thinning = thinning, 
      big_range = big_range, 
      iter_start = iter_start, 
      seed = seed + iter_start + i
      )
    )
  gc()
  # saving mcmc samples directly in parent environment
  for(i in seq(length(samples)))
  {
    # update states
    mcmc_nngp_list$states[[i]] = samples[[i]]$state
    # update records
    for(j in names(samples[[1]]$params_records))
    {
      mcmc_nngp_list$records[[i]][[j]] = abind::abind(mcmc_nngp_list$records[[i]][[j]], samples[[i]]$params_records[[j]])
    }
        
  }
  # update iterations
  mcmc_nngp_list$iterations$checkpoints = rbind(mcmc_nngp_list$iterations$checkpoints, c(iter_start + 100,  mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 2] + as.numeric(Sys.time()-t_start, unit = "mins")))
  print(mcmc_nngp_list$iterations$checkpoints)
  mcmc_nngp_list$iterations$thinning = c(mcmc_nngp_list$iterations$thinning, iter_start + which(seq(100)*thinning == floor(seq(100)*thinning)))
  # plot diagnostics
  if(plot_diags) diagnostic_plots(mcmc_nngp_list, plot_PSRF_fields = plot_PSRF_fields, burn_in = burn_in)
  parallel::stopCluster(cl = cl)
  return(mcmc_nngp_list)
}
