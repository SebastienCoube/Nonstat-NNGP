#' @export
mcmc_nngp_run_nonstationary = function(mcmc_nngp_list, 
                         burn_in = .5, seed = 1, # MCMC parameters
                         n_cores = NULL, thinning = .1, n_iterations_update = 300, #run parameters
                         n_cycles = 5, plot_diags = T, 
                         field_n_chromatic = 3, field_n_mala = 1, 
                         plot_PSRF_fields = F, debug_outfile = "debug_mcmc_nngp.txt"
                         )
{
  t_start = Sys.time()
  ######################################
  # Default parameters following model # 
  ######################################
  # number of cores
  if(is.null(n_cores)) n_cores = min(length(mcmc_nngp_list$states), parallel::detectCores()-1)
  #################
  # Sanity checks #
  #################
  # parallelization  
  if(!is.numeric(n_cores))stop("n_cores must be a positive round number")
  if((floor(n_cores)!=n_cores) |  n_cores<1)stop("n_cores must be a positive round number")
  # iterations and thinning
  if(!is.numeric(n_cycles))stop("n_cycles must be a positive round number")
  if((floor(n_cycles)!=n_cycles) |  n_cycles<1)stop("n_cycles must be a positive round number")
  if(!is.numeric(n_iterations_update))stop("n_iterations_update must be a positive round number")
  if((floor(n_iterations_update)!=n_iterations_update) |  n_iterations_update<1)stop("n_iterations_update must be a positive round number")
  if(n_iterations_update<50)stop("n_iterations_update must be bigger than or equal to 50")
  if(!is.numeric(thinning))stop("thinning is a proportion and must be between 0 and 1")
  if((thinning<0)|(thinning>1))stop("thinning is a proportion and must be between 0 and 1")
  # field
  if(!is.numeric(field_n_chromatic))stop("field_n_chromatic must be a positive round number")
  if((floor(field_n_chromatic)!=field_n_chromatic)|  field_n_chromatic<0)stop("field_n_chromatic must be a positive round number")
  if(!is.numeric(field_n_mala))stop("field_n_mala must be a positive round number")
  if((floor(field_n_mala)!=field_n_mala)|  field_n_mala<0)stop("field_n_mala must be a positive round number")
  if((field_n_chromatic==0) & (field_n_mala==0)) stop("Either field_n_chromatic or field_n_mala must be different from 0")
  gc()
  cl = parallel::makeCluster(n_cores, outfile = debug_outfile)
  mcmc_nngp_list_ = mcmc_nngp_list[-match(c("records", "states"), names(mcmc_nngp_list))]
  parallel::clusterExport(cl = cl, varlist = c("mcmc_nngp_list_", "mcmc_nngp_update_Gaussian", "seed"), envir = environment())
  rm(mcmc_nngp_list_);gc()
  #parallel::clusterEvalQ(cl = cl, expr = "library(Matrix);library(GpGp);library(Bidart);library(expm);library(MfUSampler)")
  #parallel::clusterEvalQ(cl = cl, expr = "library(Matrix, lib.loc = '/home/user/s/scoube/R_packages/');library(GpGp, lib.loc = '/home/user/s/scoube/R_packages/');library(Bidart, lib.loc = '/home/user/s/scoube/R_packages/');library(expm, lib.loc = '/home/user/s/scoube/R_packages/');library(MfUSampler, lib.loc = '/home/user/s/scoube/R_packages/')")
  i_cycle = 1
  while(i_cycle <= n_cycles)
  {
    print(paste("cycle =", i_cycle))
    iter_start = mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 1]
    # mcmc_sampling
    samples = 
    parallel::parLapply(cl = cl, X = seq(length(mcmc_nngp_list$states)), 
      fun 
       = function(i)mcmc_nngp_update_Gaussian(
        data = mcmc_nngp_list_$data, 
        hierarchical_model = mcmc_nngp_list_$hierarchical_model, 
        vecchia_approx = mcmc_nngp_list_$vecchia_approx, 
        state = mcmc_nngp_list$states[[i]], 
        n_iterations_update = n_iterations_update, thinning = thinning, 
        field_n_chromatic = field_n_chromatic, field_n_mala = field_n_mala, 
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
    mcmc_nngp_list$iterations$checkpoints = rbind(mcmc_nngp_list$iterations$checkpoints, c(iter_start + n_iterations_update,  mcmc_nngp_list$iterations$checkpoints[nrow(mcmc_nngp_list$iterations$checkpoints), 2] + as.numeric(Sys.time()-t_start, unit = "mins")))
    print(mcmc_nngp_list$iterations$checkpoints)
    t_start = Sys.time()
    mcmc_nngp_list$iterations$thinning = c(mcmc_nngp_list$iterations$thinning, iter_start + which(seq(n_iterations_update)*thinning == floor(seq(n_iterations_update)*thinning)))
    # plot diagnostics
    diagnostic_plots(mcmc_nngp_list, plot_PSRF_fields = plot_PSRF_fields, burn_in = burn_in)
    i_cycle = i_cycle+1
  }
  #parallel::stopCluster(cl = cl)
  return(mcmc_nngp_list)
}
