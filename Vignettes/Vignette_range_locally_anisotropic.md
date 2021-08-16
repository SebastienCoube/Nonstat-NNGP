This toy example shows how to model *locally anisotropic range*. It
shows :

-   how to initialize the model
-   how to diagnose convergence
-   how to compare the models (with a stationary model) using DIC
-   how to do predictions

Generate synthetic data
-----------------------

    set.seed(1)
    # Generate locations
    locs = cbind(5*runif(6000), 5*runif(6000))
    #Generate latent fields for the parameters
    latent_field_range =  
      cbind(GpGp::fast_Gp_sim(c(.5, .5, 1, 0), locs = locs), 
            GpGp::fast_Gp_sim(c(.5, .5, 1, 0), locs = locs), 
            GpGp::fast_Gp_sim(c(.5, .5, 1, 0), locs = locs))
    # observing with duplicates
    n_obs = 2 * nrow(locs)
    observation_idx = c(sample(nrow(locs), nrow(locs) , replace =F), sample(nrow(locs), n_obs - nrow(locs) , replace =T))
    observed_locs = locs[observation_idx,]
    # covariates
    X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5)))
    colnames(X) = c("binomial")
    # regression coeffs
    beta = c(100, 5)
    beta_range = matrix(c(-1.6,.0, .0, -1.6,.0, .0, 0 ,.0, .0), ncol = 3)
    # Use NNGP to generate data
    NNarray = GpGp::find_ordered_nn(locs, 10) # Nearest Neighbor Array
    sparse_chol = Bidart::compute_sparse_chol(covfun_name = "nonstationary_exponential_anisotropic", range_beta = beta_range, locs = locs, NNarray = NNarray, range_field =  latent_field_range, range_X = as.matrix(cbind(1, locs)))
    # get latent NNGP field and observed Gaussian field
    latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)
    observed_field =  as.vector(latent_field[observation_idx] +rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta)

Let’s plot the latent field for the NNGP’s range…

    Bidart::plot_ellipses(locs[seq(200),], latent_field_range[seq(200),], shrink = .01, main = "range_ellipses")

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-2-1.png)

Let’s plot the latent field…

    Bidart::plot_pointillist_painting(locs, latent_field)

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-3-1.png)

Nonstationary model
-------------------

Initialize

    mcmc_nngp_list_range = Bidart::mcmc_nngp_initialize_nonstationary (
      observed_locs = observed_locs, #spatial locations
      observed_field = c(observed_field), # Response variable
      X = X, # Covariates for the observed field
      m = 5, #number of Nearest Neighbors
      reordering = c("maxmin"), #Reordering
      covfun = "nonstationary_exponential_anisotropic", response_model = "Gaussian", # covariance model and response model
      range_range = .5, 
      log_NNGP_matern_covfun = "matern_isotropic", # covariance function for the hyperpriors
      log_NNGP_matern_smoothness = 1, # covariance function for the hyperpriors
      n_chains = 3,  # number of MCMC chains
      seed = 10
    )

    ## Setup done, 3.14248633384705 s elapsed

Run

    for(i in seq(60))
    {
      mcmc_nngp_list_range = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list_range, n_iterations_update = 100, n_cycles = 1)
    }

Check convergence. (Normally diagnostics pop automatically or are all
generated using diagnostic\_plots() and look great in Rstudio, but in
order to get decent pictures in the md lower functions are called)

    # normally one should use diagnostic_plots()  , but in order to get good pictures lower functions are used
    name = "range_log_scale"
    Bidart::plot_PSRF(PSRF = Bidart::grb_diags_field(record_arrays = lapply(mcmc_nngp_list_range$records, function(x)x[[name]]), iterations = mcmc_nngp_list_range$iterations$thinning), varname = name, individual_varnames = row.names(mcmc_nngp_list_range$states$chain_1$params[[name]]))

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    # normally one should use diagnostic_plots()  , but in order to get good pictures lower functions are used
    name = "range_beta"
    Bidart::plot_PSRF(PSRF = Bidart::grb_diags_field(record_arrays = lapply(mcmc_nngp_list_range$records, function(x)x[[name]]), iterations = mcmc_nngp_list_range$iterations$thinning), varname = name, individual_varnames = row.names(mcmc_nngp_list_range$states$chain_1$params[[name]]))

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-7-1.png)

    # normally one should use diagnostic_plots()  , but in order to get good pictures lower functions are used
    name = "scale_beta"
    Bidart::plot_PSRF(PSRF = Bidart::grb_diags_field(record_arrays = lapply(mcmc_nngp_list_range$records, function(x)x[[name]]), iterations = mcmc_nngp_list_range$iterations$thinning), varname = name, individual_varnames = row.names(mcmc_nngp_list_range$states$chain_1$params[[name]]))

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-8-1.png)

    # normally one should use diagnostic_plots()  , but in order to get good pictures lower functions are used
    name = "noise_beta"
    Bidart::plot_PSRF(PSRF = Bidart::grb_diags_field(record_arrays = lapply(mcmc_nngp_list_range$records, function(x)x[[name]]), iterations = mcmc_nngp_list_range$iterations$thinning), varname = name, individual_varnames = row.names(mcmc_nngp_list_range$states$chain_1$params[[name]]))

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-9-1.png)

    # normally one should use diagnostic_plots()  , but in order to get good pictures lower functions are used
    name = "beta"
    Bidart::plot_PSRF(PSRF = Bidart::grb_diags_field(record_arrays = lapply(mcmc_nngp_list_range$records, function(x)x[[name]]), iterations = mcmc_nngp_list_range$iterations$thinning), varname = name, individual_varnames = row.names(mcmc_nngp_list_range$states$chain_1$params[[name]]))

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-10-1.png)

Stationary model
----------------

Initialize

    mcmc_nngp_list_stat = Bidart::mcmc_nngp_initialize_nonstationary (
      observed_locs = observed_locs, #spatial locations
      observed_field = c(observed_field), # Response variable
      X = X, # Covariates for the observed field
      m = 5, #number of Nearest Neighbors
      reordering = c("maxmin"), #Reordering
      covfun = "exponential_isotropic", response_model = "Gaussian", # covariance model and response model
      n_chains = 3,  # number of MCMC chains
      seed = 10
    )

    ## Setup done, 0.616835355758667 s elapsed

Run

    for(i in seq(40))
    {
      mcmc_nngp_list_stat = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list_stat, n_iterations_update = 100, n_cycles = 1)
    }

Convergence is checked using diagnostic\_plots()

Model comparison
----------------

    print(Bidart::DIC(mcmc_nngp_list_range))

    ## [1] 36645.31

    print(Bidart::DIC(mcmc_nngp_list_stat))

    ## [1] 36709.73

Prediction
----------

    predicted_locs = as.matrix(expand.grid(seq(0, 5, .1), seq(0, 5, .1)))
    pred = Bidart::predict_latent_field(mcmc_nngp_list = mcmc_nngp_list_range, predicted_locs = predicted_locs)

    par(mfrow = c(1, 1))
    Bidart::plot_pointillist_painting(pred$predicted_locs, pred$summaries$field[1,,], cex = 1.8, main = "predicted mean")

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-15-1.png)

    Bidart::plot_pointillist_painting(pred$predicted_locs, pred$summaries$field[5,,], cex = 1.8, main = "standard deviation of the prediction")

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-15-2.png)

    Bidart::plot_ellipses(pred$predicted_locs, log_range = pred$summaries$range_field[1,,], shrink = .01, main = "range ellipses (a posteriori mean of log-range parameters)")

![](Vignette_range_locally_anisotropic_files/figure-markdown_strict/unnamed-chunk-15-3.png)
