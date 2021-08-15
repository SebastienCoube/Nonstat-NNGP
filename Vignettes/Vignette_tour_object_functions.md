This vignette gives a tour of the objects and functions of the package
Bidart, reviewing the arguments of the functions and the contents of the
objects.

Generate synthetic data
-----------------------

First, generate some synthetic data (heteroscedasticity of the latent
field and the noise too)

    set.seed(1)
    # Generate locations
    locs = cbind(5*runif(3000), 5*runif(3000))
    #Generate latent fields for the parameters
    latent_field_noise =   GpGp::fast_Gp_sim(c(.5, 1, 1, 0), locs = locs)
    latent_field_scale =   GpGp::fast_Gp_sim(c(.5, 1, 1, 0), locs = locs)
    # set number of observations
    n_obs = 2 * nrow(locs)
    # observing with duplicates
    observation_idx = c(sample(nrow(locs), nrow(locs) , replace =F), sample(nrow(locs), n_obs - nrow(locs) , replace =T))
    observed_locs = locs[observation_idx,]
    # covariates
    X = as.data.frame(cbind(rbinom(n = n_obs, size = 1, prob = .5)))
    colnames(X) = c("binomial")
    X_noise = X
    # regression coeffs
    beta = c(100, 5)
    beta_noise = c(-.5, 1)
    # get logarithm of the parameters
    log_noise_variance = as.matrix(cbind(1, X_noise)) %*% beta_noise + latent_field_noise[observation_idx]
    log_scale = as.matrix(latent_field_scale[observation_idx])
    # Use NNGP to generate data
    NNarray = GpGp::find_ordered_nn(locs, 10) # Nearest Neighbor Array
    sparse_chol = Bidart::compute_sparse_chol(covfun_name = "exponential_isotropic", range_beta = log(.1), locs = locs, NNarray = NNarray) # stationary covariance
    # get latent NNGP field and observed Gaussian field
    latent_field = GpGp::fast_Gp_sim_Linv(Linv = sparse_chol[[1]], NNarray = NNarray)
    observed_field = as.vector(exp(.5*log_scale) * latent_field[observation_idx] +exp(.5*log_noise_variance)*rnorm(n_obs)+ cbind(1, as.matrix(X))%*%beta)

Initialize the model
--------------------

### Arguments

The model object is a *list* created with the function
**Bidart::mcmc\_nngp\_initialize\_nonstationary** In order to initialize
the model, you need :

-   **observed\_locs**, a *matrix* of spatial locations with *2* columns
    (yet).
-   **observed\_field**, a *vector* of observations done at
    observed\_locs. You need to have nrow(observed\_locs) ==
    length(observed field)

You can further indicate :

-   **covfun**, a *string* for the covariance function name. Takes the
    value “exponential\_isotropic” by default (can be :
    “exponential\_isotropic”, “exponential\_sphere”,
    “exponential\_spacetime”, “exponential\_spheretime”,
    “exponential\_anisotropic”, “nonstationary\_exponential\_isotropic”,
    “nonstationary\_exponential\_isotropic\_sphere”,
    “nonstationary\_exponential\_anisotropic”,
    “nonstationary\_exponential\_anisotropic\_sphere”).
-   **m**, an *integer* for the number of parents (by default 5).
-   **reordering**, a *string* indicating the reordering. Takes the
    value “maxmin” by default (can be : “maxmin”, “random”, “coord”,
    “dist\_to\_point”, “middleout”).
-   **response\_model**, a *string* indicating the data model (for now
    and by default, “Gaussian”).
-   **n\_chains**, an *integer* indicating the number of MCMC chains (by
    default 3).
-   **seed**, an *integer* giving the seed (by default 1).
-   **X**, a *data.frame* of covariates to explain the observed
    variables (intercept is added automatically).
-   **noise\_X**, a *data.frame* of covariates to explain the
    heteroscedasticity of the noise (intercept is added automatically).
-   **noise\_range**, a *positive number* to indicate the range of the
    log-NNGP prior for the noise variance (initializes a log-NNGP latent
    field for the noise variance).
-   **scale\_X**, a *data.frame* of covariates to explain the
    heteroscedasticity of the latent field. The covariates must be
    constant within a spatial location (intercept is added
    automatically).
-   **scale\_range**, a *positive number* to indicate the range of the
    log-NNGP prior for the latent field variance (initializes a log-NNGP
    latent field for the latent field variance).

If you indicated a nonstationary covariance function
(“nonstationary\_exponential\_isotropic”,
“nonstationary\_exponential\_isotropic\_sphere”,
“nonstationary\_exponential\_anisotropic”,
“nonstationary\_exponential\_anisotropic\_sphere”), you need to precise
**either or both** :

-   **range\_X**, a *data.frame* of covariates to explain the range of
    the latent field. The covariates must be constant within a spatial
    location (intercept is added automatically).
-   **range\_range**, a *positive number* to indicate the range of the
    log-NNGP prior for the latent field variance (initializes a log-NNGP
    or matrix log-NNGP latent field for the latent field range).

If you filled either noise\_range, scale\_range, or range\_range, you
need to precise :

-   **log\_NNGP\_matern\_covfun**, a *string* indicating a covariance
    function of the Matérn family (can be : “matern\_isotropic”,
    “matern\_sphere”).
-   if you want, you can indicate **log\_NNGP\_matern\_smoothness**, a
    *positive number* for log\_NNGP\_matern\_covfun (by default 1).

Here we go, initialize a model with heteroscedastic latent field and
noise.

    mcmc_nngp_list_nonstat = Bidart::mcmc_nngp_initialize_nonstationary (
      observed_locs = observed_locs, #spatial locations
      observed_field = c(observed_field), # Response variable
      X = X, # Covariates for the observed field
      m = 10, #number of Nearest Neighbors
      reordering = c("maxmin"), #Reordering
      covfun = "exponential_isotropic", response_model = "Gaussian", # covariance model and response model
      noise_X = X_noise, # covariates for the noise
      noise_range = .5, # range for latent field of parameters, if NULL no latent field
      scale_X = NULL, # covariates for the scale
      scale_range = .5, # range for latent field of parameters, if NULL no latent field
      log_NNGP_matern_covfun = "matern_isotropic", # covariance function for the hyperpriors
      log_NNGP_matern_smoothness = 1, # covariance function for the hyperpriors
      n_chains = 3,  # number of MCMC chains
      seed = 10
    )

    ## Setup done, 1.52983450889587 s elapsed

    names(mcmc_nngp_list_nonstat)

    ## [1] "data"               "hierarchical_model" "vecchia_approx"    
    ## [4] "states"             "records"            "t_begin"           
    ## [7] "seed"               "iterations"

### What is there in the model ?

#### **data** is a *list* containing the information about the observations. It contains :

-   **observed\_locs**, the *matrix* of observed locations.
-   **observed field**, the *vector* of observations.
-   **locs**, a *matrix* of duplicate-less, re-orderd locations.
-   **covariates**, a *list* about the covariates of the latent field,
    the covariates of the scale, noise, range (even if it is just an
    intercept).

<!-- -->

    names(mcmc_nngp_list_nonstat$data)

    ## [1] "locs"           "observed_field" "observed_locs"  "covariates"

    names(mcmc_nngp_list_nonstat$data$covariates)

    ## [1] "X"       "range_X" "scale_X" "noise_X"

#### **hierarchical\_model** is a *list* containing the information about the hierarchical model and its hyperpriors. It contains :

-   **response\_model** is a *string* indicating the response model
    (only “Gaussian” yet).
-   **covfun** is a *string* indicating the covariance function.
-   **hyperprior\_covariance** is a *list* giving information about the
    log-GP and matrix log-GP priors.

<!-- -->

    names(mcmc_nngp_list_nonstat$hierarchical_model)

    ## [1] "response_model"        "covfun"                "hyperprior_covariance"

    mcmc_nngp_list_nonstat$hierarchical_model$response_model

    ## [1] "Gaussian"

    mcmc_nngp_list_nonstat$hierarchical_model$covfun

    ## [1] "exponential_isotropic"

    names(mcmc_nngp_list_nonstat$hierarchical_model$hyperprior_covariance)

    ## [1] "log_NNGP_matern_covfun" "scale_NNGP_prior"       "noise_NNGP_prior"

#### **vecchia\_approx** is a *list* containing the information about the NNGP graph. It contains among others :

-   **NNarray**, the *matrix of integers* indicating the parents in the
    DAG (cf documentation of GpGP).
-   **MRF\_adjacency\_mat**, a *sparse matrix* of adjacency for the
    moral graph.
-   **coloring**, a *vector of integers* indicating the colors on the
    moral graph.
-   **locs\_match\_matrix**, a *sparse matrix* matching un-duplicated
    locations of locs and observations from observed\_field possibly
    done at the same site. The rest is ancillary objects derived from
    these and used as shortcuts.

<!-- -->

    names(mcmc_nngp_list_nonstat$vecchia_approx)

    ##  [1] "reordering"             "n_locs"                 "n_obs"                 
    ##  [4] "locs_match"             "locs_match_matrix"      "hctam_scol"            
    ##  [7] "hctam_scol_1"           "obs_per_loc"            "NNarray"               
    ## [10] "NNarray_non_NA"         "sparse_chol_column_idx" "sparse_chol_row_idx"   
    ## [13] "MRF_adjacency_mat"      "coloring"               "duplicated_locs"

#### **states** is a *list* containing the information about the current state of each MCMC chain.

There are n\_chains elements in this list, each corresponding to a
chain.

-   **params**, is a *list* gathering the chain parameters.
-   **sparse\_chol\_and\_stuff**, is a *list* of ancillary objects
    derived from params.  
-   **momenta**, is a *list* of the HMC momenta.
-   **transition\_kernels**, is a *list* of Metropolis-Hastings or HMC
    parameters (such as proposal variance, leapfrog integration step…).

<!-- -->

    names(mcmc_nngp_list_nonstat$states$chain_1)

    ## [1] "params"                "sparse_chol_and_stuff" "momenta"              
    ## [4] "transition_kernels"

The rest is empty for now, we must run the model.

Run the model
-------------

### Arguments

In order to run the model, the function
**Bidart::mcmc\_nngp\_run\_nonstationary** takes the existing list,
appends MCMC samples, and returns a new list.

You need :

-   **mcmc\_nngp\_list**, the *list* containing the model.

You can add computational options :

-   **n\_cyles**, an *integer* indicating the number of MCMC cycles (by
    default 5)
-   **n\_iterations\_update**, an *integer* greater than 50 indicating
    the number of MCMC iteration per cycle (by default 300)
-   **n\_cores**, an *integer* indicating the number of cores (by
    default min(n\_chains, parallel::detectCores()-1))
-   **thinning**, a *number* between 0 and 1 indicating the proportion
    of kept iterations (by default 0.1)
-   **field\_n\_chromatic**, an *integer* indicating the number of
    chromatic sweeps for the latent field (by default 3)
-   **field\_n\_mala**, an *integer* indicating the number of HMC sweeps
    for the latent field (by default 1)

You can add some options about the diagnostic plots that pop up between
MCMC cycles

-   **plot\_diags**, a *boolean* to indicate if diagnostics for
    high-level parameters must be plotted (by default T)
-   **burn\_in**, a *number* between 0 and 1 to discard the first
    iterations in order to propose the diagnostics (by default 0.5)
-   **plot\_PSRF\_fields**, a *Boolean* to indicate if diagnostics must
    be plotted for the latent fields. Very costly. (by default F)
-   **debug\_outfile** = a file with messages from the parallel chains
    (by default “debug\_mcmc\_nngp.txt”)

<!-- -->

    mcmc_nngp_list_nonstat = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list_nonstat, n_iterations_update = 100, n_cycles = 1, plot_diags = F)

MCMC samples stored in the model
--------------------------------

Now that we have run the model, we can look up the rest of
mcmc\_nngp\_list.

#### **Records** is a *list* containing the saved states of the MCMC chains.

There are n\_chains elements in this list, each corresponding to a
chain. Each element is itself a *list* containing *arrays* of saved
parameters. The arrays correspond to the model parameters from the
“params” sub-list of each chain in “states”. The two first dimensions of
each arrays correspond to the original dimension of the parameters from
“params”. The samples are stacked along the third dimension.

    names(mcmc_nngp_list_nonstat$records$chain_1)

    ## [1] "range_beta"      "beta"            "noise_beta"      "noise_log_scale"
    ## [5] "noise_field"     "scale_beta"      "scale_log_scale" "scale_field"    
    ## [9] "field"

    # all MCMC samples of the first component of parameter beta
    mcmc_nngp_list_nonstat$records$chain_1$beta[1,,]

    ##  [1] 102.4682 102.5107 102.5408 102.6092 102.5844 102.5116 102.5508 102.4346
    ##  [9] 102.5471 102.5604

    # first MCMC samples of all components of parameter beta
    mcmc_nngp_list_nonstat$records$chain_1$beta[,,1]

    ## [1] 102.468154   4.988401

#### **iterations** is a *list* about the MCMC chain run.

It has 2 elements :

-   **checkpoints** is a *matrix* that is updated at each MCMC cycle.
    Its first column corresponds to the iteration, its second column
    corresponds to the time in minutes. Its first row correspond to the
    time at the end of model initialization.
-   **thinning** is a *vector of integers* that indicates which
    iterations were kept. The i-th slice from the arrays of Records
    corresponds to the MCMC state at the i-th element of thinning. So
    here the first slice of the array of records corresponds to the
    10-th iteration, the second to the 20-th, etc.

<!-- -->

    names(mcmc_nngp_list_nonstat$iterations)

    ## [1] "checkpoints" "thinning"

    mcmc_nngp_list_nonstat$iterations$checkpoints

    ##      iteration       time
    ## [1,]         0 0.02548741
    ## [2,]       100 0.83365432

    mcmc_nngp_list_nonstat$iterations$thinning

    ##  [1]  10  20  30  40  50  60  70  80  90 100

Prediction
----------

### Arguments

In order to predict at unobserved locations, the function
**Bidart::predict\_latent\_field** is used.

Its arguments are :

-   **mcmc\_nngp\_list**, a *list* corresponding to a fit model
    (diagnostic plots pop-up at each cycle).
-   **predicted\_locs**, a *matrix* of spatial locations where
    prediction needs to be done.
-   **X\_range\_pred**, a *data.frame* of covariates for the range at
    the predicted locations if covariates were specified in
    mcmc\_nngp\_initialize\_nonstationary (intercept is added
    automatically).
-   **X\_scale\_pred**, a *data.frame* of covariates for the marginal
    variance at the predicted locations if covariates were specified in
    mcmc\_nngp\_initialize\_nonstationary (intercept is added
    automatically).
-   **burn\_in**, a *number between 0 and 1* indicating the proportion
    of MCMC samples to discard (by default 0.5).
-   **n\_cores**, an *integer* indicating the number of cores.
-   **predict\_range**, a *Boolean* indicating whether each linear and
    random effects of range must be returned, their sum being returned
    though (by default F).
-   **predict\_scale**, a *Boolean* indicating whether each linear and
    random effects of marginal variance must be returned, their sum
    being returned though (by default F).

Here is just an example, remember that the model is not properly fit.

    predicted_locs = as.matrix(expand.grid(seq(0, 5, .1), seq(0, 5, .1)))
    pred = Bidart::predict_latent_field(mcmc_nngp_list = mcmc_nngp_list_nonstat, predicted_locs = predicted_locs, predict_scale = T)

### What is there in the prediction ?

#### **predicted\_samples** is a *list* containing the predictive MCMC samples.

It has n\_chains elements, one for each chain. Each of those *lists*
contains :

-   **field**, an *array* containing samples of the latent field stacked
    along the third dimension.

If the range is nonstationary, it also contains :

-   **log\_range**, an *array* containing samples of the log-range (or
    coordinates of the log range matrix) stacked along the third
    dimension.
-   if in addition predict\_range == T, samples of each linear effect
    and the random effect are included.

If the marginal variance of the latent field is nonstationary, it also
contains :

-   **log\_scale**, an *array* containing samples of the log-marginal
    variance stacked along the third dimension.
-   if in addition predict\_scale == T, samples of each linear effect
    and the random effect are included.

<!-- -->

    names(pred$predicted_samples$chain_1)

    ## [1] "scale_field"       "field"             "scale_(Intercept)"
    ## [4] "log_scale"

#### **summaries** is a *list* containing *arrays* that summary the MCMC predictive samples.

Those arrays have depth 1, and are virtually matrices. The statistics
are stacked along the first dimension of the array.

-   The first row is the **mean**
-   The second row is the **quandile 0.025**
-   The third row is the **median**
-   The fourth row is the **quandile 0.975**
-   The fift row is the **standard deviation**

<!-- -->

    pred$summaries$log_scale[,1:10,1]

    ##              [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
    ## mean    0.3311054  0.2680914  0.3008119  0.3135428  0.3061449  0.2999395
    ## q0.025 -0.2067922 -0.2279572 -0.2237408 -0.1847469 -0.1457625 -0.1824663
    ## median  0.2976207  0.2839132  0.3008806  0.3457849  0.3598187  0.2163061
    ## q0.975  0.7927676  0.6847423  0.7354910  0.8033166  0.6812568  0.7776041
    ## sd      0.3049020  0.2750104  0.3030261  0.3373110  0.2909610  0.3068144
    ##              [,7]       [,8]       [,9]       [,10]
    ## mean    0.2835988  0.2295540  0.2633066  0.37263648
    ## q0.025 -0.4937245 -0.7372746 -0.2200086 -0.07542769
    ## median  0.3675448  0.3241388  0.3039188  0.38948103
    ## q0.975  0.8286278  0.9489365  0.7715685  0.91756901
    ## sd      0.4061701  0.4852615  0.2971592  0.32460265
