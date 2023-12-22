#################################
# Options to set before the run #
#################################

## --------- IMPORTANT ---------- ##
# For some reason you have to set manually num_threads_per_chain in the function mcmc_nngp_run

# working directory
# setwd("NDVI")
# directory if you want/need to install manually libraries, put NULL else
install_the_libraries =F
# directory if you want/need to install manually libraries, put NULL else
libraries_dir = "../R_packages/"# NULL
# rounding the coordinates allows to save computation effort at
# the expense of some accuracy of the coordinates
# this parameter controls the rounding of the locations. the smaller the more rounding
rounding_factor = 1/350
# Does some informative plots prior to MCMC run
plot_stuff = F

############################
# Installation and loading #
############################
# load libraries from custom folder in case default is not writable 
# like in my SLURM cluster
if(!is.null(libraries_dir)){
  if(install_the_libraries){
    install.packages(
      c("GpGp","FNN","abind","fields","parallel","Matrix", "expm"), 
      lib = libraries_dir, 
      repos = c(CRAN = "https://cloud.r-project.org")
    )
    install.packages("Bidart_1.0.tar.gz", lib = libraries_dir)
  }
  library(GpGp, lib.loc = libraries_dir);library(FNN, lib.loc = libraries_dir);
  library(Bidart, lib.loc = libraries_dir);library(abind, lib.loc = libraries_dir);
  library(fields, lib.loc = libraries_dir);library(parallel, lib.loc = libraries_dir);
  library(Matrix, lib.loc = libraries_dir);
  }
# load libraries from default folder
if(is.null(libraries_dir)){
  if(install_the_libraries){
    install.packages(
      c("GpGp","FNN","abind","fields","parallel","Matrix"), 
      repos = c(CRAN = "https://cloud.r-project.org")
    )
    install.packages("Bidart_1.0.tar.gz")
  }
  library(GpGp);library(FNN);library(Bidart);library(abind);library(fields);library(parallel);library(Matrix)
}


########################
# Loading the data set #
########################
set.seed(1)
# load the data set
#load(file = "data_cleaned_small_expanded.RData")
load(file ="data_cleaned_small_expanded.RData")
response_variable = data_cleaned_small$NDVI
response_variable = as.vector(scale(response_variable))
locs = cbind(data_cleaned_small$scaled_x, data_cleaned_small$scaled_y)

########################
# rounding coordinates #
########################
locs  = rounding_factor * round(locs/rounding_factor)
print("total number of observations")
print(nrow(locs))
print("number of rounded locations")
print(sum(!duplicated(locs)))
#########################
# train-test separation #
#########################

# matching unique locs with all locs
unique_locs = unique(split(locs, row(locs)))
locs_match_idx = match(split(locs, row(locs)), unique_locs)
# selecting separated locs using GpGp order maxmin 
selected_locs = GpGp::order_maxmin(do.call(rbind, unique_locs))[seq(10000)]
selected_locs = which(locs_match_idx %in% selected_locs)
# representing test locs
if(plot_stuff){
  if(plot_stuff)plot(locs[selected_locs,], pch = 16, cex = .1, xlab = "", ylab="", main = "test locs")
}
# indexing data sets
train_locs = locs[-selected_locs,]
train_field = response_variable[-selected_locs]
saveRDS(train_field, "train_field.RDS")
saveRDS(response_variable[selected_locs],  "test_field.RDS",)
saveRDS(train_locs, "train_locs.RDS")
saveRDS(rounding_factor * round(locs[selected_locs,]/rounding_factor), "test_locs.RDS")


##############################
# plotting interest variable #
##############################
if(plot_stuff){
  par(mfrow = c(1,2))
  Bidart::plot_pointillist_painting(
    locs, 
    response_variable, 
    cex = .2, 
    main = "Full data"
  )
  Bidart::plot_pointillist_painting(
    rounding_factor * round(locs/rounding_factor), 
    response_variable, cex = .22, 
    main = "With rounding"
  )
  par(mfrow = c(1,1))
}

####################
# Getting PP basis #
####################

PP = Bidart::get_PP(
  observed_locs = train_locs, # spatial sites
  matern_range = .2,
  n_PP = 50, # number of knots
  m = 15 # number of NNGP parents
  )
saveRDS(PP, "PP.RDS")

if(plot_stuff){
  # comparison between PP and NNGP
  seed_vector =  rnorm(PP$n_PP + nrow(PP$unique_reordered_locs))
  par(mfrow = c(1,2))
  Bidart::plot_pointillist_painting(train_locs, Bidart::X_PP_mult_right(PP = PP, use_PP = T, Y = seed_vector[seq(PP$n_PP)]), cex = .3, main ="NNGP into PP")
  points(PP$knots, pch = 3, cex = .3)
  Bidart::plot_pointillist_painting(rbind(PP$knots, PP$unique_reordered_locs), as.vector(Matrix::solve(PP$sparse_chol, seed_vector)), cex = .3, main = "NNGP")
  par(mfrow = c(1,1))
}

###############################################################
# Nonstat aniso + heteroskedastic run   with additional noise #
###############################################################

if(length(grep("run_nonstat_aniso_heterosk_additional_noise.RDS", list.files()))==0){
  saveRDS(1, "run_nonstat_aniso_heterosk_additional_noise.RDS")
  mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
    observed_locs = train_locs, 
    observed_field = train_field + rnorm(length(train_field)), 
    X = as.data.frame(train_locs), 
    m = 10, 
    nu = 1.5, 
    range_PP = T, 
    noise_PP = T, 
    anisotropic = T, 
    sphere = F, 
    n_chains = 2, 
    PP = PP, seed = 2, 
    noise_log_scale_prior = c(-3, 5)
  )
  for(i in seq(30)){
    mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(
      mcmc_nngp_list,
      n_cores = 2, 
      num_threads_per_chain = ceiling(parallel::detectCores()/2), 
      lib.loc = libraries_dir, 
      plot_diags = T
    )
    saveRDS(mcmc_nngp_list, "run_nonstat_aniso_heterosk_additional_noise.RDS")
    pdf("diags_nonstat_aniso_heterosk_additional_noise.pdf")
    Bidart::diagnostic_plots(mcmc_nngp_list)
    dev.off()
  }
}

#######################################
# Nonstat aniso + heteroskedastic run #
#######################################


if(length(grep("run_nonstat_aniso_heterosk.RDS", list.files()))==0){
  saveRDS(1, "run_nonstat_aniso_heterosk.RDS")
  mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
    observed_locs = train_locs, 
    observed_field = train_field, 
    X = as.data.frame(train_locs), 
    m = 10, 
    nu = 1.5, 
    range_PP = T, 
    noise_PP = T, 
    anisotropic = T, 
    sphere = F, 
    n_chains = 2, 
    PP = PP, 
    noise_log_scale_prior = c(-3, 5)
  )
  for(i in seq(15)){
    mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(
      mcmc_nngp_list,
      n_cores = 2, 
      num_threads_per_chain = ceiling(parallel::detectCores()/(2)), 
      lib.loc = libraries_dir, 
      plot_diags = T
    )
    saveRDS(mcmc_nngp_list, "run_nonstat_aniso_heterosk.RDS")
    pdf("diags_nonstat_aniso_heterosk.pdf")
    Bidart::diagnostic_plots(mcmc_nngp_list)
    dev.off()
  }
}


##################
# Stat aniso run #
##################

if(length(grep("run_stat_aniso.RDS", list.files()))==0){
  saveRDS(1, "run_stat_aniso.RDS")
  mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
    observed_locs = train_locs, 
    observed_field = train_field, 
    X = as.data.frame(train_locs), 
    m = 10, 
    nu = 1.5, 
    range_PP = F, 
    noise_PP = F, 
    anisotropic = T, sphere = F, 
    n_chains = 2,  
    PP = PP, 
    noise_log_scale_prior = c(-3, 5)
  )
  for(i in seq(10)){
    mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(
      mcmc_nngp_list,
      n_cores = 2, 
      num_threads_per_chain = ceiling(parallel::detectCores()/2),
      lib.loc = libraries_dir
    )
    saveRDS(mcmc_nngp_list, "run_stat_aniso.RDS")
    pdf("diags_stat_aniso.pdf")
    Bidart::diagnostic_plots(mcmc_nngp_list)
    dev.off()
  }
}

#######################################
# Stat aniso + heteroskedasticity run #
#######################################

if(length(grep("run_stat_aniso_heterosk.RDS", list.files()))==0){
  saveRDS(1, "run_stat_aniso_heterosk.RDS")
  mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
    observed_locs = train_locs, 
    observed_field = train_field, 
    X = as.data.frame(train_locs), 
    m = 10, 
    nu = 1.5, 
    range_PP = F, 
    noise_PP = T, 
    anisotropic = T, sphere = F, 
    n_chains = 2,  
    PP = PP, 
    noise_log_scale_prior = c(-3, 5)
  )
  for(i in seq(15)){
    mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(
      mcmc_nngp_list,
      n_cores = 2, 
      num_threads_per_chain = ceiling(parallel::detectCores()/2),
      lib.loc = libraries_dir
    )
    saveRDS(mcmc_nngp_list, "run_stat_aniso_heterosk.RDS")
    pdf("diags_stat_aniso_heterosk.pdf")
    Bidart::diagnostic_plots(mcmc_nngp_list)
    dev.off()
  }
}

################
# Stat iso run #
################

if(length(grep("run_stat_iso.RDS", list.files()))==0){
  mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
    observed_locs = train_locs,  
    observed_field = train_field, 
    X = as.data.frame(train_locs),  
    m = 10,  
    nu = 1.5,  
    range_PP = F,  
    noise_PP = F,  
    anisotropic = F, sphere = F, 
    n_chains = 2, 
    PP = PP, 
    noise_log_scale_prior = c(-3, 5)
  )
  for(i in seq(10)){
    mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(
      mcmc_nngp_list,
      n_cores = 2, 
      num_threads_per_chain = ceiling(parallel::detectCores()/2),
      lib.loc = libraries_dir
    )
    saveRDS(mcmc_nngp_list, "run_stat_iso.RDS")
    pdf("diags_stat_iso.pdf")
    Bidart::diagnostic_plots(mcmc_nngp_list)
    dev.off()
  }
}

########################################
# Stat iso run with heteroskadasticity #
########################################

if(length(grep("run_stat_iso_heterosk.RDS", list.files()))==0){
  mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
    observed_locs = train_locs,  
    observed_field = train_field, 
    X = as.data.frame(train_locs),  
    m = 10,  
    nu = 1.5,  
    range_PP = F,  
    noise_PP = T,  
    anisotropic = F, sphere = F, 
    n_chains = 2, 
    PP = PP, 
    noise_log_scale_prior = c(-3, 5)
  )
  for(i in seq(10)){
    mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary_socket(
      mcmc_nngp_list,
      n_cores = 2, 
      num_threads_per_chain = ceiling(parallel::detectCores()/2),
      lib.loc = libraries_dir
    )
    saveRDS(mcmc_nngp_list, "run_stat_iso_heterosk.RDS")
    pdf("diags_stat_iso_heterosk.pdf")
    Bidart::diagnostic_plots(mcmc_nngp_list)
    dev.off()
  }
}




