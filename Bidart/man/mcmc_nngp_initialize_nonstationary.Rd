% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_nngp_initialize_nonstationary.R
\name{mcmc_nngp_initialize_nonstationary}
\alias{mcmc_nngp_initialize_nonstationary}
\title{Initializes a list containing the observations, the hierarchical model, and the MCMC chains.}
\usage{
mcmc_nngp_initialize_nonstationary(
  observed_locs = NULL,
  observed_field = NULL,
  X = NULL,
  m = 10,
  nu = 1.5,
  anisotropic = F,
  sphere = F,
  PP = NULL,
  n_chains = 2,
  noise_PP = F,
  noise_X = NULL,
  noise_beta_mean = NULL,
  noise_beta_precision = NULL,
  noise_log_scale_prior = NULL,
  scale_PP = F,
  scale_X = NULL,
  scale_beta_mean = NULL,
  scale_beta_precision = NULL,
  scale_log_scale_prior = NULL,
  range_PP = F,
  range_X = NULL,
  range_beta_mean = NULL,
  range_beta_precision = NULL,
  range_log_scale_prior = NULL,
  seed = 1
)
}
\arguments{
\item{observed_locs}{a matrix of spatial coordinates where observations are done}

\item{observed_field}{a vector of observations of the interest variable}

\item{X}{a data.frame of covariates explaining the interest variable through fixed linear effects}

\item{m}{number of nearest neighbors to do Vecchia's approximation}

\item{nu}{Matern smoothness, either 0.5 or 1.5}

\item{anisotropic}{anisotropic covariance}

\item{sphere}{Boolean, indicating lon-lat data}

\item{noise_X}{a data.frame of covariates explaining the Gaussian noise variance through fixed linear effects}

\item{noise_beta_mean}{vector indicating the prior mean for the regression coefficients explaining the Gaussian noise variance through fixed linear effects. Filled automatically if NULL. The number of rows is equal to the number of variables in X_noise after adding an intercept and expanding the factors, the number of columns is 1.}

\item{noise_beta_precision}{matrix for the prior precision for the regression coefficients explaining the Gaussian noise variance through fixed linear effects. Filled automatically if NULL.}

\item{noise_log_scale_prior}{1 times 2 matrix for the prior on the log-variance of the noise PP field.}

\item{scale_X}{a data.frame of covariates explaining the Gaussian process marginal variance through fixed linear effects}

\item{scale_beta_mean}{vector indicating the prior mean for the regression coefficients explaining the Gaussian process marginal variance through fixed linear effects. Filled automatically if NULL. The number of rows is equal to the number of variables in X_scale after adding an intercept, the number of columns is 1.}

\item{scale_beta_precision}{matrix for the prior precision for the regression coefficients explaining the Gaussian process marginal variance through fixed linear effects. Filled automatically if NULL.}

\item{scale_log_scale_prior}{1 times 2 matrix for the prior on the log-variance of the scale PP field.}

\item{range_X}{a data.frame of covariates explaining the Gaussian process range through fixed linear effects}

\item{range_beta_mean}{vector (3-columns matrix in anisotropy case) indicating the prior mean for the regression coefficients explaining the Gaussian process range through fixed linear effects. Filled automatically if NULL.The number of rows is equal to the number of variables in X_range after adding an intercept, the number of columns is 1 if isotropic function or 3 if anisotropic function.}

\item{range_beta_precision}{matrix for the prior precision for the regression coefficients explaining the Gaussian process range through fixed linear effects. Filled automatically if NULL. In anisotropic case, the matrix has size 3*n_var, each 3-block indicating the precision for a determinant-direction-direction 3-uplet: intecrept-det, intercept-dir, intercept-dir, V1-det, V1-dir, V1-dir, etc...}

\item{range_log_scale_prior}{1 times 2 matrix for the prior on the log-variance of the range PP field. 
In the case of anisotropic range, input an 3 times 2 matrix, indicating bounds for the eigenvalues of the trivariate log-variance matrix.}
}
\description{
Initializes a list containing the observations, the hierarchical model, and the MCMC chains.
}
