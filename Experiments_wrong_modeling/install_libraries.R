libraries_dir = "../R_packages"
install.packages(
#  #c("GpGp","FNN","abind","fields","parallel","Matrix", "expm"), 
#  c("Rcpp"), 
  c("RcppArmadillo"), 
  lib = libraries_dir, 
  repos = c(CRAN = "https://cloud.r-project.org")
)
library(Rcpp, lib.loc = libraries_dir)
library(RcppArmadillo, lib.loc = libraries_dir)
install.packages("Bidart_1.0.tar.gz", lib = libraries_dir)
