
######################
# nonstat with range #
######################

load("Vegetation/cleaned_data2_expanded.RData")
set.seed(1)
selected_locs = sample(seq(nrow(data_cleaned2)),20000)

# extracting locations
locs = cbind(data_cleaned2$scaled_x[selected_locs], data_cleaned2$scaled_y[selected_locs])
z = data_cleaned2$NDVI[selected_locs]
reordering = GpGp::order_maxmin(locs)
locs = locs[reordering,]
z = z[reordering]

# spatial basis
#NNarray = GpGp::find_ordered_nn(locs, 5)
#sparse_chol = Matrix::sparseMatrix(
#  i = row(NNarray)[!is.na(NNarray)],
#  j = NNarray[!is.na(NNarray)],
#  x = GpGp::vecchia_Linv(c(1, .1, 1, 0), "matern_isotropic", locs, NNarray)[!is.na(NNarray)],
#)
#PP_basis = Matrix::solve(
#  sparse_chol, 
#  diag(1, nrow(sparse_chol), 1000)
#)
#SVD = irlba::irlba(A = PP_basis, nu = 50, nv = 50)
#saveRDS(SVD, "Vegetation/KL_decomposition.RDS")
SVD = readRDS("Vegetation/KL_decomposition.RDS")

plot(SVD$d^2)
abline(h = 0)
abline(v = 20)


Bidart::plot_pointillist_painting(locs, SVD$u[,1])
Bidart::plot_pointillist_painting(locs, SVD$u[,2])
Bidart::plot_pointillist_painting(locs, SVD$u[,3])
Bidart::plot_pointillist_painting(locs, SVD$u[,4])
Bidart::plot_pointillist_painting(locs, SVD$u[,5])
Bidart::plot_pointillist_painting(locs, SVD$u[,10])
Bidart::plot_pointillist_painting(locs, SVD$u[,30])
Bidart::plot_pointillist_painting(locs, SVD$u[,40])
Bidart::plot_pointillist_painting(locs, SVD$u[,50])

par(mfrow = c(2, 1))
t1 =Sys.time()
set.seed(t1)
k = 20
Bidart::plot_pointillist_painting(locs, SVD$u[,seq(k)] %*% (SVD$d[seq(k)] * rnorm(k)))
set.seed(t1)
k = 50
Bidart::plot_pointillist_painting(locs, SVD$u[,seq(k)] %*% (SVD$d[seq(k)] * rnorm(k)))

mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
  observed_locs = locs, 
  observed_field = z, 
  m = 5, reordering = c("maxmin"), covfun = "nonstationary_exponential_anisotropic", 
  range_X = as.data.frame(SVD$u[,seq(20)])
  )

#mcmc_nngp_list = readRDS("Vegetation/run_nsr_basis_10.RDS")
for(i in seq(80)){
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_iterations_update = 100, n_cycles = 1, swap_range_scale = F, field_n_chromatic = 3, field_n_mala = 1)
}
#saveRDS(mcmc_nngp_list, "Vegetation/run_r_basis_10_nu=20.RDS")




###################################
# nonstat with range + many knots #
###################################

load("Vegetation/cleaned_data2_expanded.RData")
set.seed(1)
selected_locs = sample(seq(nrow(data_cleaned2)),20000)

# extracting locations
locs = cbind(data_cleaned2$scaled_x[selected_locs], data_cleaned2$scaled_y[selected_locs])
z = data_cleaned2$NDVI[selected_locs]
reordering = GpGp::order_maxmin(locs)
locs = locs[reordering,]
z = z[reordering]

# spatial basis
#NNarray = GpGp::find_ordered_nn(locs, 5)
#sparse_chol = Matrix::sparseMatrix(
#  i = row(NNarray)[!is.na(NNarray)],
#  j = NNarray[!is.na(NNarray)],
#  x = GpGp::vecchia_Linv(c(1, .1, 1, 0), "matern_isotropic", locs, NNarray)[!is.na(NNarray)],
#)
#PP_basis = Matrix::solve(
#  sparse_chol, 
#  diag(1, nrow(sparse_chol), 1000)
#)
#SVD = irlba::irlba(A = PP_basis, nu = 50, nv = 50)
#saveRDS(SVD, "Vegetation/KL_decomposition.RDS")
SVD = readRDS("Vegetation/KL_decomposition.RDS")

plot(SVD$d^2)
abline(h = 0)
abline(v = 20)


Bidart::plot_pointillist_painting(locs, SVD$u[,1])
Bidart::plot_pointillist_painting(locs, SVD$u[,2])
Bidart::plot_pointillist_painting(locs, SVD$u[,3])
Bidart::plot_pointillist_painting(locs, SVD$u[,4])
Bidart::plot_pointillist_painting(locs, SVD$u[,5])
Bidart::plot_pointillist_painting(locs, SVD$u[,10])
Bidart::plot_pointillist_painting(locs, SVD$u[,30])
Bidart::plot_pointillist_painting(locs, SVD$u[,40])
Bidart::plot_pointillist_painting(locs, SVD$u[,50])

par(mfrow = c(2, 1))
t1 =Sys.time()
set.seed(t1)
k = 20
Bidart::plot_pointillist_painting(locs, SVD$u[,seq(k)] %*% (SVD$d[seq(k)] * rnorm(k)))
set.seed(t1)
k = 50
Bidart::plot_pointillist_painting(locs, SVD$u[,seq(k)] %*% (SVD$d[seq(k)] * rnorm(k)))

mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
  observed_locs = locs, 
  observed_field = z, 
  m = 5, reordering = c("maxmin"), covfun = "nonstationary_exponential_anisotropic", 
  range_X = as.data.frame(SVD$u[,seq(40)])
  )

for(i in seq(80)){
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_iterations_update = 100, n_cycles = 1, swap_range_scale = F, field_n_chromatic = 3, field_n_mala = 1)
}
saveRDS(mcmc_nngp_list, "Vegetation/run_r_basis_10_nu=40.RDS")




##############
# Stationary #
##############



load("Vegetation/cleaned_data2_expanded.RData")

set.seed(1)
selected_locs = sample(seq(nrow(data_cleaned2)),20000)


# extracting locations
locs = cbind(data_cleaned2$scaled_x[selected_locs], data_cleaned2$scaled_y[selected_locs])
z = data_cleaned2$NDVI[selected_locs]
reordering = GpGp::order_maxmin(locs)
locs = locs[reordering,]
z = z[reordering]



mcmc_nngp_list = Bidart::mcmc_nngp_initialize_nonstationary(
  observed_locs = locs, 
  observed_field = z, 
  m = 5, reordering = c("maxmin"), covfun = "exponential_anisotropic2D"
)


for(i in seq(10)){
  mcmc_nngp_list = Bidart::mcmc_nngp_run_nonstationary(mcmc_nngp_list, n_iterations_update = 100, n_cycles = 1, debug_outfile = NULL, swap_range_scale = T)
}
saveRDS(mcmc_nngp_list, "Vegetation/run_stat.RDS")

Bidart::diagnostic_plots(mcmc_nngp_list, burn_in = .15)


