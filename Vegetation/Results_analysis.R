# load
load("Vegetation/cleaned_data2_expanded.RData")
run_stat = readRDS("Vegetation/run_stat.RDS")
run_nonstat = readRDS("Vegetation/run_r_basis_10_nu=20.RDS")
KL_decomposition = readRDS("Vegetation/KL_decomposition.RDS") # loading karhunen loeve decomposition

# DIC
Bidart::DIC(run_stat)
Bidart::DIC(run_nonstat)

# Predictions and estimations

#### spatial basis
###set.seed(1)
###selected_locs = sample(seq(nrow(data_cleaned2)),20000)
###locs = cbind(data_cleaned2$scaled_x[selected_locs], data_cleaned2$scaled_y[selected_locs])
###reordering = GpGp::order_maxmin(locs)
###locs = locs[reordering,]
###
###selected_test_locs = sample(seq(nrow(data_cleaned2))[-selected_locs],20000)
###z_pred = data_cleaned2$NDVI[selected_test_locs]
###test_locs = cbind(data_cleaned2$scaled_x[selected_test_locs], data_cleaned2$scaled_y[selected_test_locs])
###locs_ = rbind(locs, test_locs)
###remove(data_cleaned2);gc()
###
###NNarray = GpGp::find_ordered_nn(locs_, 5)
###sparse_chol = Matrix::sparseMatrix(
###  i = row(NNarray)[!is.na(NNarray)],
###  j = NNarray[!is.na(NNarray)],
###  x = GpGp::vecchia_Linv(c(1, .1, 1, 0), "matern_isotropic", locs_, NNarray)[!is.na(NNarray)],
###)
###PP_basis = Matrix::solve(
###  sparse_chol, 
###  diag(1, nrow(sparse_chol), 1000)
###)
###
###KL_basis = (PP_basis %*% KL_decomposition$v[,seq(20)] %*% diag(1/KL_decomposition$d[seq(20)]))[-seq(run_nonstat$vecchia_approx$n_locs),seq(20)]
###saveRDS(KL_basis, "Vegetation/KL_basis.RDS")
###remove(PP_basis);gc()
###
#### Prediction
###pred_stat = Bidart::predict_latent_field(mcmc_nngp_list = run_stat, predicted_locs = test_locs, n_cores = 3, burn_in = .1)
###saveRDS(pred_stat, "Vegetation/pred_stat.RDS")
###gc()
###pred_nonstat = Bidart::predict_latent_field(mcmc_nngp_list = run_nonstat, predicted_locs = test_locs, n_cores = 3, X_range_pred = KL_basis, predict_range = T, predict_scale = T)
###saveRDS(pred_nonstat, "Vegetation/pred_nonstat.RDS")
###gc()
###
#### Estimation
###estimation_stat = Bidart::estimate_parameters(run_stat)
###saveRDS(estimation_stat, "Vegetation/estimation_stat.RDS")
###gc()
###estimation_nonstat = Bidart::estimate_parameters(run_nonstat)
###saveRDS(estimation_nonstat, "Vegetation/estimation_nonstat.RDS")
###gc()


estimation_nonstat = readRDS("Vegetation/estimation_nonstat.RDS")

# plot ellipses


Bidart::plot_ellipses(locs = run_nonstat$data$observed_locs[seq(1000),], log_range = cbind(1, KL_decomposition$u[,seq(20)])[seq(1000),] %*% estimation_nonstat$range_beta[1,,], shrink = .01)
#Bidart::plot_ellipses(locs = run_nonstat$data$locs[seq(1000),], log_range = run_nonstat$data$covariates$range_X$X_locs[seq(1000),] %*% run_nonstat$states$chain_3$params$range_beta, shrink = .01)
#Bidart::plot_ellipses(locs = run_nonstat$data$locs[seq(1000),], log_range = run_nonstat$data$covariates$range_X$X_locs[seq(1000),] %*% run_nonstat$states$chain_2$params$range_beta, shrink = .01)
#Bidart::plot_ellipses(locs = run_nonstat$data$locs[seq(1000),], log_range = run_nonstat$data$covariates$range_X$X_locs[seq(1000),] %*% run_nonstat$states$chain_1$params$range_beta, shrink = .01)


#BHEP test
records = NULL
for(chain_idx in seq(3))
{
  for(x in run_nonstat$records[chain_idx])
  {
    records = rbind(records, 
                    do.call(cbind, 
                            list(
                              t(x$range_beta[,1,-seq(100)]),
                              t(x$range_beta[,2,-seq(100)]),
                              t(x$range_beta[,3,-seq(100)]),
                              x$noise_beta[,,-seq(100)],
                              x$scale_beta[,,-seq(100)])
                    ))
  }
}

test = mnt::test.BHEP(data = records)
test$Decision
saveRDS(test, "Vegetation/test_BHEP.RDS")

plot(c(-4, 4), c(-4, 4))
for( i in seq(ncol(records)))lines(x = quantile(x = (records[,i]-mean(records[,i]))/sd(records[,i]), probs = seq(.01, .99, .01)), y = qnorm(seq(.01, .99, .01)))
abline(a = 0, b = 1, col = 2)


qqnorm(records[,5])
qqnorm(records[,77])

plot(density(records[,1]))
plot(records[,1])
plot(records[,2])

