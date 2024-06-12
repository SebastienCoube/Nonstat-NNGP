#setwd("matern_samples_ex/")
set.seed(1)
locs = as.matrix(expand.grid(seq(0, 1, .005), seq(0, 1, .005)))
locs = locs[GpGp::order_maxmin(locs),]
PP = Bidart::get_PP(locs, .1, n_PP = 50)
Bidart::compare_PP_NNGP(PP)
range_beta = matrix(.25*rnorm(153), 51)
range_beta[1,1] = -4.5
NNarray = GpGp::find_ordered_nn(locs, 10)

sparse_chol_aniso = 
  Matrix::sparseMatrix(
    i = row(NNarray)[!is.na(NNarray)], 
    j = NNarray[!is.na(NNarray)], 
    x = 
      Bidart::compute_sparse_chol(
        range_beta = range_beta, 
        NNarray = NNarray, locs = locs, 
        anisotropic = T, 
        sphere = F,
        PP = PP, use_PP = T, 
        range_X = matrix(1, nrow(locs)), 
        compute_derivative = F, 
        nu = 1.5, 
        locs_idx = NULL, 
        num_threads = 10
      )[[1]][!is.na(NNarray)], 
    triangular = T
  )

sparse_chol_iso = 
  Matrix::sparseMatrix(
    i = row(NNarray)[!is.na(NNarray)], 
    j = NNarray[!is.na(NNarray)], 
    x = 
      Bidart::compute_sparse_chol(
        range_beta = range_beta[,1,drop=F], 
        NNarray = NNarray, locs = locs, 
        anisotropic = F, 
        sphere = F,
        PP = PP, use_PP = T, 
        range_X = matrix(1, nrow(locs)), 
        compute_derivative = F, 
        nu = 1.5, 
        locs_idx = NULL, 
        num_threads = 10
      )[[1]][!is.na(NNarray)], 
    triangular = T
  )

log_range = 
  (
    Bidart::X_PP_mult_right(
      X = matrix(1, nrow(locs)), 
      PP = PP, use_PP = T, 
      locs_idx = NULL, 
      Y = range_beta
    )
  )

# tatato = GpGp::fast_Gp_sim(covparms = c(.5, .1, 1.5, 0), locs = locs, m = 10, covfun_name = "matern_isotropic")-4
# log_range = cbind(tatato, tatato, 0)

z = rnorm(nrow(locs))
w_aniso = as.vector(Matrix::solve(sparse_chol_aniso, z))
w_iso = as.vector(Matrix::solve(sparse_chol_iso, z))
Bidart::plot_pointillist_painting(locs, log_range[,1])

png("aniso.png", width = 2000, height = 2000)
par(mar = c(0,0,0,0))
plot(0,0, xlim = c(-.04, 1.04), ylim = c(-.04, 1.04), xaxt='n', yaxt='n', xlab = "", ylab = "", type = "n",bty="n")
Bidart::plot_pointillist_painting(locs, w_aniso, cex = 1.8, add = T)
dev.off()


png("iso.png", width = 2000, height = 2000)
par(mar = c(0,0,0,0))
plot(0,0, xlim = c(-.04, 1.04), ylim = c(-.04, 1.04), xaxt='n', yaxt='n', xlab = "", ylab = "", type = "n", bty="n")
Bidart::plot_pointillist_painting(locs, w_iso, cex = 1.8, add = T)
dev.off()


estep = .1
locs_4_ellipses = as.matrix(expand.grid(seq(0, 1, estep), seq(0, 1, estep)))
ellipse_idx = match(split(locs_4_ellipses, row(locs_4_ellipses)), split(locs, row(locs)))

png("ellipses_aniso.png", width = 2000, height = 2000)
par(mar = c(0,0,0,0), lwd = 5)
plot(0,0, xlim = c(-.04, 1.04), ylim = c(-.04, 1.04), xaxt='n', yaxt='n', xlab = "", ylab = "", type = "n", bty="n")
plot_ellipses(locs_4_ellipses, log_range[ellipse_idx,], shrink = sqrt(8*nu), add = T)
dev.off()

png("ellipses_iso.png", width = 2000, height = 2000)
par(mar = c(0,0,0,0), lwd = 5)
plot(0,0, xlim = c(-.04, 1.04), ylim = c(-.04, 1.04), xaxt='n', yaxt='n', xlab = "", ylab = "", type = "n", bty="n")
plot_ellipses(locs_4_ellipses, log_range[ellipse_idx,1, drop=F], shrink = sqrt(8*nu), add = T)
dev.off()

