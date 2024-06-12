setwd("matern_samples_ex/")

set.seed(3)
locs = cbind(seq(0, 1, .005), 1)
z = rnorm(nrow(locs))

pdf("matern_samples_smoothness.pdf", width = 7)
plot(
  locs[,1], 
  t(chol(GpGp::matern_isotropic(c(1,1,.5,0), locs))) %*% z, 
  type = "l", xlab = "spatial location", ylab = "Gaussian Processes"
)
lines(
  locs[,1], 
  t(chol(GpGp::matern_isotropic(c(1,.5,1,0), locs))) %*% z, 
  col=2
)
lines(
  locs[,1], 
  t(chol(GpGp::matern_isotropic(c(1,.3,1.5,0), locs))) %*% z, 
  col=4
)
legend("topleft", fill = c(1,2,4), 
       legend = c(
         "range = 1,    smoothness = 0.5",
         "range = 0.5, smoothness = 1",
         "range = 0.3, smoothness = 1.5"
         ))
dev.off()


set.seed(3)

pdf("matern_samples_var.pdf", width = 7)
par(mar = c(2,2,.1,.1))
z = rnorm(nrow(locs))
plot(
  locs[,1], 
  t(chol(GpGp::matern_isotropic(c(10,.1,.5,0), locs))) %*% z, 
  type = "l", xlab = "spatial location", ylab = "Gaussian processes"
)
#z = rnorm(nrow(locs))
lines(
  locs[,1], 
  t(chol(GpGp::matern_isotropic(c(1,.1,.5,0), locs))) %*% z, 
  col = 2
)
#z = rnorm(nrow(locs))
lines(
  locs[,1], 
  t(chol(GpGp::matern_isotropic(c(.1,.1,.5,0), locs))) %*% z, 
  col = 4
)


legend("topleft", fill = c(1,2,4), 
       legend = c(
         "marginal var = 10",
         "marginal var = 1",
         "marginal var = 0.1"
       ))
dev.off()

set.seed(3)

pdf("matern_samples_range.pdf", width = 7)
par(mar = c(2,2,.1,.1))
z = rnorm(nrow(locs))
plot(
  locs[,1], 
  t(chol(GpGp::matern_isotropic(c(1,.01,.5,0), locs))) %*% z, 
  type = "l", xlab = "spatial location", ylab = "Gaussian processes"
)
#z = rnorm(nrow(locs))
lines(
  locs[,1], 
  t(chol(GpGp::matern_isotropic(c(1,.1,.5,0), locs))) %*% z, 
  col = 2
)
#z = rnorm(nrow(locs))
lines(
  locs[,1], 
  t(chol(GpGp::matern_isotropic(c(1,1,.5,0), locs))) %*% z, 
  col = 4
)


legend("topleft", fill = c(1,2,4), 
       legend = c(
         "range = .01",
         "range = .1",
         "range = 1"
       ))
dev.off()
