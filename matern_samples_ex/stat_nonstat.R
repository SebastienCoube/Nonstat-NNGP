setwd("matern_samples_ex/")
locs = cbind(seq(0, 1, .001), 1)
w = GpGp::fast_Gp_sim(c(1,.002, 1.5, 0), locs = locs)
par(mar = c(1,1,1,1))
pdf("stationary.pdf")
plot(locs[,1], w, type = "l", xlab = "spatial site", ylab = "process", main = "")
dev.off()
w_ = w * (1+3*pnorm((locs[,1]-.5)*8))
pdf("nonstationary.pdf")
plot(locs[,1], w_, type = "l", xlab = "spatial site", ylab = "process", main = "")
dev.off()

library(abind, )