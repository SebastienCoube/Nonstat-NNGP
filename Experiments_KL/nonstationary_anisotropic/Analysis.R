res = readRDS("Experiments_KL/nonstationary_anisotropic/res_elliptic.RDS")
res

KL = as.numeric(res[,6])
sc = res[,2]
ordering = res[,3]
m = factor(res[,4], levels = c("5", "10", "20"))
xtable::xtable(summary(lm(KL~sc +ordering*m)), digits = 2)
boxplot(KL~ordering +m)
