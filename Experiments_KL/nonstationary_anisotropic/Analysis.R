res = readRDS("Experiments_KL/nonstationary_anisotropic/res_elliptic.RDS")
res

res = res[-which(res[,4]==1),]
res


KL = as.numeric(res[,6])
sc = res[,2]
ordering = res[,3]
m = factor(res[,4], levels = c("5", "10", "20"))
xtable::xtable(summary(lm(KL~sc +ordering*m)))




res = readRDS("Experiments_KL/nonstationary_anisotropic/res_elliptic.RDS")
res
res = res[-which(res[,4]==1 | res[,2]!=.5),]
res
KL = as.numeric(res[,6])
sc = res[,2]
ordering = res[,3]
m = factor(res[,4], levels = c("5", "10", "20"))
boxplot(KL~ordering +m)
