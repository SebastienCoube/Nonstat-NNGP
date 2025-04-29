setwd("sensitivity_PP_range/")
res = as.data.frame(t(sapply(list.files()[grep("range_exp_ended", list.files())], readRDS )))
input = readRDS("range_design.RDS")

log_r = log(input[sapply(row.names(res), function(x)as.numeric(substr(x, 17, nchar(x)))),1])
plot(log_r, res$w_smooth_MSE)
plot(log_r, res$w_pred_MSE)
plot(log_r, res$log_range_MSE)
plot(log_r, res$log_range_scale)


boxplot(res$log_range_scale ~ round(log_r, 3))
boxplot(res$log_range_MSE ~ round(log_r, 3))
boxplot(res$pred_elpd ~ round(log_r, 3))
boxplot(res$obs_elpd ~ round(log_r, 3))


log_r_factor =  as.factor(round(log_r, 3))
levels(log_r_factor) = c("_log(0.05)", "_log(0.05 * sqrt(2))", "_log(0.1)", "_log(0.1 * sqrt(2)) " ,  "_log(0.2)")
log_r_factor = relevel(log_r_factor, ref = "_log(0.1)")


out_ = signif(cbind(
  summary(lm(res$log_range_scale ~log_r_factor + as.factor(res$seed.Var2)))$coefficients[seq(5), c(1, 4)],
  summary(lm(res$log_range_MSE ~log_r_factor + as.factor(res$seed.Var2)))$coefficients[seq(5), c(1, 4)],
  summary(lm(res$obs_elpd ~log_r_factor + as.factor(res$seed.Var2)))$coefficients[seq(5), c(1, 4)],
  summary(lm(res$pred_elpd ~log_r_factor + as.factor(res$seed.Var2)))$coefficients[seq(5), c(1, 4)]
), 2)
out=  as.data.frame(out_)
for(i in seq(nrow(out))){
  for(j in seq(ncol(out)/2)){
    if(out_[i, 2*j]<2e-16){
      out[i, 2*j-1] = paste("\textbf{", out_[i, 2*j-1], "}")
      out[i, 2*j] = paste("\textbf{", out_[i, 2*j], "}")
    }
  }
}

out = rbind(c(rep(c("Estimate", "p-value"), 4)), out)

for(i in seq(nrow(out))){
  for(j in seq(ncol(out))){
    out[i, j] = paste("$", out[i, j], "$")
  }
}

for(i in seq(nrow(out))){
  for(j in seq(ncol(out)-1)){
    out[i, j] = paste(out[i, j], "&")
  }
  out[i, ncol(out)] = paste(out[i, ncol(out)], "\\")
}
row.names(out) = NULL

out
