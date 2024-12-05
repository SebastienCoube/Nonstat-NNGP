
setwd("NDVI_r3/")
res = do.call(
  rbind, 
  list(
    NNGP_full_monty =        readRDS("res/NNGP_full_monty_score.RDS"),
    NNGP_aniso_heterosk =    readRDS("res/NNGP_aniso_heterosk_score.RDS"), 
    NNGP_heterosk =          readRDS("res/NNGP_heterosk_score.RDS"), 
    NNGP_vanilla =           readRDS("res/NNGP_vanilla_score.RDS"), 
    NNGP_lociso_heterosk =   readRDS("res/NNGP_lociso_heterosk_score.RDS")
  )
)
res = res[,c("train_elpd", "loo_elpd", "lump_elpd", "train_coverage", "loo_coverage", "lump_center_coverage", "loo_MSE", "lump_center_MSE", "time")]
ncol(res)
res = cbind(
  do.call(cbind, list(
    "method" = c("NNGP", "NNGP", "NNGP", "NNGP", "NNGP"),
    "estimation" = c("Bayesian", "Bayesian", "Bayesian", "Bayesian", "Bayesian"),
    "aniso" = c("yes", "yes", "no", "no", "no"),
    "nonstat.range" = c("yes", "no", "no", "no", "yes"),
    "nonstat.noise" = c("yes", "yes", "yes", "no", "yes")
  )), res
)
res = rbind(res, "local GP" =  c("local GP", "local MLE", "no", "yes", "yes", unlist(readRDS(file = "res/local_GPs_score.RDS"))))
res = rbind(res, "INLA stat" = c("INLA", "Bayesian", "no", "no", "no", unlist(readRDS(file = "res/inla_stat_score.RDS")[c("train_elpd", "loo_elpd", "lump_elpd", "train_coverage", "loo_coverage", "lump_center_coverage", "loo_MSE", "lump_center_MSE", "time")])))
res = as.data.frame(res)
for(x in c("train_elpd", "loo_elpd", "lump_elpd", "train_coverage", "lump_center_coverage", "loo_coverage", "loo_MSE", "lump_center_MSE", "time"))res[, x] = as.numeric(res[, x])
res["NNGP_full_monty", "time"] = res[1, "time"]*24
res["INLA stat", "time"] = res["INLA stat", "time"]/60
res[,"time"] = sapply(res[,"time"], function(x)paste(as.integer(x), "h", round(60*(x - as.integer(x))), "min"))
res[,c("train_elpd", "loo_elpd", "lump_elpd", "train_coverage", "loo_coverage", "lump_center_coverage", "loo_MSE", "lump_center_MSE")] = round(res[,c("train_elpd", "loo_elpd", "lump_elpd", "train_coverage", "loo_coverage", "lump_center_coverage", "loo_MSE", "lump_center_MSE")], 2)

res = rbind(res, c("INLA", "empirical Bayes", "no", "no", "yes", NA, NA,  NA, NA,  NA, NA,NA, NA,  "(crashed)"), c("INLA", "empirical Bayes", "no", "yes", "yes", NA, NA, NA, NA, NA, NA, NA, NA, "(crashed)"))

res

res = cbind(res, "min ESS" = c(
            min(unlist(Bidart::ESS(readRDS("res/NNGP_full_monty.RDS")[[1]], burn_in = .1))),
            min(unlist(Bidart::ESS(readRDS("res/NNGP_aniso_heterosk.RDS")[[1]], burn_in = .1))),
            min(unlist(Bidart::ESS(readRDS("res/NNGP_heterosk.RDS")[[1]], burn_in = .1))),
            min(unlist(Bidart::ESS(readRDS("res/NNGP_vanilla.RDS")[[1]], burn_in = .1))),
            min(unlist(Bidart::ESS(readRDS("res/NNGP_lociso_heterosk.RDS")[[1]], burn_in = .1))), 
            NA, NA, NA, NA
))


res[,"min ESS"] = round(res[,"min ESS"])
row.names(res) = NULL

sink("table_recapitulative.txt")
print(res)
sink()
xtable::xtable((res))


