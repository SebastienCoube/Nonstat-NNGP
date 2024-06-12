setwd("UCLA/range/")
seed = 1
res = array(0, c(4, 49, 5))
for(seed in seq(49)){
  experiment = readRDS(paste("range/range", seed, ".RDS", sep = ""))
  for(model in seq(length(experiment))){
    res[model, seed, ] = experiment[[model]]
  }
}
res = res[,which(res[1,,1]!="tatato"),-5]
res = array(as.numeric(res), dim(res))
extracted_info = apply(res, c(1, 3), mean )
row.names(extracted_info) = names(experiment)
print(names(experiment[[1]])[-5])
colnames(extracted_info) = c("smooth MSE", "pred MSE", "log range MSE", "time")

extracted_info[,"time"] = round(extracted_info[,"time"])
extracted_info[,-4] = round(extracted_info[,-4], 3)

xtable::xtable(extracted_info)
