setwd("UCLA/range_r3/")
seed = 1
res = array(0, c(4, 49, 9))
namez = names(readRDS(paste("range_r3/range", seed, ".RDS", sep = ""))[[1]])
for(seed in seq(49)){
  experiment = readRDS(paste("range_r3/range", seed, ".RDS", sep = ""))
  for(model in seq(length(experiment))){
    res[model, seed, ] = experiment[[model]][namez]
  }
}
res = res[,which(res[1,,1]!="tatato"),-9]
res = array(as.numeric(res), dim(res))
extracted_info = apply(res, c(1, 3), mean )
row.names(extracted_info) = names(experiment)
print(names(experiment[[1]])[-9])
colnames(extracted_info) = c("smooth MSE", "pred MSE", "log range MSE", "obs cover", "pred cover", "obs elpd", "pred elpd", "time")

extracted_info[,"time"] = round(extracted_info[,"time"])
extracted_info[,-8] = round(extracted_info[,-8], 2)

extracted_info = as.data.frame(extracted_info)
extracted_info$time = paste(extracted_info$time, c("min", "min", "sec", "h"))

xtable::xtable(extracted_info)
