setwd("UCLA/heterosk_r3/")
seed = 1
res = array(0, c(4, 49, 9))
namez = names(readRDS(paste("heterosk_r3/het", 1, ".RDS", sep = ""))[[1]])
for(seed in seq(49)){
  experiment = readRDS(paste("heterosk_r3/het", seed, ".RDS", sep = ""))
  for(model in seq(length(experiment))){
    res[model, seed, ] = experiment[[model]][namez]
  }
}

experiment = readRDS(paste("heterosk_r3/het", 1, ".RDS", sep = ""))
print(paste("there are", sum(is.na(res[1,,1])), "NAs"))
res = res[,which(!is.na(res[1,,1])),-9]
res = array(as.numeric(res), dim(res))
dim(res)
extracted_info = apply(res, c(1, 3), mean )
row.names(extracted_info) = names(experiment)
print(namez)
colnames(extracted_info) = c("smooth MSE", "pred MSE", "noise logvar MSE", "obs coverage", "pred coverage", "obs elpd", "pred elpd", "time")

extracted_info = as.data.frame(extracted_info)

extracted_info[,"time"] = round(extracted_info[,"time"])
extracted_info[,-8] = round(extracted_info[,-8], 2)
extracted_info[,8] = round(extracted_info[,8])
extracted_info[,8] = paste(extracted_info[,8], c("min", "min", "h", "s"))


xtable::xtable(extracted_info)
