
# checking diags for first runs
i = 56
res = readRDS(paste("Experiments_wrong_modeling/locally_iso/res", i, "complete.RDS", sep = ""))
Bidart::diagnostic_plots(res$run)


where_results_are = "Experiments_wrong_modeling/locally_iso/"

inputs = NULL
diags = NULL

for(result_name in list.files(where_results_are)[grep("complete", list.files(where_results_are))])
{
  print(result_name)
  res = readRDS(paste(where_results_are, result_name, sep = ""))
  inputs = rbind(inputs, sapply(res$inputs$inputs, as.character))
  diags = rbind(diags, res$performance)
}

saveRDS(diags, "Experiments_wrong_modeling/locally_iso/diags.RDS")
saveRDS(inputs, "Experiments_wrong_modeling/locally_iso/inputs.RDS")

diags = readRDS("Experiments_wrong_modeling/locally_iso/diags.RDS")
inputs = readRDS("Experiments_wrong_modeling/locally_iso/inputs.RDS")


model_color = function(model, dataset)
{
  if(all(model==dataset))return(3)
  if(all(model==dataset | (model == "S" & dataset == "N")))return(8)
  if(all(model==dataset | (model == "N" & dataset == "S")))return(4)
  return(2)
}

colnames(inputs)
model_cases = expand.grid(c("N", "S"), c("N", "S"), c("N", "S"), stringsAsFactors = F)
model_cases =  split(model_cases, row(model_cases))
inputs[inputs=="nonstat"]="N"
inputs[inputs=="stat"]="S"

model_name = function(range_model, scale_model, noise_model)
{
  if(range_model=="N"&scale_model=="N"&noise_model=="N")return(expression(paste(alpha, "+", sigma^2, "+", tau^2)))
  if(range_model=="S"&scale_model=="N"&noise_model=="N")return(expression(paste(sigma^2, "+", tau^2)))
  if(range_model=="S"&scale_model=="S"&noise_model=="N")return(expression(paste(tau^2)))
  if(range_model=="S"&scale_model=="N"&noise_model=="S")return(expression(paste(sigma^2)))
  if(range_model=="N"&scale_model=="S"&noise_model=="N")return(expression(paste(alpha, "+", tau^2)))
  if(range_model=="N"&scale_model=="S"&noise_model=="S")return(expression(paste(alpha)))
  if(range_model=="N"&scale_model=="N"&noise_model=="S")return(expression(paste(alpha, "+", sigma^2)))
  if(range_model=="S"&scale_model=="S"&noise_model=="S")return(expression(symbol("\306")))
}

model_names_ = sapply(model_cases, function(x)model_name(x[1],  x[2], x[3]))

# data type
dev.off()
for(i in unique(inputs[,4]))
{ 
for(j in unique(inputs[,5]))
{ 
for(k in unique(inputs[,6]))
{ 
  for (diag in c("smooth_field_mse", "pred_field_mse", "DIC" ))
  {
    print(i)
    print(j)
    print(k)
    pdf(paste("Experiments_wrong_modeling/locally_iso/", diag, "_", i, "_", j, "_", k, ".pdf", sep = ""), height = 4, width = 6)
    idx = which(inputs[,4]==i & inputs[,5]==j & inputs[,6]==k)
    color = sapply(model_cases, function(model)model_color(model, c(i, j, k)))
    boxplot(unlist(diags[idx,diag]) ~inputs[idx,1]+inputs[idx,2]+inputs[idx,3], main = ""#model_name(i, j, k)
            ,
            col = color, xlab = "model", ylab = diag, names = model_names_)
    #legend("topleft", fill = c(8, 2, 3, 4), legend = c("under", "wrong", "right", "over"))
    dev.off()
  }
  }
}
}
dev.off()
  
# Illustrating behavior of chains wirh over modeling


res = readRDS("Experiments_wrong_modeling/locally_iso/res2complete.RDS")
pdf("Experiments_wrong_modeling/locally_iso/range_over.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(200)], dim = c(1, 1, 200))),
  iterations = res$run$iterations$thinning[seq(200)], varname = "range_log_scale", starting_proportion = 0)
dev.off()

res = readRDS("Experiments_wrong_modeling/locally_iso/res10complete.RDS")
pdf("Experiments_wrong_modeling/locally_iso/range_just_right.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(200)], dim = c(1, 1, 200))),
  iterations = res$run$iterations$thinning[seq(200)], varname = "range_log_scale", starting_proportion = 0)
dev.off()




# range and scale log scales

where_results_are = "Experiments_wrong_modeling/locally_iso/"

log_scales = NULL

for(result_name in list.files(where_results_are)[grep("complete", list.files(where_results_are))])
{
  print(result_name)
  res = readRDS(paste(where_results_are, result_name, sep = ""))
  if((res$inputs$inputs$model_range == "nonstat")&(res$inputs$inputs$model_scale == "nonstat")
     &(res$inputs$inputs$data_noise == "stat"))
  log_scales = rbind(log_scales, c(res$estimation$summaries$range_log_scale[1], res$estimation$summaries$scale_log_scale[1], res$inputs$inputs$data_range, res$inputs$inputs$data_scale))
}
colnames(log_scales) = c("mean_range_log_scale", "mean_scale_log_scale", "range_data", "scale_data")
saveRDS(object = log_scales, "Experiments_wrong_modeling/log_scales.RDS")


log_scales = readRDS("Experiments_wrong_modeling/log_scales.RDS")
data_names = c(expression(symbol("\306")), expression(alpha), expression(sigma^2), expression(paste(alpha, "+", sigma^2)))
  
pdf("Experiments_wrong_modeling/locally_iso/alpha.pdf")
boxplot(log_scales[,1]~log_scales[,3] + log_scales[,4], #main = expression(paste("log scale of ", alpha)), 
        xlab = "data type", ylab = "log-NNGP log variance estimate", names =data_names)
dev.off()
pdf("Experiments_wrong_modeling/locally_iso/sigma.pdf")
boxplot(log_scales[,2]~log_scales[,3] + log_scales[,4], #main = expression(paste("log scale of ", sigma^2)), 
        xlab = "data type", ylab = "log-NNGP log variance estimate", names =data_names)
dev.off()



