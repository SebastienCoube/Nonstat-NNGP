
where_results_are = "Experiments_wrong_modeling/locally_aniso/"

inputs = NULL
diags = NULL

for(result_name in list.files(where_results_are)[grep("complete", list.files(where_results_are))])
{
  print(result_name)
  res = readRDS(paste(where_results_are, result_name, sep = ""))
  #if(!is.null(res$run))
  #{
  #  inputs = rbind(inputs, sapply(res$inputs$inputs, as.character))
  #  diags = rbind(diags, res$performance)
  #  diags[nrow(diags),"DIC"] = DIC(res$run)
  #  
  #}
  inputs = rbind(inputs, sapply(res$inputs$inputs, as.character))
  diags = rbind(diags, res$performance)
}

#saveRDS(diags, "Experiments_wrong_modeling/locally_aniso/diags.RDS")
#saveRDS(inputs, "Experiments_wrong_modeling/locally_aniso/inputs.RDS")

diags = readRDS("Experiments_wrong_modeling/locally_aniso/diags.RDS")
inputs = readRDS("Experiments_wrong_modeling/locally_aniso/inputs.RDS")
inputs[inputs == "nonstat_elliptic"] = "nonstat_aniso"

model_color = function(model, dataset)
{
  if(model==dataset)return(3)
  if((model=="nonstat_aniso")&(dataset=="stat"))return(4)
  if((model=="nonstat_circular")&(dataset=="stat"))return(4)
  if((model=="nonstat_aniso")&(dataset=="nonstat_circular"))return(4)
  if((model=="stat")&(dataset=="nonstat_circular"))return(8)
  if((model=="stat")&(dataset=="nonstat_aniso"))return(8)
  if((model=="nonstat_circular")&(dataset=="nonstat_aniso"))return(8)
}


model_name = function(range_model)
{
  if(range_model=="nonstat_aniso")return(expression(paste(A)))
  if(range_model=="nonstat_circular")return(expression(paste(alpha)))
  return(expression(symbol("\306")))
}

model_cases = c("nonstat_aniso", "nonstat_circular", "stat")
model_names_ = sapply(model_cases, model_name)

# data type
dev.off()
for(i in unique(inputs[,2]))
{ 
  print(i)
  for (diag in c("smooth_field_mse", "pred_field_mse", "DIC" ))
  {
    pdf(paste("Experiments_wrong_modeling/locally_aniso/", diag, "_", i, ".pdf", sep = ""), height = 4, width = 6)
    idx = which(inputs[,2]==i)
    color = sapply(model_cases, function(model)model_color(model, i))
    boxplot(unlist(diags[idx,diag]) ~inputs[idx,1], main = ""#model_name(i)
              ,
            col = color, xlab = "model", ylab = diag, names = model_names_
            )
    dev.off()
  }
}
dev.off()

# under-modelling, right modelling

res = readRDS("Experiments_wrong_modeling/locally_aniso/res3complete.RDS")
pdf("Experiments_wrong_modeling/locally_aniso/range_over_stat.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(200)], dim = c(6, 1, 200))),
  iterations = res$run$iterations$thinning[seq(200)], varname = "range_log_scale", starting_proportion = 0)
dev.off()

res = readRDS("Experiments_wrong_modeling/locally_aniso/res6complete.RDS")
pdf("Experiments_wrong_modeling/locally_aniso/range_over_circ.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(200)], dim = c(6, 1, 200))),
  iterations = res$run$iterations$thinning[seq(200)], varname = "range_log_scale", starting_proportion = 0)
dev.off()

res = readRDS("Experiments_wrong_modeling/locally_aniso/res9complete.RDS")
pdf("Experiments_wrong_modeling/locally_aniso/range_over_ell.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(200)], dim = c(6, 1, 200))),
  iterations = res$run$iterations$thinning[seq(200)], varname = "range_log_scale", starting_proportion = 0)
dev.off()



