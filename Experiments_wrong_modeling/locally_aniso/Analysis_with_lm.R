library(magrittr)
setwd("Experiments_wrong_modeling/locally_aniso/")
where_results_are = "res/aniso/"
# checking diags for first runs
i = 6+0+0
  
res = readRDS(paste(where_results_are, "res", i, "complete.RDS", sep = ""))
res$inputs$inputs
Bidart::diagnostic_plots(res$run, starting_proportion = .1, burn_in = .5)
Bidart::ESS(res$run, burn_in = .5)
res$performance





inputs = NULL
diags = NULL
cover = NULL
result_name = "res114complete.RDS"
for(result_name in list.files(where_results_are)[grep("complete", list.files(where_results_are))])
{
 print(result_name)
 res = readRDS(paste(where_results_are, result_name, sep = ""))
 
 # correction because of typo in experiment script
 true_range_var = rep(0, 6)
 if(res$inputs$inputs$data_range == "nonstat_elliptic")true_range_var[c(1,4,6)] = 1
 if(res$inputs$inputs$data_range == "nonstat_circular")true_range_var[c(1)] = 1
 
   if(res$inputs$inputs$model_range=="nonstat_elliptic"){
     range_var_quantiles = 
       apply(res$estimation$samples$range_log_scale, 3, Bidart::expmat, simplify = F) %>% 
       lapply(function(x)x[lower.tri(x, T)]) %>% 
       do.call(what = rbind) %>% 
       apply(2, quantile, probs = c(.025, .975))
   }
 if(res$inputs$inputs$model_range=="nonstat_circular"){
   range_var_quantiles = 
     apply(res$estimation$samples$range_log_scale, 3, Bidart::expmat, simplify = F) %>% 
     lapply(function(x)x[lower.tri(x, T)]) %>% 
     do.call(what = rbind) %>% 
     apply(2, quantile, probs = c(.025, .975))
   range_var_quantiles = cbind(range_var_quantiles, matrix(0, 2,5))
 }
 if(res$inputs$inputs$model_range=="stat")range_var_quantiles = matrix(0, 2,6)
 res$coverage$range_var = (true_range_var>= range_var_quantiles[1,])&(true_range_var<= range_var_quantiles[2,])
 
 
 inputs = rbind(inputs, sapply(res$inputs$inputs, as.character))
 diags = rbind(diags, res$performance)
 cover = rbind(cover, unlist(res$coverage))
}

saveRDS(diags, "diags.RDS")
saveRDS(cover, "cover.RDS")
saveRDS(inputs, "inputs.RDS")

#########################################
# Lm analysis of performance indicators #
#########################################

diags = readRDS("diags.RDS")
inputs = readRDS("inputs.RDS")
inputs[inputs=="nonstat"]="N"
inputs[inputs=="stat"]="S"

# color for a model combination


model_color = function(model, dataset)
{
  if(model==dataset)return("color{green}")
  if((model=="nonstat_elliptic")&(dataset=="S"))return("color{blue}")
  if((model=="nonstat_circular")&(dataset=="S"))return("color{blue}")
  if((model=="nonstat_elliptic")&(dataset=="nonstat_circular"))return("color{blue}")
  if((model=="S")&(dataset=="nonstat_circular"))return("color{gray}")
  if((model=="S")&(dataset=="nonstat_elliptic"))return("color{gray}")
  if((model=="nonstat_circular")&(dataset=="nonstat_elliptic"))return("color{gray}")
}

# putting types of models and data in a more simple format
data_types  = inputs[,2]
model_types = inputs[,1]

# split models following data
model_types =split(model_types, data_types)
# convert model to factor
model_types = lapply(model_types, as.factor)
# put data type as reference factor for model type
model_types = mapply(function(x, y)relevel(x, ref = y), model_types, names(model_types), SIMPLIFY = F)

# regress performance of diags on model type, for each data type
name = "smooth_field_mse" 
regression_perf_model = lapply(colnames(diags), function(name)
{
  perf = unlist(diags[,name])
  mapply(function(perf, model_type, seed)lm(perf~model_type + seed), 
         perf = split(perf, data_types), 
         seed = split(inputs[,"seed"], data_types),
         model_types, 
         SIMPLIFY = F)
})
names(regression_perf_model)=colnames(diags)



# put the coefficients in a given order (they do not come in the same order in the regressions)
#smooth_mse = regression_perf_model[[1]]
#coeff_list = smooth_mse
#coeffs =summary(coeff_list[[1]])$coefficients[,1]
put_coeffs_in_order = function(coeffs, model_cases)
{
  sapply(
    model_cases, 
    function(x){
      if(length(grep(pattern = x, names(coeffs)))==0)return(0)
      return(coeffs[grep(x, names(coeffs))])
    })
}

# the 3 model cases
model_cases = c("nonstat_elliptic", "nonstat_circular", "S")

model_case_tex = c("A", "alpha", "emptyset")

# putting regression coeffs together
coeffs = lapply(
  # getting regression coefficients for each type of data and each type of model
  lapply(regression_perf_model, function(x)lapply(x, function(x)summary(x)$coeff[,1])),
  function(x)sapply(x,function(x)put_coeffs_in_order(x, model_cases= unique(data_types)))
)
# reordering
x = coeffs[[1]]
coeffs = lapply(coeffs, function(x)x[sapply(model_cases, function(name)grep(name, row.names(x))), sapply(model_cases, function(name)grep(name, colnames(x)))])
# rounding
coeffs$smooth_field_mse=round(coeffs$smooth_field_mse, 2)
coeffs$pred_field_mse=round(coeffs$pred_field_mse, 2)
coeffs$DIC=round(coeffs$DIC)
coeffs$log_range_det_mse=round(coeffs$log_range_det_mse, 2)
coeffs$log_range_aniso_mse=round(coeffs$log_range_aniso_mse,2 )


# putting p-values coeffs together
p_vals = lapply(
  # getting regression coefficients for each type of data and each type of model
  lapply(regression_perf_model, function(x)lapply(x, function(x)summary(x)$coeff[,4])),
  function(x)sapply(x,function(x)put_coeffs_in_order(x, model_cases= unique(data_types)))
)
# reordering
p_vals = lapply(p_vals, function(x)x[sapply(model_cases, function(name)grep(name, row.names(x))), sapply(model_cases, function(name)grep(name, colnames(x)))])

# converting to little stars
get_p_val_stars = function(p_val_mat)
{
  res = p_val_mat
  res[p_val_mat<1]=""
  res[p_val_mat<.01]="*"
  res[p_val_mat<.001]="**"
  res[p_val_mat< 2e-16]="***"
  res[p_val_mat==0]=""
  res
}
get_p_val_stars(matrix(c(0, 1e-16, .0001, .005, .02, .08, .14)))
p_val_stars  = lapply(p_vals, function(x)get_p_val_stars(x))

# matrix of model-data color code
colors_matrix = matrix(mapply(model_color, expand.grid(model_cases,model_cases)[,1], expand.grid(model_cases,model_cases)[,2]), length(model_cases))

# tex tables
#MODEL IN ROW
tables_to_be_exported = mapply(ascii::paste.matrix, lapply(seq(5), function(x)colors_matrix), coeffs, p_val_stars, SIMPLIFY = F)
names(tables_to_be_exported)=names(p_val_stars)

tables_to_be_exported = lapply(tables_to_be_exported, function(x)
  {
  unname(rbind(c(" ", model_case_tex),cbind(model_case_tex, x)))
}
  )
sink("tables.txt")
for(i in seq(length(tables_to_be_exported)))
{
  cat("SLASHbegin{minipage}{.5SLASHtextwidth}")
  cat("\n")
  print(xtable::xtable(tables_to_be_exported[[i]]))
  cat(paste("SLASHsubcaption{", names(tables_to_be_exported)[i], " cover}", sep=""))
  cat("\n")
  cat(paste("SLASHlabel{tab:", names(tables_to_be_exported)[i], "_aniso_cover}", sep=""))
  cat("\n")
  cat("SLASHend{minipage}") 
  cat("\n")
}
sink()

############################
# plotting range log scale #
############################

i = 3 +9
res = readRDS(paste(where_results_are, "res", i, "complete.RDS", sep = ""))
pdf("range_over_stat.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(150)], dim = c(6, 1, 150))),
  iterations = res$run$iterations$thinning[seq(80)], varname = "range_log_scale", starting_proportion = 0)
dev.off()

i = 6 + 9
res = readRDS(paste(where_results_are, "res", i, "complete.RDS", sep = ""))
pdf("range_over_circ.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(150)], dim = c(6, 1, 150))),
  iterations = res$run$iterations$thinning[seq(80)], varname = "range_log_scale", starting_proportion = 0)
dev.off()

i = 9 + 9
res = readRDS(paste(where_results_are, "res", i, "complete.RDS", sep = ""))
pdf("range_over_ell.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(150)], dim = c(6, 1, 150))),
  iterations = res$run$iterations$thinning[seq(80)], varname = "range_log_scale", starting_proportion = 0)
dev.off()


######################
# analysing coverage #
######################
inputs = readRDS("inputs.RDS")
cover = readRDS("cover.RDS")


cover_split = split(cover, col(cover))
names(cover_split) = c(paste("range", seq(3), sep = ""), "scale", "noise", paste("range_var", seq(6), sep = ""))

cover_tables = lapply(cover_split, function(x)tapply(x, INDEX = list(inputs[,1], inputs[,2]), FUN = mean)[c(2,1,3), c(2,1,3)])
cover_tables = lapply(cover_tables, round, digits = 2)
cover_tables = lapply(cover_tables, function(x)ascii::paste.matrix(colors_matrix, x))
cover_tables = lapply(cover_tables, function(x){
  colnames(x) = c("A", "alpha", "emptyset")
  row.names(x) = c("A", "alpha", "emptyset")
  return(x)
})
sink("tables_cover.txt")
for(i in seq(length(cover_tables))){
  cat("SLASHbegin{minipage}{.3SLASHtextwidth}")
  cat("\n")
  cat("SLASHcentering")
  cat("\n")
  print(xtable::xtable(cover_tables[[i]]))
  cat(paste("SLASHsubcaption{", names(cover_tables)[i], " cover}", sep=""))
  cat("\n")
  cat(paste("SLASHlabel{tab:", names(cover_tables)[i], "_aniso_cover}", sep=""))
  cat("\n")
  cat("SLASHend{minipage}") 
  cat("\n")
  cat("\n")
}
sink()
