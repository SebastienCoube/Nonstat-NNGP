setwd("Experiments_wrong_modeling/locally_iso/")
library(magrittr)
where_results_are = "res/iso/"

#################################
# checking diags for first runs #
#################################
i =
   1+9*64
  # 8 # stat 
  # 16 # nonstat range
  # 24 + 64# nonstat scale, first one crashed 
  #32 # nonstat range +scale
  #40 # nonstat noise
  #56 # nonstat scale + noise
  #64 # nonstat full
  
res = readRDS(paste(where_results_are, "res", i, "complete.RDS", sep = ""))
res$estimation$summaries$range_beta
res$coverage$range_beta
Bidart::diagnostic_plots(res$run, starting_proportion = .2)
res$inputs$inputs

##############
# extracting #
##############

inputs = NULL
diags = NULL
cover = NULL

results = list.files(where_results_are)[grep("complete", list.files(where_results_are))][-seq(3)]
i = 1
for(i in seq(i, length(results)))
{
  result_name = results[i]
  print(result_name)
  res = readRDS(paste(where_results_are, result_name, sep = ""))
  inputs = rbind(inputs, sapply(res$inputs$inputs, as.character))
  diags = rbind(diags, res$performance)
  cover = rbind(cover, res$coverage)
}

saveRDS(diags, "diags.RDS")
saveRDS(inputs, "inputs.RDS")
saveRDS(cover, "cover.RDS")



#########################################
# Lm analysis of performance indicators #
#########################################

diags = readRDS("diags.RDS")
inputs = readRDS("inputs.RDS")
inputs[inputs=="nonstat"]="N"
inputs[inputs=="stat"]="S"

# color for a model combination
correct_model = function(model, dataset)
{
  if((model=="N")&(dataset=="S"))return("over")
  if((model=="S")&(dataset=="N"))return("under")
  return("right")
}
model_color = function(model, dataset)
{
  res = c()
  res[1]= correct_model(substr(model, 1,1), substr(dataset, 1,1))
  res[2]= correct_model(substr(model, 3,3), substr(dataset, 3,3))
  res[3]= correct_model(substr(model, 5,5), substr(dataset, 5,5))
  if(all(res=="right"))return("color{green}")
  if(all(res=="right"|res=="under"))return("color{gray}")
  if(all(res=="right"|res=="over"))return("color{blue}")
  return("color{brown}")
}

model_color("N N N", "N N N")
model_color("N N N", "N N S")
model_color("S N N", "N N N")
model_color("S N N", "N S N")

# putting types of models and data in a more simple format
data_types  = mapply(paste, inputs[,4], inputs[,5], inputs[,6])
model_types = mapply(paste, inputs[,1], inputs[,2], inputs[,3])

# split models following data
model_types =split(model_types, data_types)
# convert model to factor
model_types = lapply(model_types, as.factor)
# put data type as reference factor for model type
model_types = mapply(function(x, y)relevel(x, ref = y), model_types, names(model_types))

# regress performance of diags on model type, for each data type
regression_perf_model = lapply(colnames(diags), function(name)
{
  perf = unlist(diags[,name])
  mapply(function(perf, model_type, seed)lm(perf~model_type + seed), 
         perf = split(perf, data_types), 
         seed = split(inputs[,"seed"], data_types) %>% lapply(unlist),
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

# the 8 model cases
model_cases = expand.grid(c("N", "S"), c("N", "S"), c("N", "S"))
model_cases = split(model_cases, row(model_cases))
model_cases = lapply(model_cases, function(x)do.call(paste, x))


model_cases_tex = as.matrix(expand.grid(c("alpha", ""), c("sigma", ""), c("tau", "")))
model_cases_tex = split(model_cases_tex, row(model_cases_tex))
model_case_tex = lapply(model_cases_tex, function(x)paste(x[1], x[2], x[3]))
model_case_tex[8]="emptyset"
model_case_tex = unlist(model_case_tex)

# putting regression coeffs together
coeffs = lapply(
  # getting regression coefficients for each type of data and each type of model
  lapply(regression_perf_model, function(x)lapply(x, function(x)summary(x)$coeff[,1])),
  function(x)sapply(x,function(x)put_coeffs_in_order(x, model_cases= unique(data_types)))
)
# reordering
coeffs = lapply(coeffs, function(x)x[sapply(model_cases, function(name)grep(name, row.names(x))), sapply(model_cases, function(name)grep(name, colnames(x)))])
# rounding
coeffs = lapply(coeffs, round, digits = 2)
coeffs$DIC = round(coeffs$DIC)

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
get_p_val_stars(matrix(c(0, .0001, .005, .02, .08, .14)))
p_val_stars  = lapply(p_vals, function(x)get_p_val_stars(x))

# matrix of model-data color code
colors_matrix = matrix(mapply(model_color, expand.grid(model_cases,model_cases)[,1], expand.grid(model_cases,model_cases)[,2]), length(model_cases))

# tex tables
#MODEL IN ROW
tables_to_be_exported = mapply(ascii::paste.matrix, list(colors_matrix,colors_matrix,colors_matrix), coeffs, p_val_stars, SIMPLIFY = F)
names(tables_to_be_exported)=names(p_val_stars)

tables_to_be_exported = lapply(tables_to_be_exported, function(x)
  {
  x = cbind(model_case_tex, x)
  x = rbind(c(" ", model_case_tex), x)
  x
}
  )
sink("tables.txt")
cat("{\n")
xtable::xtable(tables_to_be_exported$smooth_field_mse)
cat("\\caption{Latent field smoothing.}\n")
cat("\\label{tab:wrong_modeling_iso_smooth}\n")
cat("}\n")

cat("{\n")
xtable::xtable(tables_to_be_exported$pred_field_mse)
cat("\\caption{Latent field prediction.}\n")
cat("\\label{tab:wrong_modeling_iso_pred}\n")
cat("}\n")


cat("{\n")
xtable::xtable(tables_to_be_exported$DIC)
cat("\\caption{DIC.}\n")
cat("\\label{tab:wrong_modeling_iso_DIC}\n")
cat("}\n")


cat("{\n")
xtable::xtable(tables_to_be_exported$log_range_mse)
cat("\\caption{Log-range smoothing MSE.}\n")
cat("\\label{tab:wrong_modeling_iso_range}\n")
cat("}\n")

cat("{\n")
xtable::xtable(tables_to_be_exported$log_scale_mse)
cat("\\caption{Log-variance of the latent field smoothing MSE.}\n")
cat("\\label{tab:wrong_modeling_iso_scale}\n")
cat("}\n")


cat("{\n")
xtable::xtable(tables_to_be_exported$log_noise_mse)
cat("\\caption{Log-variance of the noise smoothing MSE.}\n")
cat("\\label{tab:wrong_modeling_iso_noise}\n")
cat("}\n")

sink()


######################
# Plotting log scale #
######################

res = readRDS("res/iso/res2complete.RDS")
pdf("range_over.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(150)], dim = c(1, 1, 150))),
  iterations = res$run$iterations$thinning[seq(40)], varname = "range_log_scale", starting_proportion = 0)
dev.off()

res = readRDS("res/iso/res10complete.RDS")
pdf("range_just_right.pdf")
par(mfrow = c(1, 1))
Bidart::plot_log_scale(
  lapply(res$run$records, function(x)array(x$range_log_scale[,,seq(150)], dim = c(1, 1, 150))),
  iterations = res$run$iterations$thinning[seq(40)], varname = "range_log_scale", starting_proportion = 0)
dev.off()



######################
# analysing coverage #
######################

inputs = readRDS("inputs.RDS")
inputs[inputs=="nonstat"]="N"
inputs[inputs=="stat"]="S"

cover = readRDS("cover.RDS")%>% unlist() %>% matrix(nrow= nrow(readRDS("inputs.RDS")))
colnames(cover) = colnames(readRDS("cover.RDS"))

inputs_split = split(inputs, col(inputs))
names(inputs_split) = colnames(inputs)
cover_split = split(cover, col(cover))
names(cover_split) = colnames(cover)

data_types  = mapply(paste, inputs[,4], inputs[,5], inputs[,6])
model_types = mapply(paste, inputs[,1], inputs[,2], inputs[,3])


cover_tables = lapply(cover_split, function(x)tapply(x, INDEX = list(unlist(model_types), data_types), FUN = mean))
cover_tables = lapply(cover_tables, round, digits = 2)
cover_tables = lapply(cover_tables, function(x)x[sapply(model_cases, function(name)grep(name, row.names(x))), sapply(model_cases, function(name)grep(name, colnames(x)))])

cover_tables = lapply(cover_tables, function(x)ascii::paste.matrix(colors_matrix, x))
cover_tables = lapply(cover_tables, function(x)
{
  x = cbind(model_case_tex, x)
  x = rbind(c(" ", model_case_tex), x)
  x = cbind(" ", x)
  x = unname(x)
  x
}
)


sink("tables_cover.txt")
for(i in seq(length(cover_tables))){
  cat("{\n")
  print(xtable::xtable(cover_tables[[i]]), include.rownames=FALSE, include.colnames=FALSE)
  cat(paste("SLASHsubcaption{", names(cover_tables)[i], " cover}", sep=""))
  cat(paste("SLASHlabel{tab:", names(cover_tables)[i], "_iso_cover}", sep=""))
  cat("}\n")
}
sink()





log_scales = NULL

results = list.files(where_results_are)[grep("complete", list.files(where_results_are))][-seq(3)]

for(result_name in results)
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
        xlab = "Model used to generate the data", ylab = expression(paste("estimate of  ", gamma, "_", alpha)), names =data_names)
dev.off()
pdf("Experiments_wrong_modeling/locally_iso/sigma.pdf")
boxplot(log_scales[,2]~log_scales[,3] + log_scales[,4], #main = expression(paste("log scale of ", sigma^2)), 
        xlab = "Model used to generate the data", names =data_names, 
        ylab = expression(paste("estimate of  ", gamma, "_", sigma))
)
dev.off()




