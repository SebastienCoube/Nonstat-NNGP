# #' @export
#plot_ellipses = function(locs, log_range, shrink = .1, main = "ellipses", add  =F)
#{
#  if(ncol(log_range) ==3)matrices = lapply(split(log_range, row(log_range)), Bidart::expmat)
#  if(ncol(log_range) ==1)matrices = lapply(log_range, function(x)exp(x) * diag(1/sqrt(ncol(locs)), ncol(locs)))
#  if(!add)plot(locs, type = "n", xlab = "", ylab = "", main = main)
#  for(i in seq(nrow(locs)))
#  {
#    ell = ellipse::ellipse(matrices[[i]])[seq(25 * 4),] * shrink
#    ell[,1] = ell[,1]+locs[i, 1]
#    ell[,2] = ell[,2]+locs[i, 2]
#    lines(ell)
#  }
#}

#' @export
get_colors = function(x)heat.colors(100)[round((x - min(x))/(max(x)-min(x))*90)+1]
plot_pointillist_painting = function(locs, field, cex = 1, main = NULL)
{
  plot(locs, col = get_colors(field), main = main, pch = 15, cex = cex, xlab  ="", ylab = "")
}

#' @export
array_cumsum = function(x)
{
  res = array(0, dim = dim(x))
  for(i in seq(dim(x)[1]))
  {
    for(j in seq(dim(x)[2]))
    {
      res[i, j, ] = cumsum(x[i, j,])
    }
  }
  res
}

#' @export
array_multiply_3 = function(x, M)
{
  res = array(0, dim = c(dim(x)[c(1, 2)], M@Dim[2]))
  for(i in seq(dim(res)[1]))
  {
    for(j in seq(dim(res)[2]))
    {
      res[i, j, ] = as.vector(x[i, j,] %*% M)
    }
  }
  res
}

#' @export
array_multiply_2 = function(x, M)
{
  res = array(0, dim = c(dim(x)[1], dim(M)[2], dim(x)[3]))
  for(i in seq(dim(res)[1]))
  {
    for(j in seq(dim(res)[3]))
    {
      res[i, ,j ] = as.vector(M %*% x[i, ,j])
    }
  }
  res
}

#' @export
array_multiply_1 = function(x, M)
{
  res = array(0, dim = c(dim(M)[1], dim(x)[2], dim(x)[3]))
  for(i in seq(dim(res)[2]))
  {
    for(j in seq(dim(res)[3]))
    {
      res[,i ,j ] = as.vector(M %*% x[,i,j])
    }
  }
  res
}

#' @export
grb_diags_field = function(record_arrays, iterations, burn_in = .5, starting_proportion = .5)
{
  cumsums = lapply(record_arrays, array_cumsum)
  sqcumsums = lapply(record_arrays, function(x)array_cumsum(x^2))
  iter_start_idx = which(iterations > iterations[length(iterations)] * starting_proportion)[1]
  diff_array = cbind(0, seq(iter_start_idx, length(iterations))) 
  n = rep(0, nrow(diff_array))
  lower_bound_idx = 1
  for(i in seq(nrow(diff_array)))
  {
    while(iterations[lower_bound_idx] < iterations[diff_array[i, 2]]* burn_in){lower_bound_idx = lower_bound_idx+1}
    diff_array[i, 1] = lower_bound_idx
    n[i] = i + iter_start_idx -lower_bound_idx -1
  }
  diff_matrix = Matrix::sparseMatrix(
    j = c(seq(nrow(diff_array)), seq(nrow(diff_array))), 
    i = c(diff_array), 
    x = rep(c(-1, 1), each = nrow(diff_array)), 
    )
  mean_estimators = lapply(cumsums, function(x)array_multiply_3(x, M = diff_matrix %*% Matrix::Diagonal(x = 1/n)))
  sqmean_estimators = lapply(sqcumsums, function(x)array_multiply_3(x, M = diff_matrix %*% Matrix::Diagonal(x = 1/n)))
  var_estimators = mapply(function(x, y)  array_multiply_3(x - y^2, M = Matrix::Diagonal(x = n/(n-1))), x = sqmean_estimators,  y = mean_estimators, SIMPLIFY = F)
  within_mean = Reduce("+", mean_estimators)/length(mean_estimators)
  within_var = Reduce("+", var_estimators)/length(var_estimators)
  between_var = (Reduce("+", lapply(mean_estimators, function(x)x^2))/length(mean_estimators) - within_mean^2) * length(mean_estimators) / (length(mean_estimators)-1)
  #PSRF =   (length(var_estimators) +1)/(length(var_estimators)) * array_multiply_3(x = between_var/within_var, M = Matrix::Diagonal(x = 1/n))
  PSRF =   (length(var_estimators) +1)/(length(var_estimators)) * array_multiply_3(x = between_var/within_var, M = Matrix::Diagonal(x = (n-1)/n))
  for(i in seq(dim(PSRF)[3]))
  {
    PSRF[,,i] = PSRF[,,i] + (n[i]+1)/n[i]
  }
  PSRF_quantiles  =apply(PSRF, 3, function(x)quantile(x, probs = c(1, .99, .9, .5), na.rm = T))
  #PSRF_binned = array(cut(PSRF, breaks = c(-Inf, 1.1, 1.2, 1.5, +Inf)), dim(PSRF))
  list("iterations" = iterations[diff_array[,2]], "mean" = mean_estimators, "var" = var_estimators,"PSRF" = PSRF, "PSRF_quantiles" = PSRF_quantiles#, "PSRF_binned" = PSRF_binned
       )
}


#arrays =                   list(
#  array(rnorm(10000000)+10, dim = c(10, 10, 100000)),
#  array(rnorm(10000000)+10, dim = c(10, 10, 100000)),
#  array(rnorm(10000000)+10, dim = c(10, 10, 100000))
#)
#test = grb_diags_field(iterations = seq(100000)*10, record_arrays = arrays, burn_in = .5, starting_proportion = .5)



#' @export
plot_PSRF = function(PSRF, individual_varnames = NULL, varname = "")
{
  if(any(is.infinite(PSRF$PSRF)))
  {
    plot(0, 0, type = "n", main = paste("PSRF of", varname, "not represented because of infinite GRB diags"), xlab = "iterations", ylab = "PSRF quantiles")
    return()
  }
  if((prod(dim(PSRF$PSRF)[c(1, 2)]))>9)
  {
    plot(PSRF$iterations, PSRF$PSRF_quantiles[1,], ylim = c(1, max(PSRF$PSRF_quantiles[1,])), main = paste("Quantiles of PSRF of", varname), type = "l", xlab = "iterations", ylab = "PSRF quantiles")
    lines(PSRF$iterations, PSRF$PSRF_quantiles[2,], col = 2)
    lines(PSRF$iterations, PSRF$PSRF_quantiles[3,], col = 3)
    lines(PSRF$iterations, PSRF$PSRF_quantiles[4,], col = 4)
    legend(bg = "white", "topleft", legend = c("max", .99, .9, .5), fill = c(1, 2, 3, 4))
  }
  if(prod(dim(PSRF$PSRF)[c(1, 2)])<10)
  {
    plot(PSRF$iterations, PSRF$iterations, type = "n", ylim = c(1, max(PSRF$PSRF)), main = paste("PSRF of each component of", varname), xlab = "iterations", ylab = "PSRF")
    for(i in seq(1, dim(PSRF$PSRF)[1]))
    {
      for(j in seq(1, dim(PSRF$PSRF)[2]))
      {
        lines(PSRF$iterations, PSRF$PSRF[i,j,], col = i+(j-1) * dim(PSRF$PSRF)[1])
      }
    }
    if(dim(PSRF$PSRF)[2]>1 & !is.null(individual_varnames)) labs = paste(rep(individual_varnames, dim(PSRF$PSRF)[2]), rep(seq(dim(PSRF$PSRF)[2]), each = length(individual_varnames)))
    if(dim(PSRF$PSRF)[2]==1 & !is.null(individual_varnames)) labs = individual_varnames
    if(is.null(individual_varnames)) labs = paste("component", seq(prod(dim(PSRF$PSRF)[c(1, 2)])))
    legend(bg = "white", "topleft", legend = labs, fill = seq(length(labs)))
  }
  abline(h = 1.2, lty = "longdash")
  abline(h = 1)
}
#arrays =                   list(
#  array(rnorm(600), dim = c(4, 2, 100)),
#  array(rnorm(600), dim = c(4, 2, 100)),
#  array(rnorm(600), dim = c(4, 2, 100))
#)
#test = grb_diags_field(iterations = seq(100), record_arrays = arrays, burn_in = .5, starting_proportion = .5)
#plot_PSRF(test, individual_varnames = seq(4), varname = "tatato")



#' @export
plot_log_scale = function(log_scale_arrays, iterations, starting_proportion = .5, varname)
{
  kept_iterations = which(iterations>starting_proportion*iterations[length(iterations)])
  if(dim(log_scale_arrays[[1]])[1]==1 & dim(log_scale_arrays[[1]])[2]==1)
  {
    log_scale = sapply(log_scale_arrays, function(x)x[,,kept_iterations])
    plot(iterations[kept_iterations], iterations[kept_iterations], 
         type = "n", xlab = "iteration", ylab = "log scale", main = varname, 
         ylim = c(min(unlist(log_scale)), max(unlist(log_scale))))
    # loop over chains
    for(i in seq(ncol(log_scale)))
    {
      lines(iterations[kept_iterations], log_scale[,i])
    }
  }
  if(dim(log_scale_arrays[[1]])[1]==6)
  {
    log_scale = sapply(log_scale_arrays, function(x)x[,,kept_iterations])
    #marginal_logvars = lapply(log_scale_arrays, function(x)
    #  apply(x[,,kept_iterations], 2, function(x) diag(
    #    rbind(c(1/sqrt(2), 1/sqrt(2), 0), c(1/sqrt(2), -1/sqrt(2), 0), c(0, 0, 1)) %*%
    #      Bidart::symmat(x) %*%
    #      cbind(c(1/sqrt(2), 1/sqrt(2), 0), c(1/sqrt(2), -1/sqrt(2), 0), c(0, 0, 1)) 
    #    ))
    #)
    marginal_logvars = lapply(log_scale_arrays, function(x)
      apply(x[,,kept_iterations], 2, function(x) diag(
          Bidart::symmat(x)
          )
    ))
    plot(iterations[kept_iterations], iterations[kept_iterations], 
         type = "n", xlab = "iteration", ylab = "log scale", main = paste( varname, "split by component"), 
         ylim = c(min(unlist(marginal_logvars)), max(unlist(marginal_logvars))))
    for(i in seq(length(marginal_logvars)))
    {
      for(j in seq(nrow(marginal_logvars[[i]])))
      {
        lines(iterations[kept_iterations],marginal_logvars[[i]][j,], col = j)
      }
    }
    legend(bg = "white", "topleft", legend = c("Determinant", "Anisotropy", "Anisotropy"), fill = seq(3))
  }
}

#arrays =                   list(
#  array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(6, 1, 1000)),
#  array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(6, 1, 1000)),
#  array(rnorm(6000)+ rep(c(10,10,0,10,0,0), 1000), dim = c(6, 1, 1000))
#)
#
#plot_log_scale(log_scale_arrays = arrays, iterations = seq(1000), starting_proportion = .5)




#' @export
plot_beta = function(beta_arrays, iterations, starting_proportion = .5, varname, var_names = NULL)
{
  #if(dim(beta_arrays[[1]])[2]==3)beta_arrays = lapply(beta_arrays, function(x)array_multiply_2(x, matrix(c(1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2), -1/sqrt(2), 0, 0, 0, 1), 3)))
  kept_iterations = which(iterations>starting_proportion*iterations[length(iterations)])
  for(i in seq(dim(beta_arrays[[1]])[2]))
  {
    upper_window = max(sapply(beta_arrays, function(x) max(x[,i,kept_iterations])))
    lower_window = min(sapply(beta_arrays, function(x) min(x[,i,kept_iterations])))
    if(dim(beta_arrays[[1]])[2] == 3)
    {
      if(i ==1)main = paste("Determinant component of", varname)
      if(i !=1)main = paste("anisotropy component", i-1, "of", varname)
    }else main = varname
    plot(iterations[kept_iterations], iterations[kept_iterations], type = "n", main = main, 
         ylab = "parameter value", xlab = "iteration", ylim = c(lower_window, upper_window)
         )
    for(j in seq(dim(beta_arrays[[1]])[1]))
    {
      for(x in beta_arrays)lines(iterations[kept_iterations], x[j, i, kept_iterations], col = j)
    }
    if(!is.null(var_names))legend(bg = "white", "topleft", legend = var_names, fill = seq(dim(beta_arrays[[1]])[1]))
  }
}

beta_arrays = list(
  array(rnorm(400), dim = c(2, 2, 100)),
  array(rnorm(400), dim = c(2, 2, 100)),
  array(rnorm(400), dim = c(2, 2, 100))
)

#plot_beta(beta_arrays, seq(100), starting_proportion = .5, varname = "tatato", var_names = c(1, 2))
#' @export
diagnostic_plots = function(mcmc_nngp_list, plot_PSRF_fields = F, plot_PSRF_actual_parameters = F, burn_in = .5, starting_proportion = .5)
{
  records_names = names(mcmc_nngp_list$records$chain_1)
  if(!plot_PSRF_fields)records_names = records_names[-grep("field", records_names)]
  range_names = records_names[grep("range", records_names)]
  records_names = setdiff(records_names, range_names)
  noise_names = records_names[grep("noise", records_names)]
  records_names = setdiff(records_names, noise_names)
  scale_names = records_names[grep("scale", records_names)]
  records_names = setdiff(records_names, scale_names)
  response_names = records_names
  for(i in list(range_names, noise_names, scale_names, response_names))
  {
    par(mfrow = c(1, 1))
    for(j in i)
    {
      PSRF = grb_diags_field(record_arrays = lapply(mcmc_nngp_list$records, function(x)x[[j]]), iterations = mcmc_nngp_list$iterations$thinning, burn_in = burn_in, starting_proportion = starting_proportion)
      plot_PSRF(PSRF = PSRF, varname = j, individual_varnames = row.names(mcmc_nngp_list$states$chain_1$params[[j]]))
    }
  }
  for(i in list(range_names, noise_names, scale_names))
  {
    #print(i)
    #par(mfrow = c(length(grep("log_scale", i)) + dim(mcmc_nngp_list$records$chain_1 [[i[grep("beta", i)]]])[2], 1))
    par(mfrow = c(1, 1))
    if(length(grep("log_scale", i))>0)plot_log_scale(log_scale_arrays = lapply(mcmc_nngp_list$records, function(x)x[[i[grep("log_scale", i)]]]), iterations = mcmc_nngp_list$iterations$thinning, starting_proportion = starting_proportion, varname = i[grep("log_scale", i)])
    plot_beta(beta_arrays = lapply(mcmc_nngp_list$records, function(x)x[[i[grep("beta", i)]]]), iterations = mcmc_nngp_list$iterations$thinning, starting_proportion = starting_proportion, varname = i[grep("beta", i)], var_names = row.names(mcmc_nngp_list$states$chain_1$params[[i[grep("beta", i)]]]))
  }
}
