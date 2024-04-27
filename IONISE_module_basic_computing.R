library(MASS)
library(mvtnorm)
library(dplyr)
library(invgamma)

get_divisor_dist <- function(theta, dist_type){
  if (dist_type == "lognormal"){
    dd = ceiling(5*max(theta[2], theta[4]))
    #dd = ceiling(5/min(theta[2], theta[4]))
  }else{
    if (theta[1] != 1 & theta[3] != 1){
      dd = ceiling(5 / min((theta[1]-1)*theta[2], (theta[3]-1)*theta[4]))
    }else if(theta[1] != 1){
      dd = ceiling(5 / min((theta[1]-1)*theta[2], theta[4]))
    }else if(theta[3] != 1){
      dd = ceiling(5 / min(theta[2], (theta[3]-1)*theta[4]))
    }else{
      dd = ceiling(5 / min(theta[2], theta[4]))
    }
  }
  dd = dd * 5
  return(dd)
}

eval_density <- function(x, theta, dist_type, density_type){
  if (density_type == "pdf"){
    if (dist_type == "exp"){
      density_val = dexp(x, rate = 1/theta[2])
    }else if (dist_type == "gamma"){
      density_val = dgamma(x, shape = theta[1], scale = theta[2])
    }else if (dist_type == "invgamma"){
      zero_pos = which(x==0)
      density_val = dinvgamma(x, shape = theta[1], rate = theta[2])
      density_val[zero_pos] = 0
    }else if (dist_type == "lognormal"){
      density_val = dlnorm(x, meanlog = theta[1], sdlog = theta[2])
    }else if (dist_type == "weibull"){
      density_val = dweibull(x, shape = theta[1], scale = theta[2])
    }else{
      warning("The type of distribution must be one of c(exp, gamma, invgamma, lognormal, weibull)!")
      return(NA)
    }
  }else if (density_type == "cdf"){
    if (dist_type == "exp"){
      density_val = pexp(x, rate = 1/theta[2])
    }else if (dist_type == "gamma"){
      density_val = pgamma(x, shape = theta[1], scale = theta[2])
    }else if (dist_type == "invgamma"){
      density_val = pinvgamma(x, shape = theta[1], rate = theta[2])
    }else if (dist_type == "lognormal"){
      density_val = plnorm(x, meanlog = theta[1], sdlog = theta[2])
    }else if (dist_type == "weibull"){
      density_val = pweibull(x, shape = theta[1], scale = theta[2])
    }else{
      warning("The type of distribution must be one of c(exp, gamma, invgamma, lognormal, weibull)!")
      return(NA)
    }
  }else{
    warning("The type of density must be one of c(pdf, cdf)!")
    return(NA)
  }
  return(density_val)
}


numer_int <- function(xspan, yval){
  # xspan has to be evenly spaced.
  if(length(xspan) == 1){
    intval0 = 0
  }else{
    h = xspan[2] - xspan[1]
    # if the number of data points is even then it uses the trapezoidal rule.
    if(length(xspan)%%2 == 0){
      intval0 = 0
      for(jj in 2:length(xspan)){
        intval0 = intval0 + (xspan[jj] - xspan[jj-1]) * (yval[jj]+yval[jj-1])/2
      }
    }else{ # if the number of data points is odd then it uses Simpson's rule.
      intval0 = 0
      for(jj in seq(from=2,by=2,to=length(xspan))){
        intval0 = intval0 + h * (yval[jj-1] + 4*yval[jj] + yval[jj+1])/3
      }
    }
  }
  return(intval0)
}

conv_two_gamma <- function(tmax, alpha1, beta1, alpha2, beta2){
  if(tmax == 0){
    conv_val = 0
    return(conv_val)
  }else{
    xspan = seq(from = 0, to = tmax, by = tmax/1000)
    conv_val = numer_int(xspan, dgamma(xspan, shape = alpha1, scale = beta1) * dgamma(max(xspan) - xspan, shape = alpha2, scale = beta2))
  }
  return(conv_val)
}

conv_two_pdfs <- function(tmax, theta1, theta2, dist_type){
  if(tmax == 0){
    conv_val = 0
    return(conv_val)
  }else{
    xspan = seq(from = 0, to = tmax, by = tmax/1000)
    conv_val = numer_int(xspan, eval_density(xspan, theta = theta1, dist_type = dist_type, density_type = "pdf") * eval_density(max(xspan) - xspan, theta = theta2, dist_type = dist_type, density_type = "pdf"))
  }
  return(conv_val)
}

conv_gam_cdf_pdf <- function(tmax, alpha1, beta1, alpha2, beta2){
  if(tmax == 0){
    conv_val = 0
    return(conv_val)
  }else{
    xspan = seq(from = 0, to = tmax, by = tmax/1000)
    conv_val = numer_int(xspan, pgamma(xspan, shape = alpha1, scale = beta1, lower.tail = FALSE) * dgamma(max(xspan) - xspan, shape = alpha2, scale = beta2))
  }
  return(conv_val)
}

conv_cdf_pdf <- function(tmax, theta1, theta2, dist_type){
  if(tmax == 0){
    conv_val = 0
    return(conv_val)
  }else{
    xspan = seq(from = 0, to = tmax, by = tmax/1000)
    conv_val = numer_int(xspan, eval_density(xspan, theta = theta1, dist_type = dist_type, density_type = "cdf") * eval_density(max(xspan) - xspan, theta = theta2, dist_type = dist_type, density_type = "pdf"))
  }
  return(conv_val)
}

diff_midpoint_extension <- function(x, t_interv){
  # x must be an evenly-spaced time series data.
  if(length(x) == 1){
    return(0)
  }else if(length(x) == 2){
    return(rep(diff(x), 2))
  }else{
    diff_x = diff(x)/t_interv
    dL = length(diff_x)
    deriv_x = rep(NA, length(x))
    deriv_x[2:dL] = (diff_x[1:(dL-1)] + diff_x[2:dL])/2
    deriv_x[1] = 2*diff_x[1] - deriv_x[2]
    deriv_x[dL+1] = 2*diff_x[dL] - deriv_x[dL]
    return(deriv_x)
  }
}



mean_var_to_param <- function(mean, var, dist_type){
  if(dist_type == "gamma"){
    par1 = mean^2 / var
    par2 = var / mean 
  }else if(dist_type == "invgamma"){
    par1 = mean^2 / var + 2
    par2 = (par1-1) * mean
  }else if(dist_type == "lognormal"){
    par2 = sqrt(log(var/mean^2 +1))
    par1 = log(mean) - 1/2 * par2^2
  }else if(dist_type == "weibull"){
    if(1 + var/mean^2 > 30){
      warning("Cannot find the weibull parameters")
    }else{
      aspan = seq(from = 0.3, to = 100, by = 0.01)
      tmp_list = gamma(1+2/aspan)/gamma(1+1/aspan)^2
      par1 = aspan[max(which(tmp_list > 1+var/mean^2))]
      par2 = mean / gamma(1+1/par1)
    }
  }else if(dist_type == "exp"){
    par1 = 1
    par2 = mean
    warning("Due to the lack of DOF, only use the given mean.")
  }
  return(list(par1=par1, par2=par2))
}

param_to_mean_var <- function(par1, par2, dist_type){
  if(dist_type == "gamma"){
    mean1 = par1*par2
    var1 = par1*par2^2
  }else if(dist_type == "invgamma"){
    mean1 = par2 / (par1 - 1)
    var1 = par2^2 / ((par1 - 1)^2 * (par1 - 1))
  }else if(dist_type == "lognormal"){
    mean1 = exp(par1 + par2^2/2)
    var1 = (exp(par2^2) - 1) * exp(2*par1 + par2^2)
  }else if(dist_type == "weibull"){
    mean1 = par2 * gamma(1 + 1/par1)
    var1 = par2^2 * (gamma(1 + 2/par1) - gamma(1 + 1/par1)^2)
  }else if(dist_type == "exp"){
    mean1 = par2
    var1 = par2^2
  }
  return(list(mean=mean1, var=var1))
}


