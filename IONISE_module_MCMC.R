library(deSolve)
library(MASS)
library(mvtnorm)
library(dplyr)
library(invgamma)
library(ramcmc)


MCMC_function <- function(data_R, y_init, param_set, initial_phase, estim_infect,
                          prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                          nrepeat, tun_beta, dist_type, model_type, lik_type){
  if(model_type == "Extended_SEIR"){
    if(lik_type == "gaussian"){
      mcmc.result = MCMC_function_extended(data_daily = data_R, y_init, param_set, estim_infect,
                                           prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                           nrepeat, tun_beta, dist_type)
    }else{
      mcmc.result = MCMC_function_extended_pois(data_daily = data_R, y_init, param_set, estim_infect,
                                                prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                                nrepeat, tun_beta, dist_type)
    }
    
  }else{
    S_init = y_init$S_init
    E_init = y_init$E_init
    I_init = y_init$I_init
    R_init = y_init$R_init
    
    latent_mean = param_to_mean_var(par1 = param_set$k_tauL, par2 = param_set$s_tauL, dist_type = dist_type)$mean
    latent_var = param_to_mean_var(par1 = param_set$k_tauL, par2 = param_set$s_tauL, dist_type = dist_type)$var
    infect_mean = param_to_mean_var(par1 = param_set$k_tauI, par2 = param_set$s_tauI, dist_type = dist_type)$mean
    infect_var = param_to_mean_var(par1 = param_set$k_tauI, par2 = param_set$s_tauI, dist_type = dist_type)$var
    
    if(lik_type == "gaussian"){
      mcmc.result = MCMC_function_SEIR(data_R, S_init, E_init, I_init, R_init, 
                                       initial_phase, estim_infect,
                                       latent_mean, latent_var, infect_mean, infect_var,
                                       prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                       nrepeat, tun_beta, dist_type)
    }else{
      mcmc.result = MCMC_function_SEIR_pois(data_R, S_init, E_init, I_init, R_init, 
                                            initial_phase, estim_infect,
                                            latent_mean, latent_var, infect_mean, infect_var,
                                            prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                            nrepeat, tun_beta, dist_type)
    }
  }
  return(mcmc.result)
}

MCMC_function_SEIR <- function(data_R, S_init, E_init, I_init, R_init, 
                               initial_phase, estim_infect,
                               latent_mean, latent_var, infect_mean, infect_var,
                               prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                               nrepeat, tun_beta, dist_type, filename = "tmp.RData"){
  
  if(estim_infect){
    estim_TF = c(1, 0, 1)
  }else{
    estim_TF = c(1, 0, 0)
  }
  
  latent_shape = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type)$par1
  latent_scale = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type)$par2
  infect_shape = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type)$par1
  infect_scale = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type)$par2
  
  delayparam1 = c(latent_shape, latent_scale) # latent period (tauL) distribution parameter
  delayparam2 = c(infect_shape, infect_scale) # infectious period (tauI) distribution parameter
  
  # initial values of the parameters for MCMC sampling.
  beta_init = 0.2
  delayparam1_init = delayparam1
  delayparam2_init = delayparam2
  
  if(is.na(E_init)){
    E_init = I_init
  }
  
  daily_R = data_R[,2]
  R_true = cumsum(c(R_init, daily_R))
  
  maxT = length(R_true)-1
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1]
  
  if (dist_type == "lognormal"){
    lowbnds0 = c(-10, 0.1) # the lower bounds for the parameters of the two delays.  
  }else{
    lowbnds0 = c(1.1, 0.1) # the lower bounds for the parameters of the two delays.  
  }
  
  if(dist_type == "exp"){
    assume_exp = TRUE # FALSE: incorporating gamma distributions for latent and infectious periods.  
  }else{
    assume_exp = FALSE # FALSE: incorporating gamma distributions for latent and infectious periods.  
  }
  
  ## Initialize current R trajectory because we do not want to start with the true trajectory at the beginning.
  theta_estim = matrix(NA, nrow = nrepeat, ncol = 5) ## matrix for sampled parameters for all iterations
  ## 1st~5th columns: beta, shape1, scale1, shape2, scale2 
  theta_estim[1,] = c(beta_init, delayparam1_init, delayparam2_init)
  
  dd = get_divisor_dist(theta_estim[1,2:5], dist_type = dist_type)
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
  
  if (assume_exp == TRUE){
    init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
    theta_ode <- c(beta = theta_estim[1,1], s1 = theta_estim[1,3]^(-1), s2 = theta_estim[1,5]^(-1), N = sum(init))
    sol = ode(y = init, times=ext_tspan, func=SEIR11, parms = theta_ode)
    R_curr = sol[subseq0, dim(sol)[2]]
    R_cand = R_curr
  }else{
    if(initial_phase){
      sim_cand_ext = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = theta_estim[1,], y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type)
    }else{
      sim_cand_ext = mean_trajectory_SEIR_hist(timespan = ext_tspan, theta = theta_estim[1,], y_init = c(S_init, E_init, I_init, R_init))
    }
    R_curr = sim_cand_ext$Rt[subseq0]
    R_cand = R_curr
  }
  
  noiseType = 1
  noise_multi = 1
  noise_const = 0.1
  
  # prior parameters
  pri_beta = c(prior_mean_beta^2 / prior_var_beta, prior_mean_beta / prior_var_beta)
  pri_shape1 = c(1,1) * 0.001
  pri_scale1 = c(1,1) * 0.001
  pri_shape2 = c(prior_mean_infect_shape^2 / prior_var_infect_shape, prior_mean_infect_shape / prior_var_infect_shape)
  pri_scale2 = c(prior_mean_infect_scale^2 / prior_var_infect_scale, prior_mean_infect_scale / prior_var_infect_scale)
  
  S_delay1 <- diag(2)
  S_delay2 <- diag(2)
  tun_delay1 <- c(1.0, 1);
  tun_delay2 <- c(1.0, 1);
  
  # record whether an update for the parameters occur at each iteration.
  count_delay1 = rep(0, nrepeat)
  count_delay2 = rep(0, nrepeat)
  count_beta = rep(0, nrepeat)
  
  AR_list = rep(0, nrepeat)
  
  
  for(kk in 2:nrepeat){
    if(kk %% 100 == 2 || kk %% max(round(nrepeat/100), 1) == 0){
      cat(paste("Iteration: ", kk," out of ",nrepeat, ". ", round(kk/nrepeat * 100,1),"% done.", sep = ""))
      cat("\n")
    }
    theta_tmp = theta_estim[kk-1, ]
    
    ## STEP 1 --  Metropolis-Hasting algorithm for the transmission rate, beta
    
    repeat{
      beta_cand = rnorm(n = 1, mean = theta_tmp[1], sd = tun_beta)
      if(beta_cand > 0){
        break
      }
    }
    
    if (estim_TF[1] == 1){
      # the below three lines are needed to ensure numerical accuracy of the Heun's method.
      dd = get_divisor_dist(theta_tmp[2:5], dist_type = dist_type)
      ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
      subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
      
      if (assume_exp == TRUE){
        init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
        theta_ode <- c(beta = beta_cand, s1 = theta_tmp[3]^(-1), s2 = theta_tmp[5]^(-1), N = sum(init))
        sol = ode(y = init, times=ext_tspan, func=SEIR11, parms = theta_ode)
        R_cand = sol[subseq0, dim(sol)[2]]
      }else{
        if(initial_phase){
          sim_cand_ext = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = c(beta_cand, theta_tmp[2:5]), y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type)
        }else{
          sim_cand_ext = mean_trajectory_SEIR_hist(timespan = ext_tspan, theta = c(beta_cand, theta_tmp[2:5]), y_init = c(S_init, E_init, I_init, R_init))
        }
        R_cand = sim_cand_ext$Rt[subseq0]
      }
      
      if(noiseType == 1){
        sd_list_curr = noise_multi * sqrt(R_true - R_init) + noise_const
        sd_list_cand = sd_list_curr
      }else if(noiseType == 2){
        sd_list_curr = noise_multi * sqrt(diff(R_true)) + noise_const
        sd_list_curr = c(sd_list_curr[1] , sd_list_curr)
        sd_list_cand = sd_list_curr
      }else if(noiseType == 3){
        sd_list_curr = noise_multi * sqrt(R_curr - R_init) + noise_const
        sd_list_cand = noise_multi * sqrt(R_cand - R_init) + noise_const
      }else{
        sd_list_curr = noise_multi * sqrt(diff(R_curr)) + noise_const
        sd_list_curr = c(sd_list_curr[1] , sd_list_curr)
        sd_list_cand = noise_multi * sqrt(diff(R_cand)) + noise_const
        sd_list_cand = c(sd_list_cand[1] , sd_list_cand)
      }
      
      loglik_curr = sum(dnorm(R_curr, mean = R_true, sd = sd_list_curr, log = TRUE))
      loglik_cand = sum(dnorm(R_cand, mean = R_true, sd = sd_list_cand, log = TRUE))
      
      logpri_curr = dgamma(theta_tmp[1], shape = pri_beta[1], rate = pri_beta[2],log = TRUE)
      logpri_cand = dgamma(beta_cand, shape = pri_beta[1], rate = pri_beta[2],log = TRUE)
      
      logprop_curr = - pnorm(beta_cand, mean = 0, sd = tun_beta, log.p = TRUE)
      logprop_cand = - pnorm(theta_tmp[1], mean = 0, sd = tun_beta, log.p = TRUE)
      
      accept_ratio = exp(loglik_cand + logpri_cand + logprop_curr - (loglik_curr + logpri_curr + logprop_cand))
      
      u0 = runif(1)
      if (u0 < accept_ratio){
        R_curr = R_cand
        count_beta[kk] = 1
        theta_tmp[1] = beta_cand
      }
    }
    
    # cat("s-2 ")
    ## STEP 2 -- Metropolis Hastings sampling for the shape and scale parameters for the latent period distribution ====================
    
    if (assume_exp == TRUE){
      if(estim_TF[2] == 1){
        if(initial_phase){
          p.update = MH.delay1.ABC.exp(P = theta_tmp[3], S = S_delay1[2,2], rep = kk, R_true = R_true, R_curr = R_curr, sd_list = sd_list,
                                       tun = tun_delay1, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam2 = theta_tmp[4:5],
                                       pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0,
                                       noiseType = noiseType, noise_multi=noise_multi, noise_const=noise_const)
        }else{
          p.update = MH.delay1.ABC.exp.hist(P = theta_tmp[3], S = S_delay1[2,2], rep = kk, R_true = R_true, R_curr = R_curr, sd_list = sd_list,
                                            tun = tun_delay1, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam2 = theta_tmp[4:5],
                                            pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0,
                                            noiseType = noiseType, noise_multi=noise_multi, noise_const=noise_const)
        }
        theta_tmp[3] = p.update$P
        S_delay1[2,2] = p.update$S
        count_delay1[kk] = p.update$count
        R_curr = p.update$R_new
      }
    }else{
      if(estim_TF[2] == 1){
        if(initial_phase){
          p.update = MH.delay1.ABC(P = theta_tmp[2:3], S = S_delay1, rep = kk, R_true = R_true, R_curr = R_curr, sd_list = sd_list,
                                   tun = tun_delay1, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam2 = theta_tmp[4:5],
                                   pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0,
                                   noiseType = noiseType, noise_multi=noise_multi, noise_const=noise_const, dist_type = dist_type)
        }else{
          p.update = MH.delay1.ABC.hist(P = theta_tmp[2:3], S = S_delay1, rep = kk, R_true = R_true, R_curr = R_curr, sd_list = sd_list,
                                        tun = tun_delay1, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam2 = theta_tmp[4:5],
                                        pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0,
                                        noiseType = noiseType, noise_multi=noise_multi, noise_const=noise_const)
        }
        theta_tmp[2:3] = p.update$P
        S_delay1 = p.update$S
        count_delay1[kk] = p.update$count
        R_curr = p.update$R_new
      }
    }
    
    # cat("s-3 ")
    ## STEP 3 -- Metropolis Hastings sampling for the shape and scale parameters for the infectious period distribution====================
    if(assume_exp){
      if(estim_TF[3] == 1){
        if(initial_phase){
          p.update = MH.delay2.ABC.exp(P = theta_tmp[5], S = S_delay2[2,2], rep = kk, R_true = R_true, R_curr = R_curr, sd_list = sd_list,
                                       tun = tun_delay2, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam1 = theta_tmp[2:3],
                                       pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0,
                                       noiseType = noiseType, noise_multi=noise_multi, noise_const=noise_const)
        }else{
          p.update = MH.delay2.ABC.exp.hist(P = theta_tmp[5], S = S_delay2[2,2], rep = kk, R_true = R_true, R_curr = R_curr, sd_list = sd_list,
                                            tun = tun_delay2, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam1 = theta_tmp[2:3],
                                            pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0,
                                            noiseType = noiseType, noise_multi=noise_multi, noise_const=noise_const)
        }
        theta_tmp[5] = p.update$P
        S_delay2[2,2] = p.update$S
        count_delay2[kk] = p.update$count
        R_curr = p.update$R_new
      }
    }else{
      if(estim_TF[3] == 1){
        if(initial_phase){
          p.update = MH.delay2.ABC(P = theta_tmp[4:5], S = S_delay2, rep = kk, R_true = R_true, R_curr = R_curr, sd_list = sd_list,
                                   tun = tun_delay2, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam1 = theta_tmp[2:3],
                                   pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0,
                                   noiseType = noiseType, noise_multi=noise_multi, noise_const=noise_const, dist_type = dist_type)
        }else{
          p.update = MH.delay2.ABC.hist(P = theta_tmp[4:5], S = S_delay2, rep = kk, R_true = R_true, R_curr = R_curr, sd_list = sd_list,
                                        tun = tun_delay2, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam1 = theta_tmp[2:3],
                                        pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0,
                                        noiseType = noiseType, noise_multi=noise_multi, noise_const=noise_const)
        }
        theta_tmp[4:5] = p.update$P
        S_delay2 = p.update$S
        count_delay2[kk] = p.update$count
        R_curr = p.update$R_new
      }
    }
    theta_estim[kk,] = theta_tmp
    
    if(kk == 100){
      save(list = ls(all.names = TRUE), file = filename)  
      #save.image(file = filename)
    }
    
    if(kk %% 2000 == 0){
      #save(list = ls(all.names = TRUE), file = filename)  
      #save.image(file = filename)
    }      
    
  }
  return(theta_estim)
}



MCMC_function_SEIR_pois <- function(data_R, S_init, E_init, I_init, R_init, 
                                    initial_phase, estim_infect,
                                    latent_mean, latent_var, infect_mean, infect_var, 
                                    prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                    nrepeat, tun_beta, dist_type, filename = "tmp.RData"){
  
  if(estim_infect){
    estim_TF = c(1, 0, 1)
  }else{
    estim_TF = c(1, 0, 0)
  }
  
  latent_shape = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type)$par1
  latent_scale = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type)$par2
  infect_shape = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type)$par1
  infect_scale = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type)$par2
  
  delayparam1 = c(latent_shape, latent_scale) # latent period (tauL) distribution parameter
  delayparam2 = c(infect_shape, infect_scale) # infectious period (tauI) distribution parameter
  
  # initial values of the parameters for MCMC sampling.
  beta_init = 0.2
  delayparam1_init = delayparam1
  delayparam2_init = delayparam2
  
  if(is.na(E_init)){
    E_init = I_init
  }
  
  daily_R = data_R[,2]
  R_true = round(cumsum(c(R_init, daily_R)))
  #R_true = round(daily_R)
  
  maxT = length(R_true)-1
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1]
  
  if (dist_type == "lognormal"){
    lowbnds0 = c(-10, 0.1) # the lower bounds for the parameters of the two delays.  
  }else{
    lowbnds0 = c(1.1, 0.1) # the lower bounds for the parameters of the two delays.  
  }
  
  if(dist_type == "exp"){
    assume_exp = TRUE # FALSE: incorporating gamma distributions for latent and infectious periods.  
  }else{
    assume_exp = FALSE # FALSE: incorporating gamma distributions for latent and infectious periods.  
  }
  
  
  ## Initialize current R trajectory because we do not want to start with the true trajectory at the beginning.
  theta_estim = matrix(NA, nrow = nrepeat, ncol = 5) ## matrix for sampled parameters for all iterations
  ## 1st~5th columns: beta, shape1, scale1, shape2, scale2 
  theta_estim[1,] = c(beta_init, delayparam1_init, delayparam2_init)
  
  dd = get_divisor_dist(theta_estim[1,2:5], dist_type = dist_type)
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
  
  if (assume_exp == TRUE){
    init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
    theta_ode <- c(beta = theta_estim[1,1], s1 = theta_estim[1,3]^(-1), s2 = theta_estim[1,5]^(-1), N = sum(init))
    sol = ode(y = init, times=ext_tspan, func=SEIR11, parms = theta_ode)
    R_curr = sol[subseq0, dim(sol)[2]]
    R_cand = R_curr
  }else{
    if(initial_phase){
      sim_cand_ext = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = theta_estim[1,], y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type)
    }else{
      sim_cand_ext = mean_trajectory_SEIR_hist(timespan = ext_tspan, theta = theta_estim[1,], y_init = c(S_init, E_init, I_init, R_init))
    }
    R_curr = sim_cand_ext$Rt[subseq0]
    R_cand = R_curr
  }
  
  # prior parameters
  pri_beta = c(prior_mean_beta^2 / prior_var_beta, prior_mean_beta / prior_var_beta)
  pri_shape1 = c(1,1) * 0.001
  pri_scale1 = c(1,1) * 0.001
  pri_shape2 = c(prior_mean_infect_shape^2 / prior_var_infect_shape, prior_mean_infect_shape / prior_var_infect_shape)
  pri_scale2 = c(prior_mean_infect_scale^2 / prior_var_infect_scale, prior_mean_infect_scale / prior_var_infect_scale)
  
  S_delay1 <- diag(2)
  S_delay2 <- diag(2)
  tun_delay1 <- c(1.0, 1);
  tun_delay2 <- c(1.0, 1);
  
  # record whether an update for the parameters occur at each iteration.
  count_delay1 = rep(0, nrepeat)
  count_delay2 = rep(0, nrepeat)
  count_beta = rep(0, nrepeat)
  
  AR_list = rep(0, nrepeat)
  
  
  for(kk in 2:nrepeat){
    if(kk %% 100 == 2 || kk %% max(round(nrepeat/100), 1) == 0){
      cat(paste("Iteration: ", kk," out of ",nrepeat, ". ", round(kk/nrepeat * 100,1),"% done.", sep = ""))
      cat("\n")
    }
    theta_tmp = theta_estim[kk-1, ]
    
    ## STEP 1 --  Metropolis-Hasting algorithm for the transmission rate, beta
    
    repeat{
      beta_cand = rnorm(n = 1, mean = theta_tmp[1], sd = tun_beta)
      if(beta_cand > 0){
        break
      }
    }
    
    if (estim_TF[1] == 1){
      # the below three lines are needed to ensure numerical accuracy of the Heun's method.
      #dd = get_divisor_Gamma(theta_tmp[2:5])
      dd = get_divisor_dist(theta_tmp[2:5], dist_type = dist_type)
      ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
      subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
      
      if (assume_exp == TRUE){
        init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
        theta_ode <- c(beta = beta_cand, s1 = theta_tmp[3]^(-1), s2 = theta_tmp[5]^(-1), N = sum(init))
        sol = ode(y = init, times=ext_tspan, func=SEIR11, parms = theta_ode)
        R_cand = sol[subseq0, dim(sol)[2]]
      }else{
        if(initial_phase){
          sim_cand_ext = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = c(beta_cand, theta_tmp[2:5]), y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type)
        }else{
          sim_cand_ext = mean_trajectory_SEIR_hist(timespan = ext_tspan, theta = c(beta_cand, theta_tmp[2:5]), y_init = c(S_init, E_init, I_init, R_init))
        }
        R_cand = sim_cand_ext$Rt[subseq0]
      }
      
      loglik_curr = sum(dpois(R_true, lambda = R_curr, log = TRUE))
      loglik_cand = sum(dpois(R_true, lambda = R_cand, log = TRUE))
      
      logpri_curr = dgamma(theta_tmp[1], shape = pri_beta[1], rate = pri_beta[2], log = TRUE)
      logpri_cand = dgamma(beta_cand, shape = pri_beta[1], rate = pri_beta[2], log = TRUE)
      
      logprop_curr = - pnorm(beta_cand, mean = 0, sd = tun_beta, log.p = TRUE)
      logprop_cand = - pnorm(theta_tmp[1], mean = 0, sd = tun_beta, log.p = TRUE)
      
      accept_ratio = exp(loglik_cand + logpri_cand + logprop_curr - (loglik_curr + logpri_curr + logprop_cand))
      
      u0 = runif(1)
      if (u0 < accept_ratio){
        R_curr = R_cand
        count_beta[kk] = 1
        theta_tmp[1] = beta_cand
      }
    }
    
    # cat("s-2 ")
    ## STEP 2 -- Metropolis Hastings sampling for the shape and scale parameters for the latent period distribution ====================
    
    if (assume_exp == TRUE){
      if(estim_TF[2] == 1){
        if(initial_phase){
          p.update = MH.delay1.ABC.exp.pois(P = theta_tmp[3], S = S_delay1[2,2], rep = kk, R_true = R_true, R_curr = R_curr, 
                                            tun = tun_delay1, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam2 = theta_tmp[4:5],
                                            pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0)
        }else{
          p.update = MH.delay1.ABC.exp.pois.hist(P = theta_tmp[3], S = S_delay1[2,2], rep = kk, R_true = R_true, R_curr = R_curr, 
                                                 tun = tun_delay1, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam2 = theta_tmp[4:5],
                                                 pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0)
        }
        theta_tmp[3] = p.update$P
        S_delay1[2,2] = p.update$S
        count_delay1[kk] = p.update$count
        R_curr = p.update$R_new
      }
    }else{
      if(estim_TF[2] == 1){
        if(initial_phase){
          p.update = MH.delay1.ABC.pois(P = theta_tmp[2:3], S = S_delay1, rep = kk, R_true = R_true, R_curr = R_curr, 
                                        tun = tun_delay1, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam2 = theta_tmp[4:5],
                                        pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0, dist_type = dist_type)
        }else{
          p.update = MH.delay1.ABC.pois.hist(P = theta_tmp[2:3], S = S_delay1, rep = kk, R_true = R_true, R_curr = R_curr,
                                             tun = tun_delay1, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam2 = theta_tmp[4:5],
                                             pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0)
        }
        theta_tmp[2:3] = p.update$P
        S_delay1 = p.update$S
        count_delay1[kk] = p.update$count
        R_curr = p.update$R_new
      }
    }
    
    # cat("s-3 ")
    ## STEP 3 -- Metropolis Hastings sampling for the shape and scale parameters for the infectious period distribution====================
    if(assume_exp){
      if(estim_TF[3] == 1){
        if(initial_phase){
          p.update = MH.delay2.ABC.exp.pois(P = theta_tmp[5], S = S_delay2[2,2], rep = kk, R_true = R_true, R_curr = R_curr, 
                                            tun = tun_delay2, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam1 = theta_tmp[2:3],
                                            pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0)
        }else{
          p.update = MH.delay2.ABC.exp.pois.hist(P = theta_tmp[5], S = S_delay2[2,2], rep = kk, R_true = R_true, R_curr = R_curr, 
                                                 tun = tun_delay2, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam1 = theta_tmp[2:3],
                                                 pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0)
        }
        theta_tmp[5] = p.update$P
        S_delay2[2,2] = p.update$S
        count_delay2[kk] = p.update$count
        R_curr = p.update$R_new
      }
    }else{
      if(estim_TF[3] == 1){
        if(initial_phase){
          p.update = MH.delay2.ABC.pois(P = theta_tmp[4:5], S = S_delay2, rep = kk, R_true = R_true, R_curr = R_curr,
                                        tun = tun_delay2, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam1 = theta_tmp[2:3],
                                        pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0, dist_type = dist_type)
        }else{
          p.update = MH.delay2.ABC.pois.hist(P = theta_tmp[4:5], S = S_delay2, rep = kk, R_true = R_true, R_curr = R_curr, 
                                             tun = tun_delay2, S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta_param = theta_tmp[1], delayparam1 = theta_tmp[2:3],
                                             pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0)
        }
        theta_tmp[4:5] = p.update$P
        S_delay2 = p.update$S
        count_delay2[kk] = p.update$count
        R_curr = p.update$R_new
      }
    }
    theta_estim[kk,] = theta_tmp
    
    if(kk == 20){
      save(list = ls(all.names = TRUE), file = filename)  
      #save.image(file = filename)
    }
    
    if(kk %% 2000 == 0){
      # save(list = ls(all.names = TRUE), file = filename)  
      # save.image(file = filename)
    }    
    
  }
  return(theta_estim)
}

MCMC_function_extended <- function(data_daily, y_init, param_set, estim_infect,
                                   prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                   nrepeat, tun_beta, dist_type, filename = "tmp.RData"){
  
  if(estim_infect){
    estim_TF = c(1, 0, 1)
  }else{
    estim_TF = c(1, 0, 0)
  }
  
  latent_shape = param_set$k_tauL
  latent_scale = param_set$s_tauL
  infect_shape = param_set$k_tauI
  infect_scale = param_set$s_tauI
  
  dd_lower = get_divisor_dist(c(param_set$k_tauLh, param_set$s_tauLh, param_set$k_tauIA, param_set$s_tauIA), dist_type = dist_type)
  
  delayparam1 = c(latent_shape, latent_scale) # latent period (tauL) distribution parameter
  delayparam2 = c(infect_shape, infect_scale) # infectious period (tauI) distribution parameter
  
  # initial values of the parameters for MCMC sampling.
  beta_init = param_set$beta
  delayparam1_init = delayparam1
  delayparam2_init = delayparam2
  
  daily_C = data_daily[,2]
  cumul_init = y_init$C_init
  cumul_true = cumsum(c(cumul_init, daily_C))
  
  maxT = length(cumul_true)-1
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1]
  
  if (dist_type == "lognormal"){
    lowbnds0 = c(-10, 0.1) # the lower bounds for the parameters of the two delays.  
  }else{
    lowbnds0 = c(1.1, 0.1) # the lower bounds for the parameters of the two delays.  
  }
  
  ## Initialize current R trajectory because we do not want to start with the true trajectory at the beginning.
  theta_estim = matrix(NA, nrow = nrepeat, ncol = 5) ## matrix for sampled parameters for all iterations
  ## 1st~5th columns: beta, shape1, scale1, shape2, scale2 
  theta_estim[1,] = c(beta_init, delayparam1_init, delayparam2_init)
  
  dd = max(dd_lower, get_divisor_dist(theta_estim[1,2:5], dist_type = dist_type))
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
  
  sim_cand_ext = mean_trajectory_LARGE_dist_input(timespan = ext_tspan, param_set = param_set, y_init = y_init, dist_type = dist_type)
  cumul_curr = sim_cand_ext$cumul_Ct[subseq0]
  cumul_cand = cumul_curr
  
  noise_multi = 1
  noise_const = 0.1
  
  # prior parameters
  pri_beta = c(prior_mean_beta^2 / prior_var_beta, prior_mean_beta / prior_var_beta)
  pri_shape1 = c(1,1) * 0.001
  pri_scale1 = c(1,1) * 0.001
  pri_shape2 = c(prior_mean_infect_shape^2 / prior_var_infect_shape, prior_mean_infect_shape / prior_var_infect_shape)
  pri_scale2 = c(prior_mean_infect_scale^2 / prior_var_infect_scale, prior_mean_infect_scale / prior_var_infect_scale)
  
  S_delay1 <- diag(2)
  S_delay2 <- diag(2)
  tun_delay1 <- c(1.0, 1);
  tun_delay2 <- c(1.0, 1);
  
  # record whether an update for the parameters occur at each iteration.
  count_delay1 = rep(0, nrepeat)
  count_delay2 = rep(0, nrepeat)
  count_beta = rep(0, nrepeat)
  
  AR_list = rep(0, nrepeat)
  
  param_set_tmp = param_set
  param_set_cand = param_set_tmp
  
  for(kk in 2:nrepeat){
    if(kk %% 100 == 2 || kk %% max(round(nrepeat/100), 1) == 0){
      cat(paste("Iteration: ", kk," out of ",nrepeat, ". ", round(kk/nrepeat * 100,1),"% done.", sep = ""))
      cat("\n")
    }
    theta_tmp = theta_estim[kk-1, ]
    param_set_tmp$beta = theta_tmp[1]
    param_set_tmp$k_tauL = theta_tmp[2]
    param_set_tmp$s_tauL = theta_tmp[3]
    param_set_tmp$k_tauI = theta_tmp[4]
    param_set_tmp$s_tauI = theta_tmp[5]
    param_set_cand = param_set_tmp
    
    ## STEP 1 --  Metropolis-Hasting algorithm for the transmission rate, beta
    
    repeat{
      beta_cand = rnorm(n = 1, mean = theta_tmp[1], sd = tun_beta)
      if(beta_cand > 0){
        break
      }
    }
    param_set_cand$beta = beta_cand
    
    if (estim_TF[1] == 1){
      # the below three lines are needed to ensure numerical accuracy of the Heun's method.
      dd = max(dd_lower, get_divisor_dist(c(param_set_cand$k_tauL, param_set_cand$s_tauL, param_set_cand$k_tauI, param_set_cand$s_tauI), dist_type = dist_type))
      ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
      subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
      
      sim_cand_ext = mean_trajectory_LARGE_dist_input(timespan = ext_tspan, param_set = param_set_cand, y_init = y_init, dist_type = dist_type)
      cumul_cand = sim_cand_ext$cumul_Ct[subseq0]
      
      sd_list_curr = noise_multi * sqrt(cumul_true - cumul_init) + noise_const
      sd_list_cand = sd_list_curr
      
      loglik_curr = sum(dnorm(cumul_curr, mean = cumul_true, sd = sd_list_curr, log = TRUE))
      loglik_cand = sum(dnorm(cumul_cand, mean = cumul_true, sd = sd_list_cand, log = TRUE))
      
      logpri_curr = dgamma(param_set_tmp$beta, shape = pri_beta[1], rate = pri_beta[2],log = TRUE)
      logpri_cand = dgamma(param_set_cand$beta, shape = pri_beta[1], rate = pri_beta[2],log = TRUE)
      
      logprop_curr = - pnorm(param_set_cand$beta, mean = 0, sd = tun_beta, log.p = TRUE)
      logprop_cand = - pnorm(param_set_tmp$beta, mean = 0, sd = tun_beta, log.p = TRUE)
      
      accept_ratio = exp(loglik_cand + logpri_cand + logprop_curr - (loglik_curr + logpri_curr + logprop_cand))
      
      u0 = runif(1)
      if (u0 < accept_ratio){
        cumul_curr = cumul_cand
        count_beta[kk] = 1
        theta_tmp[1] = param_set_cand$beta
        param_set_tmp = param_set_cand
      }
    }
    
    # cat("s-2 ")
    ## STEP 2 -- Metropolis Hastings sampling for the shape and scale parameters for the latent period distribution ====================
    if(estim_TF[2] == 1){
      p.update = MH.delay1.ABC.extended(param_set = param_set_tmp, S = S_delay1, rep = kk, cumul_true = cumul_true, cumul_curr = cumul_curr, tun = tun_delay1, y_init = y_init,
                                        pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0, dd_lower,
                                        noise_multi = noise_multi, noise_const = noise_const, dist_type = dist_type)
      
      param_set_tmp = p.update$param_set
      theta_tmp[2] = param_set_tmp$k_tauL
      theta_tmp[3] = param_set_tmp$s_tauL
      S_delay1 = p.update$S
      count_delay1[kk] = p.update$count
      cumul_curr = p.update$cumul_new
    }
    
    # cat("s-3 ")
    ## STEP 3 -- Metropolis Hastings sampling for the shape and scale parameters for the infectious period distribution====================
    if(estim_TF[3] == 1){
      p.update = MH.delay2.ABC.extended(param_set = param_set_tmp, S = S_delay2, rep = kk, cumul_true = cumul_true, cumul_curr = cumul_curr, tun = tun_delay2, y_init = y_init,
                                        pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0, dd_lower,
                                        noise_multi = noise_multi, noise_const = noise_const, dist_type = dist_type)
      
      param_set_tmp = p.update$param_set
      theta_tmp[4] = param_set_tmp$k_tauI
      theta_tmp[5] = param_set_tmp$s_tauI
      S_delay2 = p.update$S
      count_delay2[kk] = p.update$count
      cumul_curr = p.update$cumul_new
    }
    
    theta_estim[kk,] = theta_tmp
    
    if(kk == 5){
      save(list = ls(all.names = TRUE), file = filename)  
      #save.image(file = filename)
    }
    
    if(kk %% 2000 == 0 || kk == nrepeat){
      # save(list = ls(all.names = TRUE), file = filename)  
      #save.image(file = filename)
    }      
    
  }
  return(theta_estim)
}


MCMC_function_extended_pois <- function(data_daily, y_init, param_set, estim_infect,
                                        prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                        nrepeat, tun_beta, dist_type, filename = "tmp.RData"){
  
  if(estim_infect){
    estim_TF = c(1, 0, 1)
  }else{
    estim_TF = c(1, 0, 0)
  }
  
  latent_shape = param_set$k_tauL
  latent_scale = param_set$s_tauL
  infect_shape = param_set$k_tauI
  infect_scale = param_set$s_tauI
  
  dd_lower = get_divisor_dist(c(param_set$k_tauLh, param_set$s_tauLh, param_set$k_tauIA, param_set$s_tauIA), dist_type = dist_type)
  
  delayparam1 = c(latent_shape, latent_scale) # latent period (tauL) distribution parameter
  delayparam2 = c(infect_shape, infect_scale) # infectious period (tauI) distribution parameter
  
  # initial values of the parameters for MCMC sampling.
  beta_init = param_set$beta
  delayparam1_init = delayparam1
  delayparam2_init = delayparam2
  
  daily_C = data_daily[,2]
  cumul_init = y_init$C_init
  cumul_true = round(cumsum(c(cumul_init, daily_C)))
  
  maxT = length(cumul_true)-1
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1]
  
  if (dist_type == "lognormal"){
    lowbnds0 = c(-10, 0.1) # the lower bounds for the parameters of the two delays.  
  }else{
    lowbnds0 = c(1.1, 0.1) # the lower bounds for the parameters of the two delays.  
  }
  
  ## Initialize current R trajectory because we do not want to start with the true trajectory at the beginning.
  theta_estim = matrix(NA, nrow = nrepeat, ncol = 5) ## matrix for sampled parameters for all iterations
  ## 1st~5th columns: beta, shape1, scale1, shape2, scale2 
  theta_estim[1,] = c(beta_init, delayparam1_init, delayparam2_init)
  
  dd = max(dd_lower, get_divisor_dist(theta_estim[1,2:5], dist_type = dist_type))
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
  
  sim_cand_ext = mean_trajectory_LARGE_dist_input(timespan = ext_tspan, param_set = param_set, y_init = y_init, dist_type = dist_type)
  cumul_curr = sim_cand_ext$cumul_Ct[subseq0]
  cumul_cand = cumul_curr
  
  # prior parameters
  pri_beta = c(prior_mean_beta^2 / prior_var_beta, prior_mean_beta / prior_var_beta)
  pri_shape1 = c(1,1) * 0.001
  pri_scale1 = c(1,1) * 0.001
  pri_shape2 = c(prior_mean_infect_shape^2 / prior_var_infect_shape, prior_mean_infect_shape / prior_var_infect_shape)
  pri_scale2 = c(prior_mean_infect_scale^2 / prior_var_infect_scale, prior_mean_infect_scale / prior_var_infect_scale)
  
  S_delay1 <- diag(2)
  S_delay2 <- diag(2)
  tun_delay1 <- c(1.0, 1);
  tun_delay2 <- c(1.0, 1);
  
  # record whether an update for the parameters occur at each iteration.
  count_delay1 = rep(0, nrepeat)
  count_delay2 = rep(0, nrepeat)
  count_beta = rep(0, nrepeat)
  
  AR_list = rep(0, nrepeat)
  
  param_set_tmp = param_set
  param_set_cand = param_set_tmp
  
  for(kk in 2:nrepeat){
    if(kk %% 100 == 2 || kk %% max(round(nrepeat/100), 1) == 0){
      cat(paste("Iteration: ", kk," out of ",nrepeat, ". ", round(kk/nrepeat * 100,1),"% done.", sep = ""))
      cat("\n")
    }
    theta_tmp = theta_estim[kk-1, ]
    param_set_tmp$beta = theta_tmp[1]
    param_set_tmp$k_tauL = theta_tmp[2]
    param_set_tmp$s_tauL = theta_tmp[3]
    param_set_tmp$k_tauI = theta_tmp[4]
    param_set_tmp$s_tauI = theta_tmp[5]
    param_set_cand = param_set_tmp
    
    ## STEP 1 --  Metropolis-Hasting algorithm for the transmission rate, beta
    repeat{
      beta_cand = rnorm(n = 1, mean = theta_tmp[1], sd = tun_beta)
      if(beta_cand > 0){
        break
      }
    }
    param_set_cand$beta = beta_cand
    if (estim_TF[1] == 1){
      # the below three lines are needed to ensure numerical accuracy of the Heun's method.
      dd = max(dd_lower, get_divisor_dist(c(param_set_cand$k_tauL, param_set_cand$s_tauL, param_set_cand$k_tauI, param_set_cand$s_tauI), dist_type = dist_type))
      ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
      subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
      
      sim_cand_ext = mean_trajectory_LARGE_dist_input(timespan = ext_tspan, param_set = param_set_cand, y_init = y_init, dist_type = dist_type)
      cumul_cand = sim_cand_ext$cumul_Ct[subseq0]
      
      loglik_curr = sum(dpois(cumul_true, lambda = cumul_curr, log = TRUE))
      loglik_cand = sum(dpois(cumul_true, lambda = cumul_cand, log = TRUE))
      
      logpri_curr = dgamma(param_set_tmp$beta, shape = pri_beta[1], rate = pri_beta[2],log = TRUE)
      logpri_cand = dgamma(param_set_cand$beta, shape = pri_beta[1], rate = pri_beta[2],log = TRUE)
      
      logprop_curr = - pnorm(param_set_cand$beta, mean = 0, sd = tun_beta, log.p = TRUE)
      logprop_cand = - pnorm(param_set_tmp$beta, mean = 0, sd = tun_beta, log.p = TRUE)
      
      accept_ratio = exp(loglik_cand + logpri_cand + logprop_curr - (loglik_curr + logpri_curr + logprop_cand))
      
      u0 = runif(1)
      if (u0 < accept_ratio){
        cumul_curr = cumul_cand
        count_beta[kk] = 1
        theta_tmp[1] = param_set_cand$beta
        param_set_tmp = param_set_cand
      }
    }
    
    # cat("s-2 ")
    ## STEP 2 -- Metropolis Hastings sampling for the shape and scale parameters for the latent period distribution ====================
    if(estim_TF[2] == 1){
      p.update = MH.delay1.ABC.extended.pois(param_set = param_set_tmp, S = S_delay1, rep = kk, cumul_true = cumul_true, cumul_curr = cumul_curr, tun = tun_delay1, y_init = y_init,
                                             pri.shape = pri_shape1, pri.scale = pri_scale1, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0, dd_lower, dist_type = dist_type)
      
      param_set_tmp = p.update$param_set
      theta_tmp[2] = param_set_tmp$k_tauL
      theta_tmp[3] = param_set_tmp$s_tauL
      S_delay1 = p.update$S
      count_delay1[kk] = p.update$count
      cumul_curr = p.update$cumul_new
    }
    
    # cat("s-3 ")
    ## STEP 3 -- Metropolis Hastings sampling for the shape and scale parameters for the infectious period distribution====================
    if(estim_TF[3] == 1){
      p.update = MH.delay2.ABC.extended.pois(param_set = param_set_tmp, S = S_delay2, rep = kk, cumul_true = cumul_true, cumul_curr = cumul_curr, tun = tun_delay2, y_init = y_init,
                                             pri.shape = pri_shape2, pri.scale = pri_scale2, maxT = maxT, flatpri = FALSE, lowbnds = lowbnds0, dd_lower, dist_type = dist_type)
      
      param_set_tmp = p.update$param_set
      theta_tmp[4] = param_set_tmp$k_tauI
      theta_tmp[5] = param_set_tmp$s_tauI
      S_delay2 = p.update$S
      count_delay2[kk] = p.update$count
      cumul_curr = p.update$cumul_new
    }
    
    theta_estim[kk,] = theta_tmp
    
    if(kk == 5){
      save(list = ls(all.names = TRUE), file = filename)  
      #save.image(file = filename)
    }
    
    if(kk %% 2000 == 0 || kk == nrepeat){
      # save(list = ls(all.names = TRUE), file = filename)  
      #save.image(file = filename)
    }      
    
  }
  return(theta_estim)
}


