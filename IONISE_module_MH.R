library(deSolve)
library(MASS)
library(mvtnorm)
library(dplyr)
library(invgamma)
library(ramcmc)

MH.delay1.ABC <- function(P,S,rep, R_true, R_curr, sd_list, tun, S_init, E_init, I_init, R_init, 
                          beta_param, delayparam2, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), 
                          noiseType, noise_multi, noise_const, dist_type) {
  count = 0
  
  
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1] 
  dd = get_divisor_dist(c(delayparam2, P.star), dist_type = dist_type)
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = c(beta_param, P.star, delayparam2), y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type)
  R_cand = sim_cand_ext$Rt[subseq0]
  
  
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
  
  l_lik_curr = sum(dnorm(R_curr, mean = R_true, sd = sd_list_curr, log = TRUE))
  l_lik_cand = sum(dnorm(R_cand, mean = R_true, sd = sd_list_cand, log = TRUE))
  
  l.prior1.st <- dgamma(P.star[1], shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
  
}



MH.delay1.ABC.pois <- function(P,S,rep, R_true, R_curr, tun, S_init, E_init, I_init, R_init, 
                               beta_param, delayparam2, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), dist_type) {
  count = 0
  
  
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1] 
  dd = get_divisor_dist(c(delayparam2, P.star), dist_type = dist_type)
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = c(beta_param, P.star, delayparam2), y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type)
  R_cand = sim_cand_ext$Rt[subseq0]
  
  l_lik_curr = sum(dpois(R_true, lambda = R_curr, log = TRUE))
  l_lik_cand = sum(dpois(R_true, lambda = R_cand, log = TRUE))
  
  l.prior1.st <- dgamma(P.star[1], shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
  
}


MH.delay2.ABC <- function(P,S,rep, R_true, R_curr, sd_list, tun, S_init, E_init, I_init, R_init, 
                          beta_param, delayparam1, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), 
                          noiseType, noise_multi, noise_const, dist_type = dist_type) {
  count = 0
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1]
  dd = get_divisor_dist(c(delayparam1, P.star), dist_type = dist_type)
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = c(beta_param, delayparam1, P.star), y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type)
  R_cand = sim_cand_ext$Rt[subseq0]
  
  
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
  
  l_lik_curr = sum(dnorm(R_curr, mean = R_true, sd = sd_list_curr, log = TRUE))
  l_lik_cand = sum(dnorm(R_cand, mean = R_true, sd = sd_list_cand, log = TRUE))
  
  
  l.prior1.st <- dgamma(P.star[1], shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}

MH.delay2.ABC.pois <- function(P,S,rep, R_true, R_curr, tun, S_init, E_init, I_init, R_init, 
                               beta_param, delayparam1, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), dist_type) {
  count = 0
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1]
  dd = get_divisor_dist(c(delayparam1, P.star), dist_type = dist_type)
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = c(beta_param, delayparam1, P.star), y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type)
  R_cand = sim_cand_ext$Rt[subseq0]
  
  l_lik_curr = sum(dpois(R_true, lambda = R_curr, log = TRUE))
  l_lik_cand = sum(dpois(R_true, lambda = R_cand, log = TRUE))
  
  l.prior1.st <- dgamma(P.star[1], shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}






MH.delay1.ABC.exp <- function(P,S,rep, R_true, R_curr, sd_list, tun, S_init, E_init, I_init, R_init, 
                              beta_param, delayparam2, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), 
                              noiseType, noise_multi, noise_const){
  count = 0
  repeat{
    u = rnorm(1, 0 ,tun[2])
    P.star = P + S*u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  
  init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
  theta_ode <- c(beta = beta_param, s1 = P.star^(-1), s2 = delayparam2[2]^(-1), N = sum(init))
  sol = ode(y = init, times=timespan0, func=SEIR11, parms = theta_ode)
  R_cand = sol[, dim(sol)[2]]
  
  
  
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
  
  l_lik_curr = sum(dnorm(R_curr, mean = R_true, sd = sd_list_curr, log = TRUE))
  l_lik_cand = sum(dnorm(R_cand, mean = R_true, sd = sd_list_cand, log = TRUE))
  l.prior2.st <- dgamma(P.star, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior2.st - l.prior2
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}



MH.delay1.ABC.exp.pois <- function(P,S,rep, R_true, R_curr, tun, S_init, E_init, I_init, R_init, 
                                   beta_param, delayparam2, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0)){
  count = 0
  repeat{
    u = rnorm(1, 0 ,tun[2])
    P.star = P + S*u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  
  init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
  theta_ode <- c(beta = beta_param, s1 = P.star^(-1), s2 = delayparam2[2]^(-1), N = sum(init))
  sol = ode(y = init, times=timespan0, func=SEIR11, parms = theta_ode)
  R_cand = sol[, dim(sol)[2]]
  
  l_lik_curr = sum(dpois(R_true, lambda = R_curr, log = TRUE))
  l_lik_cand = sum(dpois(R_true, lambda = R_cand, log = TRUE))
  
  l.prior2.st <- dgamma(P.star, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior2.st - l.prior2
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}

MH.delay2.ABC.exp <- function(P,S,rep, R_true, R_curr, sd_list, tun, S_init, E_init, I_init, R_init,
                              beta_param, delayparam1, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), 
                              noiseType, noise_multi, noise_const) {
  count = 0
  repeat{
    u = rnorm(1, 0 ,tun[2])
    P.star = P + S*u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  
  init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
  theta_ode <- c(beta = beta_param, s1 = delayparam1[2]^(-1), s2 = P.star^(-1), N = sum(init))
  sol = ode(y = init, times=timespan0, func=SEIR11, parms = theta_ode)
  R_cand = sol[, dim(sol)[2]]
  
  
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
  
  l_lik_curr = sum(dnorm(R_curr, mean = R_true, sd = sd_list_curr, log = TRUE))
  l_lik_cand = sum(dnorm(R_cand, mean = R_true, sd = sd_list_cand, log = TRUE))
  
  l.prior2.st <- dgamma(P.star, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior2.st - l.prior2
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}

MH.delay2.ABC.exp.pois <- function(P,S,rep, R_true, R_curr, tun, S_init, E_init, I_init, R_init,
                                   beta_param, delayparam1, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0)) {
  count = 0
  repeat{
    u = rnorm(1, 0 ,tun[2])
    P.star = P + S*u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  
  init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
  theta_ode <- c(beta = beta_param, s1 = delayparam1[2]^(-1), s2 = P.star^(-1), N = sum(init))
  sol = ode(y = init, times=timespan0, func=SEIR11, parms = theta_ode)
  R_cand = sol[, dim(sol)[2]]
  
  l_lik_curr = sum(dpois(R_true, lambda = R_curr, log = TRUE))
  l_lik_cand = sum(dpois(R_true, lambda = R_cand, log = TRUE))
  
  l.prior2.st <- dgamma(P.star, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior2.st - l.prior2
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}


MH.delay1.ABC.hist <- function(P,S,rep, R_true, R_curr, sd_list, tun, S_init, E_init, I_init, R_init, 
                               beta_param, delayparam2, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), 
                               noiseType, noise_multi, noise_const) {
  count = 0
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1] 
  dd = get_divisor_dist(c(delayparam2, P.star), dist_type = "gamma")
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_SEIR_hist(timespan = ext_tspan, theta = c(beta_param, P.star, delayparam2), y_init = c(S_init, E_init, I_init, R_init))
  R_cand = sim_cand_ext$Rt[subseq0]
  
  
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
  
  l_lik_curr = sum(dnorm(R_curr, mean = R_true, sd = sd_list_curr, log = TRUE))
  l_lik_cand = sum(dnorm(R_cand, mean = R_true, sd = sd_list_cand, log = TRUE))
  
  l.prior1.st <- dgamma(P.star[1], shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
  
}


MH.delay1.ABC.pois.hist <- function(P,S,rep, R_true, R_curr, tun, S_init, E_init, I_init, R_init, 
                                    beta_param, delayparam2, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0)) {
  count = 0
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1] 
  dd = get_divisor_dist(c(delayparam2, P.star), dist_type = "gamma")
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_SEIR_hist(timespan = ext_tspan, theta = c(beta_param, P.star, delayparam2), y_init = c(S_init, E_init, I_init, R_init))
  R_cand = sim_cand_ext$Rt[subseq0]
  
  
  l_lik_curr = sum(dpois(R_true, lambda = R_curr, log = TRUE))
  l_lik_cand = sum(dpois(R_true, lambda = R_cand, log = TRUE))
  
  l.prior1.st <- dgamma(P.star[1], shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
  
}


MH.delay2.ABC.hist <- function(P,S,rep, R_true, R_curr, sd_list, tun, S_init, E_init, I_init, R_init, 
                               beta_param, delayparam1, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), 
                               noiseType, noise_multi, noise_const) {
  count = 0
  
  lowbnds = c(1.1, 0.1)
  
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1]
  dd = get_divisor_dist(c(delayparam1, P.star), dist_type = "gamma")
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_SEIR_hist(timespan = ext_tspan, theta = c(beta_param, delayparam1, P.star), y_init = c(S_init, E_init, I_init, R_init))
  R_cand = sim_cand_ext$Rt[subseq0]
  
  
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
  
  l_lik_curr = sum(dnorm(R_curr, mean = R_true, sd = sd_list_curr, log = TRUE))
  l_lik_cand = sum(dnorm(R_cand, mean = R_true, sd = sd_list_cand, log = TRUE))
  
  l.prior1.st <- dgamma(P.star[1], shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}



MH.delay2.ABC.pois.hist <- function(P,S,rep, R_true, R_curr, tun, S_init, E_init, I_init, R_init, 
                                    beta_param, delayparam1, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0)){
  count = 0
  
  lowbnds = c(1.1, 0.1)
  
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1]
  dd = get_divisor_dist(c(delayparam1, P.star), dist_type = "gamma")
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_SEIR_hist(timespan = ext_tspan, theta = c(beta_param, delayparam1, P.star), y_init = c(S_init, E_init, I_init, R_init))
  R_cand = sim_cand_ext$Rt[subseq0]
  
  
  l_lik_curr = sum(dpois(R_true, lambda = R_curr, log = TRUE))
  l_lik_cand = sum(dpois(R_true, lambda = R_cand, log = TRUE))
  
  l.prior1.st <- dgamma(P.star[1], shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}


MH.delay1.ABC.exp.hist <- function(P,S,rep, R_true, R_curr, sd_list, tun, S_init, E_init, I_init, R_init, 
                                   beta_param, delayparam2, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), 
                                   noiseType, noise_multi, noise_const) {
  count = 0
  repeat{
    u = rnorm(1, 0 ,tun[2])
    P.star = P + S*u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star>lowbnds[2]){
      break
    }
  }
  
  
  timespan0 = 0:maxT
  
  init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
  theta_ode <- c(beta = beta_param, s1 = P.star^(-1), s2 = delayparam2[2]^(-1), N = sum(init))
  sol = ode(y = init, times=timespan0, func=SEIR11, parms = theta_ode)
  R_cand = sol[, dim(sol)[2]]
  
  
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
  
  l_lik_curr = sum(dnorm(R_curr, mean = R_true, sd = sd_list_curr, log = TRUE))
  l_lik_cand = sum(dnorm(R_cand, mean = R_true, sd = sd_list_cand, log = TRUE))
  
  l.prior2.st <- dgamma(P.star, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior2.st - l.prior2
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}

MH.delay1.ABC.exp.pois.hist <- function(P,S,rep, R_true, R_curr, tun, S_init, E_init, I_init, R_init, 
                                        beta_param, delayparam2, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0)) {
  count = 0
  repeat{
    u = rnorm(1, 0 ,tun[2])
    P.star = P + S*u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star>lowbnds[2]){
      break
    }
  }
  
  
  timespan0 = 0:maxT
  
  init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
  theta_ode <- c(beta = beta_param, s1 = P.star^(-1), s2 = delayparam2[2]^(-1), N = sum(init))
  sol = ode(y = init, times=timespan0, func=SEIR11, parms = theta_ode)
  R_cand = sol[, dim(sol)[2]]
  
  
  l_lik_curr = sum(dpois(R_true, lambda = R_curr, log = TRUE))
  l_lik_cand = sum(dpois(R_true, lambda = R_cand, log = TRUE))
  
  l.prior2.st <- dgamma(P.star, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior2.st - l.prior2
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}



MH.delay2.ABC.exp.hist <- function(P,S,rep, R_true, R_curr, sd_list, tun, S_init, E_init, I_init, R_init,
                                   beta_param, delayparam1, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), 
                                   noiseType, noise_multi, noise_const) {
  count = 0
  repeat{
    u = rnorm(1, 0 ,tun[2])
    P.star = P + S*u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  
  init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
  theta_ode <- c(beta = beta_param, s1 = delayparam1[2]^(-1), s2 = P.star^(-1), N = sum(init))
  sol = ode(y = init, times=timespan0, func=SEIR11, parms = theta_ode)
  R_cand = sol[, dim(sol)[2]]
  
  
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
  
  l_lik_curr = sum(dnorm(R_curr, mean = R_true, sd = sd_list_curr, log = TRUE))
  l_lik_cand = sum(dnorm(R_cand, mean = R_true, sd = sd_list_cand, log = TRUE))
  
  l.prior2.st <- dgamma(P.star, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior2.st - l.prior2
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}

MH.delay2.ABC.exp.pois.hist <- function(P,S,rep, R_true, R_curr, tun, S_init, E_init, I_init, R_init,
                                        beta_param, delayparam1, pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0)) {
  count = 0
  repeat{
    u = rnorm(1, 0 ,tun[2])
    P.star = P + S*u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star>lowbnds[2]){
      break
    }
  }
  
  timespan0 = 0:maxT
  
  init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
  theta_ode <- c(beta = beta_param, s1 = delayparam1[2]^(-1), s2 = P.star^(-1), N = sum(init))
  sol = ode(y = init, times=timespan0, func=SEIR11, parms = theta_ode)
  R_cand = sol[, dim(sol)[2]]
  
  l_lik_curr = sum(dpois(R_true, lambda = R_curr, log = TRUE))
  l_lik_cand = sum(dpois(R_true, lambda = R_cand, log = TRUE))
  
  l.prior2.st <- dgamma(P.star, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  l.prior2    <- dgamma(P     , shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l_lik_cand - l_lik_curr
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }else{
    logMH <- (l_lik_cand - l_lik_curr + l.prior2.st - l.prior2
              + log(pnorm(P - lowbnds[2], sd = S*sqrt(tun[2])))
              - log(pnorm(P.star - lowbnds[2], sd = S*sqrt(tun[2]))))
  }
  
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star
    count = 1
    R_curr = R_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count, R_new = R_curr))
}


MH.delay1.ABC.extended <- function(param_set, S, rep, cumul_true, cumul_curr, tun, y_init,
                                   pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), dd_lower,
                                   noise_multi, noise_const, dist_type){
  
  cumul_init = y_init$C_init
  count = 0
  delayparam = c(param_set$k_tauL, param_set$s_tauL)
  param_set_cand = param_set
  repeat{
    u = mvrnorm(1, c(0,0), diag(c(tun[1],tun[2])))
    delayparam.cand = delayparam + S%*%u
    if(delayparam.cand[1] > lowbnds[1] && delayparam.cand[2] > lowbnds[2]){
      break
    }
  }
  param_set_cand$k_tauL = delayparam.cand[1]
  param_set_cand$s_tauL = delayparam.cand[2]
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1] 
  dd = max(dd_lower, get_divisor_dist(c(param_set_cand$k_tauL, param_set_cand$s_tauL, param_set_cand$k_tauI, param_set_cand$s_tauI), dist_type = dist_type))
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_LARGE_dist_input(timespan = ext_tspan, param_set = param_set_cand, y_init = y_init, dist_type = dist_type)
  cumul_cand = sim_cand_ext$cumul_Ct[subseq0]
  
  sd_list_curr = noise_multi * sqrt(cumul_true - cumul_init) + noise_const
  sd_list_cand = sd_list_curr
  
  loglik_curr = sum(dnorm(cumul_curr, mean = cumul_true, sd = sd_list_curr, log = TRUE))
  loglik_cand = sum(dnorm(cumul_cand, mean = cumul_true, sd = sd_list_cand, log = TRUE))
  
  logpri_k_curr = dgamma(param_set$k_tauL, shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  logpri_k_cand = dgamma(param_set_cand$k_tauL, shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  logpri_s_curr = dgamma(param_set$s_tauL, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  logpri_s_cand = dgamma(param_set_cand$s_tauL, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (loglik_cand - loglik_curr
              + log(pmvnorm(upper = c(param_set$k_tauL - lowbnds[1], param_set$s_tauL - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(param_set_cand$k_tauL - lowbnds[1], param_set_cand$s_tauL - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (loglik_cand - loglik_curr + logpri_k_cand - logpri_k_curr + logpri_s_cand - logpri_s_curr
              + log(pmvnorm(upper = c(param_set$k_tauL - lowbnds[1], param_set$s_tauL - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(param_set_cand$k_tauL - lowbnds[1], param_set_cand$s_tauL - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min(exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    param_set = param_set_cand
    count = 1
    cumul_curr = cumul_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(param_set = param_set, S=S.update, count=count, cumul_new = cumul_curr))
}


MH.delay1.ABC.extended.pois <- function(param_set, S, rep, cumul_true, cumul_curr, tun, y_init,
                                        pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), dd_lower, dist_type){
  
  cumul_init = y_init$C_init
  count = 0
  delayparam = c(param_set$k_tauL, param_set$s_tauL)
  param_set_cand = param_set
  repeat{
    u = mvrnorm(1, c(0,0), diag(c(tun[1],tun[2])))
    delayparam.cand = delayparam + S%*%u
    if(delayparam.cand[1] > lowbnds[1] && delayparam.cand[2] > lowbnds[2]){
      break
    }
  }
  param_set_cand$k_tauL = delayparam.cand[1]
  param_set_cand$s_tauL = delayparam.cand[2]
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1] 
  dd = max(dd_lower, get_divisor_dist(c(param_set_cand$k_tauL, param_set_cand$s_tauL, param_set_cand$k_tauI, param_set_cand$s_tauI), dist_type = dist_type))
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_LARGE_dist_input(timespan = ext_tspan, param_set = param_set_cand, y_init = y_init, dist_type = dist_type)
  cumul_cand = sim_cand_ext$cumul_Ct[subseq0]
  
  loglik_curr = sum(dpois(cumul_true, lambda = cumul_curr, log = TRUE))
  loglik_cand = sum(dpois(cumul_true, lambda = cumul_cand, log = TRUE))
  
  logpri_k_curr = dgamma(param_set$k_tauL, shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  logpri_k_cand = dgamma(param_set_cand$k_tauL, shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  logpri_s_curr = dgamma(param_set$s_tauL, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  logpri_s_cand = dgamma(param_set_cand$s_tauL, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (loglik_cand - loglik_curr
              + log(pmvnorm(upper = c(param_set$k_tauL - lowbnds[1], param_set$s_tauL - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(param_set_cand$k_tauL - lowbnds[1], param_set_cand$s_tauL - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (loglik_cand - loglik_curr + logpri_k_cand - logpri_k_curr + logpri_s_cand - logpri_s_curr
              + log(pmvnorm(upper = c(param_set$k_tauL - lowbnds[1], param_set$s_tauL - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(param_set_cand$k_tauL - lowbnds[1], param_set_cand$s_tauL - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min(exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    param_set = param_set_cand
    count = 1
    cumul_curr = cumul_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(param_set = param_set, S=S.update, count=count, cumul_new = cumul_curr))
}


MH.delay2.ABC.extended <- function(param_set, S, rep, cumul_true, cumul_curr, tun, y_init,
                                   pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), dd_lower,
                                   noise_multi, noise_const, dist_type){
  
  cumul_init = y_init$C_init
  count = 0
  delayparam = c(param_set$k_tauI, param_set$s_tauI)
  param_set_cand = param_set
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    delayparam.cand = delayparam + S%*%u
    if(delayparam.cand[1] > lowbnds[1] && delayparam.cand[2] > lowbnds[2]){
      break
    }
  }
  param_set_cand$k_tauI = delayparam.cand[1]
  param_set_cand$s_tauI = delayparam.cand[2]
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1] 
  dd = max(dd_lower, get_divisor_dist(c(param_set_cand$k_tauL, param_set_cand$s_tauL, param_set_cand$k_tauI, param_set_cand$s_tauI), dist_type = dist_type))
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_LARGE_dist_input(timespan = ext_tspan, param_set = param_set_cand, y_init = y_init, dist_type = dist_type)
  cumul_cand = sim_cand_ext$cumul_Ct[subseq0]
  
  sd_list_curr = noise_multi * sqrt(cumul_true - cumul_init) + noise_const
  sd_list_cand = sd_list_curr
  
  loglik_curr = sum(dnorm(cumul_curr, mean = cumul_true, sd = sd_list_curr, log = TRUE))
  loglik_cand = sum(dnorm(cumul_cand, mean = cumul_true, sd = sd_list_cand, log = TRUE))
  
  logpri_k_curr = dgamma(param_set$k_tauI, shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  logpri_k_cand = dgamma(param_set_cand$k_tauI, shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  logpri_s_curr = dgamma(param_set$s_tauI, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  logpri_s_cand = dgamma(param_set_cand$s_tauI, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (loglik_cand - loglik_curr
              + log(pmvnorm(upper = c(param_set$k_tauI - lowbnds[1], param_set$s_tauI - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(param_set_cand$k_tauI - lowbnds[1], param_set_cand$s_tauI - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (loglik_cand - loglik_curr + logpri_k_cand - logpri_k_curr + logpri_s_cand - logpri_s_curr
              + log(pmvnorm(upper = c(param_set$k_tauI - lowbnds[1], param_set$s_tauI - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(param_set_cand$k_tauI - lowbnds[1], param_set_cand$s_tauI - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min(exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    param_set = param_set_cand
    count = 1
    cumul_curr = cumul_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(param_set = param_set, S=S.update, count=count, cumul_new = cumul_curr))
}

MH.delay2.ABC.extended.pois <- function(param_set, S, rep, cumul_true, cumul_curr, tun, y_init,
                                        pri.shape, pri.scale, maxT, flatpri = FALSE, lowbnds = c(0,0), dd_lower, dist_type){
  
  cumul_init = y_init$C_init
  count = 0
  delayparam = c(param_set$k_tauI, param_set$s_tauI)
  param_set_cand = param_set
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    delayparam.cand = delayparam + S%*%u
    if(delayparam.cand[1] > lowbnds[1] && delayparam.cand[2] > lowbnds[2]){
      break
    }
  }
  param_set_cand$k_tauI = delayparam.cand[1]
  param_set_cand$s_tauI = delayparam.cand[2]
  
  timespan0 = 0:maxT
  h = timespan0[2] - timespan0[1] 
  dd = max(dd_lower, get_divisor_dist(c(param_set_cand$k_tauL, param_set_cand$s_tauL, param_set_cand$k_tauI, param_set_cand$s_tauI), dist_type = dist_type))
  ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
  subseq0 = round(seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd))
  
  sim_cand_ext = mean_trajectory_LARGE_dist_input(timespan = ext_tspan, param_set = param_set_cand, y_init = y_init, dist_type = dist_type)
  cumul_cand = sim_cand_ext$cumul_Ct[subseq0]
  
  loglik_curr = sum(dpois(cumul_true, lambda = cumul_curr, log = TRUE))
  loglik_cand = sum(dpois(cumul_true, lambda = cumul_cand, log = TRUE))
  
  logpri_k_curr = dgamma(param_set$k_tauI, shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  logpri_k_cand = dgamma(param_set_cand$k_tauI, shape = pri.shape[1], rate = pri.shape[2], log = TRUE)
  logpri_s_curr = dgamma(param_set$s_tauI, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  logpri_s_cand = dgamma(param_set_cand$s_tauI, shape = pri.scale[1], rate = pri.scale[2], log = TRUE)
  
  if(flatpri){
    logMH <- (loglik_cand - loglik_curr
              + log(pmvnorm(upper = c(param_set$k_tauI - lowbnds[1], param_set$s_tauI - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(param_set_cand$k_tauI - lowbnds[1], param_set_cand$s_tauI - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (loglik_cand - loglik_curr + logpri_k_cand - logpri_k_curr + logpri_s_cand - logpri_s_curr
              + log(pmvnorm(upper = c(param_set$k_tauI - lowbnds[1], param_set$s_tauI - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(param_set_cand$k_tauI - lowbnds[1], param_set_cand$s_tauI - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  alpha = min(exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    param_set = param_set_cand
    count = 1
    cumul_curr = cumul_cand
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-1)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(param_set = param_set, S=S.update, count=count, cumul_new = cumul_curr))
}
