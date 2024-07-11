Gillespie_SEIR_delayed_GI <- function(S_init = 30, E_init = 0, I_init = 1, R_init = 0, beta.param, 
                                      maxT = 10, time.interval = 1, delayparam1 = c(0.3,1), delayparam2 = c(0.3,1),
                                      absorbing = TRUE){
  # Stoichiometric Matrix.
  # S_init, I_init:, R_init: the initial values for the numbers of S, I, and R.
  # beta.param: infectious rate (beta.paramSI) 
  # maxT: the maximum time of observation 
  # time.interval: interval between observation time.
  
  StoiMatrix = matrix(data = c(-1,1,0,0,0,-1,1,0,0,0,-1,1), nrow = 4, ncol = 3)
  SpecMatrix = as.matrix(c(S_init, E_init, I_init, R_init)) # The vector representing the number of the species.
  SpecMatrix.discrete = matrix(SpecMatrix, nrow = 4, ncol = floor(maxT/time.interval+1))
  timevector.discrete = seq(from = 0, by = time.interval, to = maxT)
  SpecMatrix.discrete[,1] = c(S_init, E_init, I_init, R_init)
  
  #Infectious_time_matrix = matrix(data = c(0, 1, 0), nrow = 3, ncol = 1)
  
  if(I_init > 0){
    I_labels = 1:I_init  
  }else{
    I_labels = c()
  }
  if(E_init > 0){
    E_labels = (I_init+1):(I_init + E_init)
  }else{
    E_labels = c()
  }
  if(E_init > 0){
    R_labels = (I_init + E_init + 1):(I_init + E_init + R_init)
  }else{
    R_labels = c()
  }
  
  Curr_S_label = (I_init + E_init + R_init) + 1
  Ntot = S_init + E_init + I_init + R_init # the total number of population.
  
  
  Exposed_time_matrix = matrix(data = 0, nrow = 3, ncol = (E_init + I_init))
  
  for(ii in 1:I_init){
    Exposed_time_matrix[2,ii] = ii
    #Exposed_time_matrix[3,ii] = -rgamma(n = 1, shape = delayparam1[1], scale = delayparam1[2])
  }
  if(E_init > 0){
    for(ii in (I_init+1):(I_init+E_init)){
      Exposed_time_matrix[2,ii] = ii
    }
  }
  timevector = c(0)
  
  stackTime1 <- rgamma(n = E_init, shape = delayparam1[1], scale = delayparam1[2])
  stackTime2 <- rgamma(n = I_init, shape = delayparam2[1], scale = delayparam2[2])
  
  current.time= 0 
  
  while (current.time < maxT){
    
    propensities = c(beta.param*SpecMatrix[1,ncol(SpecMatrix)]*SpecMatrix[3, ncol(SpecMatrix)] / Ntot)
    # cat(propensities)
    if (sum(SpecMatrix[1:3, ncol(SpecMatrix)]) == 0){
      # cat("BREAK!!")
      break
    }
    
    #stackTime1 = sort(stackTime1)
    #stackTime2 = sort(stackTime2)
    
    if(sum(propensities) == 0){
      time.gap = 0
    }else{
      time.gap = rexp(n = 1, rate = sum(propensities))
    }
    current.time = current.time + time.gap
    
    if(!(length(stackTime1) == 0 & length(stackTime2) == 0)){
      minStack <- min(stackTime1, stackTime2)
    }else{
      minStack <- Inf
    }
    
    # cat(current.time)
    # cat("\n")
    
    if (current.time >= maxT){
      break
    }
    
    if (current.time < minStack & sum(propensities) != 0){
      # contact S+I -> E+I
      # current.reaction.index = sum(cumsum(propensities)/sum(propensities) < r1) + 1
      current.reaction.index = 1 # the chosen event must always be S -> E reaction because the other two reactions can occur only by completing delays. 
      SpecMatrix = cbind(SpecMatrix, SpecMatrix[, ncol(SpecMatrix)] + StoiMatrix[,1])
      timevector = c(timevector, current.time)  
      stackTime1 = c(stackTime1, current.time + rgamma(n=1, shape = delayparam1[1], scale = delayparam1[2]))
      
      E_labels = c(E_labels, Curr_S_label)
      # if(length(I_labels) > 0){
      Infector_id = sample(I_labels, 1)
      # }
      Exposed_time_matrix = cbind(Exposed_time_matrix, c(Infector_id, Curr_S_label, current.time))
      
      Curr_S_label = Curr_S_label + 1
      
    }else if(minStack <= maxT){
      if (min(stackTime1) < min(stackTime2)){
        # transition: E -> I
        SpecMatrix = cbind(SpecMatrix, SpecMatrix[, ncol(SpecMatrix)] + StoiMatrix[, 2])
        timevector = c(timevector, minStack)  
        current.time = minStack
        
        EtoI_id = which(stackTime1 == min(stackTime1))
        I_labels = c(I_labels, E_labels[EtoI_id])
        E_labels = E_labels[-EtoI_id]
        
        #stackTime1 = stackTime1[-1]
        stackTime1 = stackTime1[-EtoI_id]
        stackTime2 = c(stackTime2, current.time + rgamma(n=1, shape = delayparam2[1], scale = delayparam2[2]))
      }else{
        # transition: I -> R
        SpecMatrix = cbind(SpecMatrix, SpecMatrix[, ncol(SpecMatrix)] + StoiMatrix[, 3])
        timevector = c(timevector, minStack)  
        current.time = minStack 
        
        ItoR_id = which(stackTime2 == min(stackTime2))
        R_labels = c(R_labels, R_labels[ItoR_id])
        I_labels = I_labels[-ItoR_id]
        
        stackTime2 = stackTime2[-ItoR_id]
        #stackTime2 = stackTime2[-1]
      }
    }else{ # this condition represents that minStack, current.time > maxT or sum(propensities) == 0
      if (absorbing == TRUE | SpecMatrix[1,ncol(SpecMatrix)] == 0){
        break
      }else{
        # cat("Avoid absorbing!")
        SpecMatrix = cbind(SpecMatrix, SpecMatrix[, ncol(SpecMatrix)] + StoiMatrix[,1])
        timevector = c(timevector, current.time) 
        stackTime1 = c(stackTime1, current.time + rgamma(n=1, shape = delayparam1[1], scale = delayparam1[2]))
      }
    }
    SpecMatrix.discrete[, (ceiling(current.time/time.interval)+1):floor(maxT/time.interval + 1)] = SpecMatrix[, ncol(SpecMatrix)]
  }
  
  pairmat = Exposed_time_matrix
  
  pairmat[3,] = floor(pairmat[3,])
  total_exp_num = dim(pairmat)[2]
  
  GI_mat = matrix(NA, nrow = 4, ncol = total_exp_num)
  GI_mat[1:3, 1:total_exp_num] = pairmat
  
  for(exposure_id in 1:total_exp_num){
    infector = pairmat[1, exposure_id]
    infectee = pairmat[2, exposure_id]
    secondary_time = pairmat[3, exposure_id]
    if(infector == 0){
      primary_time = NA
    }else{
      primary_id = which(infector == pairmat[2, ])
      primary_time = pairmat[3, primary_id]  
    }
    GI_mat[4, exposure_id] = secondary_time - primary_time
  }
  
  result.list = list("timevector" = timevector, "SpecMatrix" = SpecMatrix, 
                     "SpecMatrix.discrete" = SpecMatrix.discrete, 
                     "timevector.discrete" = timevector.discrete,
                     "Gen_time_mat" = GI_mat)
  return(result.list)
  
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

S_init = 1000
E_init = 0
I_init = 1
R_init = 0
beta.param = 0.6
maxT = 40

delayparam1 = c(5, 1.2)
delayparam2 = c(5, 1.2)

output1 = Gillespie_SEIR_delayed_GI(S_init = S_init, E_init = E_init, I_init = I_init, R_init = R_init, beta.param = beta.param, 
                                    maxT = maxT, time.interval = 1, delayparam1 = delayparam1, 
                                    delayparam2 = delayparam2, absorbing = F)
GI_mat = output1$Gen_time_mat
# Every column of GI_mat represents one infection.
# The 1st row of GI_mat is the infector ID.
# The 2nd row of GI_mat is the infectee ID.
# The 3rd row of GI_mat is the time of infection.
# The 4th row of GI_mat is the corresponding generation time.

