library(deSolve)
library(MASS)
library(mvtnorm)
library(dplyr)
library(invgamma)
library(ramcmc)

SEIR11 <- function(t, s, theta) {
  with(as.list(c(s, theta)), {
    dS <- -beta*S*I1/N
    dE1 <- beta*S*I1/N - s1 * E1
    dI1 <- s1 * E1 - s2 * I1
    dR <- s2 * I1
    return(list(c(dS, dE1, dI1, dR)))
  })
}

mean_trajectory_SEIR <- function(timespan, theta, y_init){
  beta.param = theta[1]
  alpha1 = theta[2]
  beta1 = theta[3]
  alpha2 = theta[4]
  beta2 = theta[5]
  
  Ntot = sum(y_init)
  
  h = timespan[2] - timespan[1]
  St = rep(NA, length(timespan))
  Et = rep(NA, length(timespan))
  It = rep(NA, length(timespan))
  Rt = rep(NA, length(timespan))
  conv_val = rep(NA, length(timespan))
  
  # hs = rep(NA, length(timespan))
  # hs_tilde = rep(NA, length(timespan))
  St[1] = y_init[1]
  Et[1] = y_init[2]
  It[1] = y_init[3]
  Rt[1] = y_init[4]
  
  for(ii in 1:length(timespan)){
    conv_val[ii] = conv_two_gamma(tmax = timespan[ii], alpha1, beta1, alpha2, beta2)
  }
  
  for(ii in 1:(length(timespan)-1)){
    
    term1 = beta.param * St[ii] * It[ii] / Ntot
    term2 = beta.param * St[1] * It[1] / Ntot * pgamma(timespan[ii], shape = alpha1, scale = beta1) + 
      numer_int(timespan[1:ii], beta.param * diff_midpoint_extension(St[1:ii] * It[1:ii], h) / Ntot * pgamma(max(timespan[1:ii]) - timespan[1:ii], shape = alpha1, scale = beta1)) + 
      Et[1] * dgamma(timespan[ii], shape = alpha1, scale = beta1)
    term3 = numer_int(timespan[1:ii], beta.param * St[1:ii] * It[1:ii] / Ntot * conv_val[ii:1]) +  It[1] * dgamma(timespan[ii], shape = alpha2, scale = beta2) + Et[1] * conv_val[ii]
    
    term1 = max(term1, 0) # adjusted to avoid negative values 
    term2 = max(term2, 0) # adjusted to avoid negative values 
    term3 = max(term3, 0) # adjusted to avoid negative values 
    
    S_increm = - term1
    E_increm = term1 - term2
    I_increm = term2 - term3
    R_increm = term3
    
    St_tilde = St[ii] + h*S_increm
    Et_tilde = Et[ii] + h*E_increm
    It_tilde = It[ii] + h*I_increm
    Rt_tilde = Rt[ii] + h*R_increm
    
    term1_tilde = beta.param * St_tilde * It_tilde / Ntot
    term2_tilde = beta.param * St[1] * It[1] / Ntot * pgamma(timespan[ii+1], shape = alpha1, scale = beta1) + 
      numer_int(timespan[1:(ii+1)], beta.param * diff_midpoint_extension(c(St[1:ii], St_tilde) * c(It[1:ii], It_tilde), h) / Ntot * pgamma(max(timespan[1:(ii+1)]) - timespan[1:(ii+1)], shape = alpha1, scale = beta1)) + 
      Et[1] * dgamma(timespan[ii+1], shape = alpha1, scale = beta1)
    term3_tilde = numer_int(timespan[1:(ii+1)], beta.param * c(St[1:ii], St_tilde) * c(It[1:ii], It_tilde) / Ntot * conv_val[(ii+1):1])  + It[1] * dgamma(timespan[ii+1], shape = alpha2, scale = beta2) + Et[1] * conv_val[ii+1]
    
    term1_tilde = max(term1_tilde, 0)
    term2_tilde = max(term2_tilde, 0)
    term3_tilde = max(term3_tilde, 0) # adjusted to avoid negative values 
    
    S_increm_tilde = - term1_tilde
    E_increm_tilde = term1_tilde - term2_tilde
    I_increm_tilde = term2_tilde - term3_tilde
    R_increm_tilde = term3_tilde
    
    St[ii+1] = St[ii] + h/2 * (S_increm + S_increm_tilde)
    Et[ii+1] = Et[ii] + h/2 * (E_increm + E_increm_tilde)
    It[ii+1] = It[ii] + h/2 * (I_increm + I_increm_tilde)
    Rt[ii+1] = Rt[ii] + h/2 * (R_increm + R_increm_tilde)
    
  }
  return(list("St" = St, "Et" = Et, "It" = It, "Rt" = Rt))
}


mean_trajectory_SEIR_dist_input <- function(timespan, theta, y_init, dist_type){
  beta.param = theta[1]
  alpha1 = theta[2]
  beta1 = theta[3]
  alpha2 = theta[4]
  beta2 = theta[5]
  
  Ntot = sum(y_init)
  
  h = timespan[2] - timespan[1]
  St = rep(NA, length(timespan))
  Et = rep(NA, length(timespan))
  It = rep(NA, length(timespan))
  Rt = rep(NA, length(timespan))
  conv_val = rep(NA, length(timespan))
  
  # hs = rep(NA, length(timespan))
  # hs_tilde = rep(NA, length(timespan))
  St[1] = y_init[1]
  Et[1] = y_init[2]
  It[1] = y_init[3]
  Rt[1] = y_init[4]
  
  for(ii in 1:length(timespan)){
    conv_val[ii] = conv_two_pdfs(tmax = timespan[ii], theta1 = c(alpha1, beta1), theta2 = c(alpha2, beta2), dist_type = dist_type)
  }
  
  for(ii in 1:(length(timespan)-1)){
    
    term1 = beta.param * St[ii] * It[ii] / Ntot
    term2 = beta.param * St[1] * It[1] / Ntot * eval_density(x = timespan[ii], theta = c(alpha1, beta1), dist_type = dist_type, density_type = "cdf") + 
      numer_int(timespan[1:ii], beta.param * diff_midpoint_extension(St[1:ii] * It[1:ii], h) / Ntot * eval_density(x = max(timespan[1:ii]) - timespan[1:ii], theta = c(alpha1, beta1), dist_type = dist_type, density_type = "cdf")) + 
      Et[1] * eval_density(x = timespan[ii], theta = c(alpha1, beta1), dist_type = dist_type, density_type = "pdf")
    term3 = numer_int(timespan[1:ii], beta.param * St[1:ii] * It[1:ii] / Ntot * conv_val[ii:1]) +  It[1] * eval_density(x = timespan[ii], theta = c(alpha2, beta2), dist_type = dist_type, density_type = "pdf") + Et[1] * conv_val[ii]
    
    term1 = max(term1, 0) # adjusted to avoid negative values 
    term2 = max(term2, 0) # adjusted to avoid negative values 
    term3 = max(term3, 0) # adjusted to avoid negative values 
    
    S_increm = - term1
    E_increm = term1 - term2
    I_increm = term2 - term3
    R_increm = term3
    
    St_tilde = St[ii] + h*S_increm
    Et_tilde = Et[ii] + h*E_increm
    It_tilde = It[ii] + h*I_increm
    Rt_tilde = Rt[ii] + h*R_increm
    
    term1_tilde = beta.param * St_tilde * It_tilde / Ntot
    term2_tilde = beta.param * St[1] * It[1] / Ntot * eval_density(x = timespan[ii+1], theta = c(alpha1, beta1), dist_type = dist_type, density_type = "cdf") + 
      numer_int(timespan[1:(ii+1)], beta.param * diff_midpoint_extension(c(St[1:ii], St_tilde) * c(It[1:ii], It_tilde), h) / Ntot * eval_density(x = max(timespan[1:(ii+1)]) - timespan[1:(ii+1)], theta = c(alpha1, beta1), dist_type = dist_type, density_type = "cdf")) + 
      Et[1] * eval_density(x = timespan[ii+1], theta = c(alpha1, beta1), dist_type = dist_type, density_type = "pdf")
    term3_tilde = numer_int(timespan[1:(ii+1)], beta.param * c(St[1:ii], St_tilde) * c(It[1:ii], It_tilde) / Ntot * conv_val[(ii+1):1])  + It[1] * eval_density(x = timespan[ii+1], theta = c(alpha2, beta2), dist_type = dist_type, density_type = "pdf") + Et[1] * conv_val[ii+1]
    
    term1_tilde = max(term1_tilde, 0) # adjusted to avoid negative values 
    term2_tilde = max(term2_tilde, 0) # adjusted to avoid negative values 
    term3_tilde = max(term3_tilde, 0) # adjusted to avoid negative values 
    
    S_increm_tilde = - term1_tilde
    E_increm_tilde = term1_tilde - term2_tilde
    I_increm_tilde = term2_tilde - term3_tilde
    R_increm_tilde = term3_tilde
    
    St[ii+1] = St[ii] + h/2 * (S_increm + S_increm_tilde)
    Et[ii+1] = Et[ii] + h/2 * (E_increm + E_increm_tilde)
    It[ii+1] = It[ii] + h/2 * (I_increm + I_increm_tilde)
    Rt[ii+1] = Rt[ii] + h/2 * (R_increm + R_increm_tilde)
    
  }
  return(list("St" = St, "Et" = Et, "It" = It, "Rt" = Rt))
}


mean_trajectory_LARGE_dist_input <- function(timespan, param_set, y_init, dist_type = "gamma"){
  beta = param_set$beta
  h_var = param_set$h_var
  e_asym = param_set$e_asym
  e_vacc = param_set$e_vacc
  vacc = param_set$vacc
  rec = param_set$rec
  death = param_set$death
  p_A = param_set$p_A
  
  k_tauL = param_set$k_tauL
  s_tauL = param_set$s_tauL
  k_tauLh = param_set$k_tauLh
  s_tauLh = param_set$s_tauLh
  k_tauI = param_set$k_tauI
  s_tauI = param_set$s_tauI
  k_tauIA = param_set$k_tauIA
  s_tauIA = param_set$s_tauIA
  
  Stoi_mat = matrix(nrow = 11, ncol = 17, byrow = TRUE, data = c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,
                                                                 0,0,0,0,+1,-1,-1,0,0,0,0,0,0,0,0,0,0,
                                                                 +1,+1,0,0,0,+1,0,0,0,-1,-1,0,0,0,0,0,0,
                                                                 0,0,+1,+1,0,0,+1,0,0,0,0,-1,-1,0,0,0,0,
                                                                 0,0,0,0,0,0,0,0,0,+1,0,0,0,-1,0,0,0,
                                                                 0,0,0,0,0,0,0,0,0,0,0,+1,0,0,-1,0,0,
                                                                 0,0,0,0,0,0,0,0,0,0,+1,0,0,0,0,-1,0,
                                                                 0,0,0,0,0,0,0,0,0,0,0,0,+1,0,0,0,-1,
                                                                 0,0,0,0,0,0,0,-1,-1,0,0,0,0,+1,+1,0,0,
                                                                 0,0,0,0,0,0,0,+1,0,0,0,0,0,0,0,+1,+1,
                                                                 0,0,0,0,0,0,0,0,+1,0,0,0,0,0,0,0,0))
  
  h = timespan[2] - timespan[1]
  St = rep(NA, length(timespan))
  Vt = rep(NA, length(timespan))
  Et = rep(NA, length(timespan))
  Eht = rep(NA, length(timespan))
  It = rep(NA, length(timespan))
  Iht = rep(NA, length(timespan))
  At = rep(NA, length(timespan))
  Aht = rep(NA, length(timespan))
  Ct = rep(NA, length(timespan))
  Rt = rep(NA, length(timespan))
  Dt = rep(NA, length(timespan))
  Nt = rep(NA, length(timespan))
  cumul_Ct = rep(NA, length(timespan))
  cumul_asympt = rep(NA, length(timespan))
  
  St[1] = y_init$S_init
  Vt[1] = y_init$V_init
  Et[1] = y_init$E_init
  Eht[1] = y_init$Eh_init
  It[1] = y_init$I_init
  Iht[1] = y_init$Ih_init
  At[1] = y_init$A_init
  Aht[1] = y_init$Ah_init
  Ct[1] = y_init$C_init
  Rt[1] = y_init$R_init
  Dt[1] = y_init$D_init
  Nt[1] = sum(St[1], Vt[1], Et[1], Eht[1], It[1], Iht[1], At[1], Aht[1], Ct[1], Rt[1], Dt[1]) - Dt[1]
  cumul_Ct[1] = Ct[1]
  cumul_asympt[1] = At[1] + Aht[1]
  
  conv_tauL_tauI = rep(NA, length(timespan))
  conv_tauL_tauIA = rep(NA, length(timespan))
  conv_tauLh_tauI = rep(NA, length(timespan))
  conv_tauLh_tauIA = rep(NA, length(timespan))
  
  for(ii in 1:length(timespan)){
    conv_tauL_tauI[ii] = conv_two_pdfs(tmax = timespan[ii], theta1 = c(k_tauL, s_tauL), theta2 = c(k_tauI, s_tauI), dist_type = dist_type)
    conv_tauL_tauIA[ii] = conv_two_pdfs(tmax = timespan[ii], theta1 = c(k_tauL, s_tauL), theta2 = c(k_tauIA, s_tauIA), dist_type = dist_type)
    conv_tauLh_tauI[ii] = conv_two_pdfs(tmax = timespan[ii], theta1 = c(k_tauLh, s_tauLh), theta2 = c(k_tauI, s_tauI), dist_type = dist_type)
    conv_tauLh_tauIA[ii] = conv_two_pdfs(tmax = timespan[ii], theta1 = c(k_tauLh, s_tauLh), theta2 = c(k_tauIA, s_tauIA), dist_type = dist_type)
  }
  
  for(ii in 1:(length(timespan)-1)){
    
    term1 = beta * St[ii] * It[ii] / Nt[ii]
    term2 = e_asym * beta * St[ii] * At[ii] / Nt[ii]
    term3 = (1+h_var) * beta * St[ii] * Iht[ii] / Nt[ii]
    term4 = (1+h_var) * e_asym * beta * St[ii] * Aht[ii] / Nt[ii]
    term5 = vacc * St[ii]
    term6 = e_vacc * beta * Vt[ii] * It[ii] / Nt[ii]
    term7 = (1+h_var) * e_vacc * beta * Vt[ii] * Iht[ii] / Nt[ii]
    term8 = rec * Ct[ii]
    term9 = death * Ct[ii]
    
    fromE = (beta / Nt[1]) * (St[1] * It[1] + e_vacc * Vt[1] * It[1] + e_asym * St[1] * At[1]) * eval_density(x = timespan[ii], theta = c(k_tauL, s_tauL), dist_type = dist_type, density_type = "cdf") + 
      numer_int(xspan = timespan[1:ii], yval = diff_midpoint_extension(x = (beta / Nt[1:ii]) * (St[1:ii] * It[1:ii] + e_vacc * Vt[1:ii] * It[1:ii] + e_asym * St[1:ii] * At[1:ii]), t_interv = h) * eval_density(x = max(timespan[1:ii]) - timespan[1:ii], theta = c(k_tauL, s_tauL), dist_type = dist_type, density_type = "cdf")) +
      Et[1] * eval_density(x = timespan[ii], theta = c(k_tauL, s_tauL), dist_type = dist_type, density_type = "pdf")
    fromEh = ((1+h_var) * beta / Nt[1]) * (St[1] * Iht[1] + e_vacc * Vt[1] * Iht[1] + e_asym * St[1] * Aht[1]) * eval_density(x = timespan[ii], theta = c(k_tauLh, s_tauLh), dist_type = dist_type, density_type = "cdf") + 
      numer_int(xspan = timespan[1:ii], yval = diff_midpoint_extension(x = ((1+h_var) * beta / Nt[1:ii]) * (St[1:ii] * Iht[1:ii] + e_vacc * Vt[1:ii] * Iht[1:ii] + e_asym * St[1:ii] * Aht[1:ii]), t_interv = h) * eval_density(x = max(timespan[1:ii]) - timespan[1:ii], theta = c(k_tauLh, s_tauLh), dist_type = dist_type, density_type = "cdf")) +
      Eht[1] * eval_density(x = timespan[ii], theta = c(k_tauL, s_tauL), dist_type = dist_type, density_type = "pdf")
    fromItoC = (1-p_A) * (numer_int(xspan = timespan[1:ii], yval = (beta / Nt[1:ii]) * (St[1:ii] * It[1:ii] + e_vacc * Vt[1:ii] * It[1:ii] + e_asym * St[1:ii] * At[1:ii]) * conv_tauL_tauI[ii:1]) + Et[1] * conv_tauL_tauI[ii]) + 
      It[1] * eval_density(x = timespan[ii], theta = c(k_tauI, s_tauI), dist_type = dist_type, density_type = "pdf")
    fromIhtoC = (1-p_A) * (numer_int(xspan = timespan[1:ii], yval = ((1+h_var) * beta / Nt[1:ii]) * (St[1:ii] * Iht[1:ii] + e_vacc * Vt[1:ii] * Iht[1:ii] + e_asym * St[1:ii] * Aht[1:ii]) * conv_tauLh_tauI[ii:1]) + Eht[1] * conv_tauLh_tauI[ii]) + 
      Iht[1] * eval_density(x = timespan[ii], theta = c(k_tauI, s_tauI), dist_type = dist_type, density_type = "pdf")
    fromAtoR = p_A * (numer_int(xspan = timespan[1:ii], yval = (beta / Nt[1:ii]) * (St[1:ii] * It[1:ii] + e_vacc * Vt[1:ii] * It[1:ii] + e_asym * St[1:ii] * At[1:ii]) * conv_tauL_tauIA[ii:1]) + Et[1] * conv_tauL_tauIA[ii]) + 
      At[1] * eval_density(x = timespan[ii], theta = c(k_tauIA, s_tauIA), dist_type = dist_type, density_type = "pdf")
    fromAhtoR = p_A * (numer_int(xspan = timespan[1:ii], yval = ((1+h_var) * beta / Nt[1:ii]) * (St[1:ii] * Iht[1:ii] + e_vacc * Vt[1:ii] * Iht[1:ii] + e_asym * St[1:ii] * Aht[1:ii]) * conv_tauLh_tauIA[ii:1]) + Eht[1] * conv_tauLh_tauIA[ii]) + 
      Aht[1] * eval_density(x = timespan[ii], theta = c(k_tauIA, s_tauIA), dist_type = dist_type, density_type = "pdf")
    
    term10 = (1-p_A) * fromE
    term11 = p_A * fromE
    term12 = (1-p_A) * fromEh
    term13 = p_A * fromEh
    term14 = fromItoC
    term15 = fromIhtoC
    term16 = fromAtoR
    term17 = fromAhtoR
    
    term1 = max(term1, 0) # adjusted to avoid negative values 
    term2 = max(term2, 0) # adjusted to avoid negative values 
    term3 = max(term3, 0) # adjusted to avoid negative values 
    term4 = max(term4, 0) # adjusted to avoid negative values 
    term5 = max(term5, 0) # adjusted to avoid negative values 
    term6 = max(term6, 0) # adjusted to avoid negative values 
    term7 = max(term7, 0) # adjusted to avoid negative values 
    term8 = max(term8, 0) # adjusted to avoid negative values 
    term9 = max(term9, 0) # adjusted to avoid negative values 
    term10 = max(term10, 0) # adjusted to avoid negative values 
    term11 = max(term11, 0) # adjusted to avoid negative values 
    term12 = max(term12, 0) # adjusted to avoid negative values 
    term13 = max(term13, 0) # adjusted to avoid negative values 
    term14 = max(term14, 0) # adjusted to avoid negative values 
    term15 = max(term15, 0) # adjusted to avoid negative values 
    term16 = max(term16, 0) # adjusted to avoid negative values 
    term17 = max(term17, 0) # adjusted to avoid negative values 
    
    
    increm_mat = Stoi_mat %*% c(term1, 
                                term2, 
                                term3, 
                                term4, 
                                term5, 
                                term6, 
                                term7, 
                                term8, 
                                term9, 
                                term10, 
                                term11, 
                                term12, 
                                term13, 
                                term14, 
                                term15, 
                                term16, 
                                term17)
    
    S_increm = increm_mat[1]
    V_increm = increm_mat[2]
    E_increm = increm_mat[3]
    Eh_increm = increm_mat[4]
    I_increm = increm_mat[5]
    Ih_increm = increm_mat[6]
    A_increm = increm_mat[7]
    Ah_increm = increm_mat[8]
    C_increm = increm_mat[9]
    R_increm = increm_mat[10]
    D_increm = increm_mat[11]
    cumul_C_increm = + term14 + term15
    cumul_asymp_increm = + term16 + term17
    
    # S_increm = - term1 - term2 - term3 - term4 - term5
    # V_increm = + term5 - term6 - term7
    # E_increm = + term1 + term2 + term6 - term10 - term11
    # Eh_increm = + term3 + term4 + term7 - term12 - term13
    # I_increm = + term10 - term14
    # Ih_increm = + term12 - term15
    # A_increm = + term11 - term16
    # Ah_increm = + term13 - term17
    # C_increm = - term8 - term9 + term14 + term15
    # R_increm = + term8 + term16 + term17
    # D_increm = + term9
    
    
    S_tilde = St[ii] + h*S_increm
    V_tilde = Vt[ii] + h*V_increm
    E_tilde = Et[ii] + h*E_increm
    Eh_tilde = Eht[ii] + h*Eh_increm
    I_tilde = It[ii] + h*I_increm
    Ih_tilde = Iht[ii] + h*Ih_increm
    A_tilde = At[ii] + h*A_increm
    Ah_tilde = Aht[ii] + h*Ah_increm
    C_tilde = Ct[ii] + h*C_increm
    R_tilde = Rt[ii] + h*R_increm
    D_tilde = Dt[ii] + h*D_increm
    N_tilde = sum(S_tilde, V_tilde, E_tilde, Eh_tilde, I_tilde, Ih_tilde, A_tilde, Ah_tilde, C_tilde, R_tilde, D_tilde) - D_tilde
    
    St_ext = c(St[1:ii], S_tilde)
    Vt_ext = c(Vt[1:ii], V_tilde)
    Et_ext = c(Et[1:ii], E_tilde)
    Eht_ext = c(Eht[1:ii], Eh_tilde)
    It_ext = c(It[1:ii], I_tilde)
    Iht_ext = c(Iht[1:ii], Ih_tilde)
    At_ext = c(At[1:ii], A_tilde)
    Aht_ext = c(Aht[1:ii], Ah_tilde)
    Ct_ext = c(Ct[1:ii], C_tilde)
    Rt_ext = c(Rt[1:ii], R_tilde)
    Dt_ext = c(Dt[1:ii], D_tilde)
    Nt_ext = c(Nt[1:ii], N_tilde)
    
    term1_tilde = beta * S_tilde * I_tilde / N_tilde
    term2_tilde = e_asym * beta * S_tilde * A_tilde / N_tilde
    term3_tilde = (1+h_var) * beta * S_tilde * Ih_tilde / N_tilde
    term4_tilde = (1+h_var) * e_asym * beta * S_tilde * Ah_tilde / N_tilde
    term5_tilde = vacc * S_tilde
    term6_tilde = e_vacc * beta * V_tilde * I_tilde / N_tilde
    term7_tilde = (1+h_var) * e_vacc * beta * V_tilde * Ih_tilde / N_tilde
    term8_tilde = rec * C_tilde
    term9_tilde = death * C_tilde
    
    fromE_tilde = (beta / Nt[1]) * (St[1] * It[1] + e_vacc * Vt[1] * It[1] + e_asym * St[1] * At[1]) * eval_density(x = timespan[ii+1], theta = c(k_tauL, s_tauL), dist_type = dist_type, density_type = "cdf") + 
      numer_int(xspan = timespan[1:(ii+1)], yval = diff_midpoint_extension(x = (beta / Nt_ext) * (St_ext * It_ext + e_vacc * Vt_ext * It_ext + e_asym * St_ext * At_ext), t_interv = h) * eval_density(x = max(timespan[1:(ii+1)]) - timespan[1:(ii+1)], theta = c(k_tauL, s_tauL), dist_type = dist_type, density_type = "cdf")) +
      Et[1] * eval_density(x = timespan[ii+1], theta = c(k_tauL, s_tauL), dist_type = dist_type, density_type = "pdf")
    fromEh_tilde = ((1+h_var) * beta / Nt[1]) * (St[1] * Iht[1] + e_vacc * Vt[1] * Iht[1] + e_asym * St[1] * Aht[1]) * eval_density(x = timespan[ii+1], theta = c(k_tauLh, s_tauLh), dist_type = dist_type, density_type = "cdf") + 
      numer_int(xspan = timespan[1:(ii+1)], yval = diff_midpoint_extension(x = ((1+h_var) * beta / Nt_ext) * (St_ext * Iht_ext + e_vacc * Vt_ext * Iht_ext + e_asym * St_ext * Aht_ext), t_interv = h) * eval_density(x = max(timespan[1:(ii+1)]) - timespan[1:(ii+1)], theta = c(k_tauLh, s_tauLh), dist_type = dist_type, density_type = "cdf")) +
      Eht[1] * eval_density(x = timespan[ii+1], theta = c(k_tauL, s_tauL), dist_type = dist_type, density_type = "pdf")
    fromItoC_tilde = (1-p_A) * (numer_int(xspan = timespan[1:(ii+1)], yval = (beta / Nt_ext) * (St_ext * It_ext + e_vacc * Vt_ext * It_ext + e_asym * St_ext * At_ext) * conv_tauL_tauI[(ii+1):1]) + Et[1] * conv_tauL_tauI[ii+1]) + 
      It[1] * eval_density(x = timespan[ii+1], theta = c(k_tauI, s_tauI), dist_type = dist_type, density_type = "pdf")
    fromIhtoC_tilde = (1-p_A) * (numer_int(xspan = timespan[1:(ii+1)], yval = ((1+h_var) * beta / Nt_ext) * (St_ext * Iht_ext + e_vacc * Vt_ext * Iht_ext + e_asym * St_ext * Aht_ext) * conv_tauLh_tauI[(ii+1):1]) + Eht[1] * conv_tauLh_tauI[ii+1]) + 
      Iht[1] * eval_density(x = timespan[ii+1], theta = c(k_tauI, s_tauI), dist_type = dist_type, density_type = "pdf")
    fromAtoR_tilde = p_A * (numer_int(xspan = timespan[1:(ii+1)], yval = (beta / Nt_ext) * (St_ext * It_ext + e_vacc * Vt_ext * It_ext + e_asym * St_ext * At_ext) * conv_tauL_tauIA[(ii+1):1]) + Et[1] * conv_tauL_tauIA[ii+1]) + 
      At[1] * eval_density(x = timespan[ii+1], theta = c(k_tauIA, s_tauIA), dist_type = dist_type, density_type = "pdf")
    fromAhtoR_tilde = p_A * (numer_int(xspan = timespan[1:(ii+1)], yval = ((1+h_var) * beta / Nt_ext) * (St_ext * Iht_ext + e_vacc * Vt_ext * Iht_ext + e_asym * St_ext * Aht_ext) * conv_tauLh_tauIA[(ii+1):1]) + Eht[1] * conv_tauLh_tauIA[ii+1]) + 
      Aht[1] * eval_density(x = timespan[ii+1], theta = c(k_tauIA, s_tauIA), dist_type = dist_type, density_type = "pdf")
    
    term10_tilde = (1-p_A) * fromE_tilde
    term11_tilde = p_A * fromE_tilde
    term12_tilde = (1-p_A) * fromEh_tilde
    term13_tilde = p_A * fromEh_tilde
    term14_tilde = fromItoC_tilde
    term15_tilde = fromIhtoC_tilde
    term16_tilde = fromAtoR_tilde
    term17_tilde = fromAhtoR_tilde
    
    
    term1_tilde = max(term1_tilde, 0) # adjusted to avoid negative values 
    term2_tilde = max(term2_tilde, 0) # adjusted to avoid negative values 
    term3_tilde = max(term3_tilde, 0) # adjusted to avoid negative values 
    term4_tilde = max(term4_tilde, 0) # adjusted to avoid negative values 
    term5_tilde = max(term5_tilde, 0) # adjusted to avoid negative values 
    term6_tilde = max(term6_tilde, 0) # adjusted to avoid negative values 
    term7_tilde = max(term7_tilde, 0) # adjusted to avoid negative values 
    term8_tilde = max(term8_tilde, 0) # adjusted to avoid negative values 
    term9_tilde = max(term9_tilde, 0) # adjusted to avoid negative values 
    term10_tilde = max(term10_tilde, 0) # adjusted to avoid negative values 
    term11_tilde = max(term11_tilde, 0) # adjusted to avoid negative values 
    term12_tilde = max(term12_tilde, 0) # adjusted to avoid negative values 
    term13_tilde = max(term13_tilde, 0) # adjusted to avoid negative values 
    term14_tilde = max(term14_tilde, 0) # adjusted to avoid negative values 
    term15_tilde = max(term15_tilde, 0) # adjusted to avoid negative values 
    term16_tilde = max(term16_tilde, 0) # adjusted to avoid negative values 
    term17_tilde = max(term17_tilde, 0) # adjusted to avoid negative values 
    
    increm_mat_tilde = Stoi_mat %*% c(term1_tilde, 
                                      term2_tilde, 
                                      term3_tilde, 
                                      term4_tilde, 
                                      term5_tilde, 
                                      term6_tilde, 
                                      term7_tilde, 
                                      term8_tilde, 
                                      term9_tilde, 
                                      term10_tilde, 
                                      term11_tilde, 
                                      term12_tilde, 
                                      term13_tilde, 
                                      term14_tilde, 
                                      term15_tilde, 
                                      term16_tilde, 
                                      term17_tilde)
    
    S_increm_tilde = increm_mat_tilde[1]
    V_increm_tilde = increm_mat_tilde[2]
    E_increm_tilde = increm_mat_tilde[3]
    Eh_increm_tilde = increm_mat_tilde[4]
    I_increm_tilde = increm_mat_tilde[5]
    Ih_increm_tilde = increm_mat_tilde[6]
    A_increm_tilde = increm_mat_tilde[7]
    Ah_increm_tilde = increm_mat_tilde[8]
    C_increm_tilde = increm_mat_tilde[9]
    R_increm_tilde = increm_mat_tilde[10]
    D_increm_tilde = increm_mat_tilde[11]
    cumul_C_increm_tilde = + term14_tilde + term15_tilde
    cumul_asymp_increm_tilde = + term16_tilde + term17_tilde
    
    # S_increm_tilde = - term1_tilde - term2_tilde - term3_tilde - term4_tilde - term5_tilde
    # V_increm_tilde = + term5_tilde - term6_tilde - term7_tilde
    # E_increm_tilde = + term1_tilde + term2_tilde + term6_tilde - term10_tilde - term11_tilde
    # Eh_increm_tilde = + term3_tilde + term4_tilde + term7_tilde - term12_tilde - term13_tilde
    # I_increm_tilde = + term10_tilde - term14_tilde
    # Ih_increm_tilde = + term12_tilde - term15_tilde
    # A_increm_tilde = + term11_tilde - term16_tilde
    # Ah_increm_tilde = + term13_tilde - term17_tilde
    # C_increm_tilde = - term8_tilde - term9_tilde + term14_tilde + term15_tilde
    # R_increm_tilde = + term8_tilde + term16_tilde + term17_tilde
    # D_increm_tilde = + term9_tilde
    
    St[ii+1] = St[ii] + h/2 * (S_increm + S_increm_tilde)
    Vt[ii+1] = Vt[ii] + h/2 * (V_increm + V_increm_tilde)
    Et[ii+1] = Et[ii] + h/2 * (E_increm + E_increm_tilde)
    Eht[ii+1] = Eht[ii] + h/2 * (Eh_increm + Eh_increm_tilde)
    It[ii+1] = It[ii] + h/2 * (I_increm + I_increm_tilde)
    Iht[ii+1] = Iht[ii] + h/2 * (Ih_increm + Ih_increm_tilde)
    At[ii+1] = At[ii] + h/2 * (A_increm + A_increm_tilde)
    Aht[ii+1] = Aht[ii] + h/2 * (Ah_increm + Ah_increm_tilde)
    Ct[ii+1] = Ct[ii] + h/2 * (C_increm + C_increm_tilde)
    Rt[ii+1] = Rt[ii] + h/2 * (R_increm + R_increm_tilde)
    Dt[ii+1] = Dt[ii] + h/2 * (D_increm + D_increm_tilde)
    Nt[ii+1] = sum(St[ii+1], Vt[ii+1], Et[ii+1], Eht[ii+1], It[ii+1], Iht[ii+1], At[ii+1], Aht[ii+1], Ct[ii+1], Rt[ii+1], Dt[ii+1]) - Dt[ii+1]
    cumul_Ct[ii+1] = cumul_Ct[ii] + h/2 * (cumul_C_increm + cumul_C_increm_tilde)
    cumul_asympt[ii+1] = cumul_asympt[ii] + h/2 * (cumul_asymp_increm + cumul_asymp_increm_tilde)
  }
  return(list("St" = St, "Vt" = Vt, "Et" = Et, "Eht" = Eht, "It" = It, "Iht" = Iht, 
              "At" = At, "Aht" = Aht, "Ct" = Ct, "Rt" = Rt, "Dt" = Dt, "Nt" = Nt, "cumul_Ct" = cumul_Ct, "cumul_asympt" = cumul_asympt))
}


mean_trajectory_SEIR_hist <- function(timespan, theta, y_init){
  beta.param = theta[1]
  alpha1 = theta[2]
  beta1 = theta[3]
  alpha2 = theta[4]
  beta2 = theta[5]
  
  Ntot = sum(y_init)
  
  h = timespan[2] - timespan[1]
  xspan = timespan
  St = rep(NA, length(xspan))
  Et = rep(NA, length(xspan))
  It = rep(NA, length(xspan))
  Rt = rep(NA, length(xspan))
  conv_val = rep(NA, length(xspan))
  conv_val2 = rep(NA, length(xspan))
  
  # hs = rep(NA, length(xspan))
  # hs_tilde = rep(NA, length(xspan))
  St[1] = y_init[1]
  Et[1] = y_init[2]
  It[1] = y_init[3]
  Rt[1] = y_init[4]
  
  for(ii in 1:length(xspan)){
    conv_val[ii] = conv_two_gamma(tmax = xspan[ii], alpha1, beta1, alpha2, beta2)
  }
  
  for(ii in 1:length(xspan)){
    conv_val2[ii] = conv_gam_cdf_pdf(tmax = xspan[ii], alpha1, beta1, alpha2, beta2) / (alpha1 * beta1) 
  }
  
  for(ii in 1:(length(xspan)-1)){
    
    term1 = beta.param * St[ii] * It[ii] / Ntot
    S_increm = - term1
    term2 = beta.param * St[1] * It[1] / Ntot * pgamma(xspan[ii], shape = alpha1, scale = beta1) + 
      numer_int(xspan[1:ii], beta.param * diff_midpoint_extension(St[1:ii] * It[1:ii], xspan[2]-xspan[1]) / Ntot * pgamma(max(xspan[1:ii]) - xspan[1:ii], shape = alpha1, scale = beta1)) + 
      Et[1] * pgamma(xspan[ii], shape = alpha1, scale = beta1, lower.tail = FALSE) / (alpha1 * beta1) ### edit this part. 
    E_increm = term1 - term2
    term3 = numer_int(xspan[1:ii], beta.param * St[1:ii] * It[1:ii] / Ntot * conv_val[ii:1]) + 
      It[1] * pgamma(xspan[ii], shape = alpha2, scale = beta2, lower.tail = FALSE) / (alpha2 * beta2) + Et[1] * conv_val2[ii]  ### edit this part. 
    term3 = max(term3, 0) # adjusted to avoid negative values 
    I_increm = term2 - term3
    R_increm = term3
    
    St_tilde = St[ii] + h*S_increm
    Et_tilde = Et[ii] + h*E_increm
    It_tilde = It[ii] + h*I_increm
    Rt_tilde = Rt[ii] + h*R_increm
    
    term1_tilde = beta.param * St_tilde * It_tilde / Ntot
    S_increm_tilde = - term1_tilde
    term2_tilde = beta.param * St[1] * It[1] / Ntot * pgamma(xspan[ii+1], shape = alpha1, scale = beta1) + 
      numer_int(xspan[1:(ii+1)], beta.param * diff_midpoint_extension(c(St[1:ii], St_tilde) * c(It[1:ii], It_tilde), xspan[2] - xspan[1]) / Ntot * pgamma(max(xspan[1:(ii+1)]) - xspan[1:(ii+1)], shape = alpha1, scale = beta1)) + 
      Et[1] * pgamma(xspan[ii+1], shape = alpha1, scale = beta1, lower.tail = FALSE) / (alpha1 * beta1) 
    E_increm_tilde = term1_tilde - term2_tilde
    term3_tilde = numer_int(xspan[1:(ii+1)], beta.param * c(St[1:ii], St_tilde) * c(It[1:ii], It_tilde) / Ntot * conv_val[(ii+1):1])  + 
      It[1] * pgamma(xspan[ii+1], shape = alpha2, scale = beta2, lower.tail = FALSE) / (alpha2 * beta2) + Et[1] * conv_val2[ii+1] 
    I_increm_tilde = term2_tilde - term3_tilde
    term3_tilde = max(term3_tilde, 0) # adjusted to avoid negative values 
    R_increm_tilde = term3_tilde
    
    St[ii+1] = St[ii] + h/2 * (S_increm + S_increm_tilde)
    Et[ii+1] = Et[ii] + h/2 * (E_increm + E_increm_tilde)
    It[ii+1] = It[ii] + h/2 * (I_increm + I_increm_tilde)
    Rt[ii+1] = Rt[ii] + h/2 * (R_increm + R_increm_tilde)
    
  }
  return(list("St" = St, "Et" = Et, "It" = It, "Rt" = Rt))
}

