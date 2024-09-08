setwd("/usrdirectory")
source("IONISE_source.R") # load the library of source codes needed to perform the package.

### ===================== CODE for reproducing FIGURE 3 =================== ###
set.seed(12) # set the random seed

args1 = 1 # From 1 to 12, each number represent one specific setting in Figure 3
flg1 = 0
set_comb = matrix(NA, nrow = 12, ncol = 5)

for(estim_infect in c(1,0)){ 
    for(dt_data in c(2)){
        for(dt_infer in c(1,2)){
            for(R0_set in c(1,2,4)){
                for(maxT in c(43)){
                    flg1 = flg1 +1
                    set_comb[flg1, ] = c(estim_infect, dt_data, dt_infer, R0_set, maxT)
                    
                }
            }
        } 
    }
} 


estim_infect = set_comb[args1, 1]
dt_data = set_comb[args1, 2]
dt_infer = set_comb[args1, 3]
R0_set = set_comb[args1, 4]
maxT = set_comb[args1, 5]

dist_type_list = c("exp", "gamma", "invgamma", "lognormal", "weibull")
dist_type_data = dist_type_list[dt_data]
dist_type_infer = dist_type_list[dt_infer]

# estimate beta, tau1, or tau2. if estim_TF = c(1,1,0) then the beta and tau1 are estimated while the tau2 is fixed as its true value.

# data_R <- read.csv("input_data_example.csv")
# data_R <- read.csv("input_data.csv")

####  ========================= Initialize variables ===========================
####  =========== Users need to specify variables in this section ==============

initial_phase = TRUE # TRUE if the data comes from the initial phase of the disease spread

Ntot = 1e+07
E_init = 0
I_init = 1
R_init = 0
S_init = Ntot - (E_init + I_init + R_init)

latent_mean = 5.48
latent_var = 7.40

# If estim_infect is TRUE then infect_mean and infect_var are used as the initial values for MCMC sampling.
# estim_infect = TRUE # TRUE if the user wants to estimate the infectious period distribution 
infect_mean = 6
infect_var = 1.2

latent_par1 = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type_data)$par1
latent_par2 = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type_data)$par2
infect_par1 = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type_data)$par1
infect_par2 = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type_data)$par2

theta_init = c(R0_set / infect_mean, latent_par1, latent_par2, infect_par1, infect_par2)

timespan0 = 0:maxT
h = timespan0[2] - timespan0[1]
dd = get_divisor_dist(theta = c(latent_par1,latent_par2,infect_par1,infect_par2), dist_type = dist_type_data)
ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
if (dist_type_data == "exp"){
    init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
    theta_ode <- c(beta = theta_init[1], s1 = theta_init[3]^(-1), s2 = theta_init[5]^(-1), N = sum(init))
    sol = ode(y = init, times=ext_tspan, func=SEIR11, parms = theta_ode)
    S_true = sol[subseq0,2]
    E_true = sol[subseq0,3]
    I_true = sol[subseq0,4]
    R_true = sol[subseq0,5]
}else{
    if(initial_phase){
        sim_cand = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = theta_init, y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type_data)
    }else{
        warning("Only initial phase is valid!")
    }
    S_true = sim_cand$St[subseq0]
    E_true = sim_cand$Et[subseq0]
    I_true = sim_cand$It[subseq0]
    R_true = sim_cand$Rt[subseq0]
}


data_R = cbind(rep(NA, length(diff(R_true))), diff(R_true))

# Specify the mean and variance for prior distributions of beta and latent period distribution parameters
prior_mean_beta = 1
prior_mean_infect_shape = 1
prior_mean_infect_scale = 1
prior_var_beta = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_shape = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_scale = 10^6 # large prior variance indicates the non-informative prior

# Specify parameters for MCMC sampling
# the expected run time highly depends on the number of iterations for the MCMC algorithm, the number of data points, and other settings.
burn = 10000 #The length of the burn-in period, which determines the number of initial MCMC samples discarded to ensure that the remaining MCMC samples are independent from their initial value. 
thinning = 1/100 # The thinning rate of the MCMC iteration. 
jump = round(1/thinning) # The reciprocal of the thinning rate.
effnum = 1000 #The number of effective MCMC iterations. Its default value is 1000. 
nrepeat = burn + effnum * jump # The total number of MCMC iterations is burn + effnum * jump.
tun_beta = 0.1 # tuning parameter; proposal variance for the MH algorithm sampling the transmission rate, beta.
selrow = seq(from = burn + jump, to = nrepeat, by = jump)

####  ===========================  Perform MCMC  ===============================
mcmc.result <- MCMC_function_SEIR_pois(data_R, S_init, E_init, I_init, R_init, 
                                       initial_phase, estim_infect,
                                       latent_mean, latent_var, infect_mean, infect_var, 
                                       prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                       nrepeat, tun_beta, dist_type = dist_type_infer, filename = "test.RData")

mcmc.result.proc = mcmc.result[selrow,]



### ===================== CODE for reproducing FIGURES 4 & 5 =================== ###
set.seed(12) # set the random seed

args1 = 1 # From 1 to 12, each number represent one specific setting in Figures 4 & 5.
flg1 = 0
set_comb = matrix(NA, nrow = 12, ncol = 5)

for(estim_infect in c(1)){ 
    for(dt_data in c(2)){
        for(dt_infer in c(1,2)){
            for(period in c(2,3)){
                for(k_infer in c(0.5, 1, 1.5)){
                    flg1 = flg1 +1
                    set_comb[flg1, ] = c(estim_infect, dt_data, dt_infer, period, k_infer)
                }
            }
        } 
    }
} 

estim_infect = set_comb[args1, 1]
dt_data = set_comb[args1, 2]
dt_infer = set_comb[args1, 3]
period = set_comb[args1, 4]
k_infer = set_comb[args1, 5]

dist_type_list = c("exp", "gamma", "invgamma", "lognormal", "weibull")
dist_type_data = dist_type_list[dt_data]
dist_type_infer = dist_type_list[dt_infer]

####  ========================= Initialize variables ===========================
####  =========== Users need to specify variables in this section ==============

initial_phase = FALSE # TRUE if the data comes from the initial phase of the disease spread

Ntot = 1e+07
if (period == 2){
    maxT = 19
    E_init_data = round(10*5.5)
    E_init_infer = round(k_infer * E_init_data)
    I_init = 25
    R_init = 1607
    infect_var = 5*(6/5)^2
    beta.param = 1.22
}else if(period == 3){
    maxT = 25
    E_init_data = round(41*5.5)
    E_init_infer = round(k_infer * E_init_data)
    I_init = 173
    R_init = 6082
    infect_var = 1.25*(6/1.25)^2
    beta.param = 0.37
}
S_init = Ntot - (E_init_data + I_init + R_init)

latent_mean = 5.48
latent_var = 7.40

infect_mean = 6

latent_par1 = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type_data)$par1
latent_par2 = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type_data)$par2
infect_par1 = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type_data)$par1
infect_par2 = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type_data)$par2

theta_init = c(beta.param, latent_par1, latent_par2, infect_par1, infect_par2)

timespan0 = 0:maxT
h = timespan0[2] - timespan0[1]
dd = get_divisor_dist(theta = c(latent_par1,latent_par2,infect_par1,infect_par2), dist_type = dist_type_data)
ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)

if (dist_type_data == "exp"){
    init = c(S=S_init, E1 = E_init_data, I1 = I_init, R = R_init)
    theta_ode <- c(beta = theta_init[1], s1 = theta_init[3]^(-1), s2 = theta_init[5]^(-1), N = sum(init))
    sol = ode(y = init, times=ext_tspan, func=SEIR11, parms = theta_ode)
    S_true = sol[subseq0,2]
    E_true = sol[subseq0,3]
    I_true = sol[subseq0,4]
    R_true = sol[subseq0,5]
}else{
    if(initial_phase){
        sim_cand = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = theta_init, y_init = c(S_init, E_init_data, I_init, R_init), dist_type = dist_type_data)
    }else{
        sim_cand = mean_trajectory_SEIR_hist(timespan = ext_tspan, theta = theta_init, y_init = c(S_init, E_init_data, I_init, R_init))
    }
    S_true = sim_cand$St[subseq0]
    E_true = sim_cand$Et[subseq0]
    I_true = sim_cand$It[subseq0]
    R_true = sim_cand$Rt[subseq0]
}

data_R = cbind(rep(NA, length(diff(R_true))), diff(R_true))


# Specify the mean and variance for prior distributions of beta and latent period distribution parameters
prior_mean_beta = 1
prior_mean_infect_shape = 1
prior_mean_infect_scale = 1
prior_var_beta = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_shape = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_scale = 10^6 # large prior variance indicates the non-informative prior

# Specify parameters for MCMC sampling
# the expected run time highly depends on the number of iterations for the MCMC algorithm, the number of data points, and other settings.
burn = 10000 #The length of the burn-in period, which determines the number of initial MCMC samples discarded to ensure that the remaining MCMC samples are independent from their initial value. 
thinning = 1/100 # The thinning rate of the MCMC iteration. 
jump = round(1/thinning) # The reciprocal of the thinning rate.
effnum = 1000 #The number of effective MCMC iterations. Its default value is 1000. 
nrepeat = burn + effnum * jump # The total number of MCMC iterations is burn + effnum * jump.
tun_beta = 0.1 # tuning parameter; proposal variance for the MH algorithm sampling the transmission rate, beta.
selrow = seq(from = burn + jump, to = nrepeat, by = jump)

####  ===========================  Perform MCMC  ===============================
mcmc.result <- MCMC_function_SEIR_pois(data_R, S_init, E_init, I_init, R_init, 
                                       initial_phase, estim_infect,
                                       latent_mean, latent_var, infect_mean, infect_var, 
                                       prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                       nrepeat, tun_beta, dist_type = dist_type_infer, filename = "test.RData")

mcmc.result.proc = mcmc.result[selrow,]



### ===================== CODE for reproducing FIGURE 6 =================== ###
set.seed(12) # set the random seed

args1 = 1 # From 1 to 12, each number represent one specific setting in Figure 6.
flg1 = 0
set_comb = matrix(NA, nrow = 12, ncol = 3)
for(dt_data in c(1,2,4,5)){
    for(dt_infer in c(dt_data)){
        for(R0_set in c(1,2,4)){
            flg1 = flg1 +1
            set_comb[flg1, ] = c(dt_data, dt_infer, R0_set)    
        }
    } 
}

estim_infect = 1
dt_data = set_comb[args1, 1]
dt_infer = set_comb[args1, 2]
R0_set = set_comb[args1, 3]
maxT = 43

dist_type_list = c("exp", "gamma", "invgamma", "lognormal", "weibull")
dist_type_data = dist_type_list[dt_data]
dist_type_infer = dist_type_list[dt_infer]

# estimate beta, tau1, or tau2. if estim_TF = c(1,1,0) then the beta and tau1 are estimated while the tau2 is fixed as its true value.

# data_R <- read.csv("input_data_example.csv")
# data_R <- read.csv("input_data.csv")

####  ========================= Initialize variables ===========================
####  =========== Users need to specify variables in this section ==============

initial_phase = TRUE # TRUE if the data comes from the initial phase of the disease spread

S_init = 1e+07
E_init = 0
I_init = 1
R_init = 0

latent_mean = 5.48
latent_var = 7.40

# If estim_infect is TRUE then infect_mean and infect_var are used as the initial values for MCMC sampling.
# estim_infect = TRUE # TRUE if the user wants to estimate the infectious period distribution 
infect_mean = 6
infect_var = 10

latent_par1 = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type_data)$par1
latent_par2 = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type_data)$par2
infect_par1 = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type_data)$par1
infect_par2 = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type_data)$par2

theta_init = c(R0_set / infect_mean, latent_par1, latent_par2, infect_par1, infect_par2)

timespan0 = 0:maxT
h = timespan0[2] - timespan0[1]
dd = get_divisor_dist(theta = c(latent_par1,latent_par2,infect_par1,infect_par2), dist_type = dist_type_data)
ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
if (dist_type_data == "exp"){
    init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
    theta_ode <- c(beta = theta_init[1], s1 = theta_init[3]^(-1), s2 = theta_init[5]^(-1), N = sum(init))
    sol = ode(y = init, times=ext_tspan, func=SEIR11, parms = theta_ode)
    S_true = sol[subseq0,2]
    E_true = sol[subseq0,3]
    I_true = sol[subseq0,4]
    R_true = sol[subseq0,5]
}else{
    if(initial_phase){
        sim_cand = mean_trajectory_SEIR_dist_input(timespan = ext_tspan, theta = theta_init, y_init = c(S_init, E_init, I_init, R_init), dist_type = dist_type_data)
    }else{
        warning("Only initial phase is valid!")
    }
    S_true = sim_cand$St[subseq0]
    E_true = sim_cand$Et[subseq0]
    I_true = sim_cand$It[subseq0]
    R_true = sim_cand$Rt[subseq0]
}


data_R = cbind(rep(NA, length(diff(R_true))), diff(R_true))

# Specify the mean and variance for prior distributions of beta and latent period distribution parameters
prior_mean_beta = 1
prior_mean_infect_shape = 1
prior_mean_infect_scale = 1
prior_var_beta = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_shape = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_scale = 10^6 # large prior variance indicates the non-informative prior

# Specify parameters for MCMC sampling
# the expected run time highly depends on the number of iterations for the MCMC algorithm, the number of data points, and other settings.
burn = 10000 #The length of the burn-in period, which determines the number of initial MCMC samples discarded to ensure that the remaining MCMC samples are independent from their initial value. 
thinning = 1/100 # The thinning rate of the MCMC iteration. 
jump = round(1/thinning) # The reciprocal of the thinning rate.
effnum = 1000 #The number of effective MCMC iterations. Its default value is 1000. 
nrepeat = burn + effnum * jump # The total number of MCMC iterations is burn + effnum * jump.
tun_beta = 0.1 # tuning parameter; proposal variance for the MH algorithm sampling the transmission rate, beta.
selrow = seq(from = burn + jump, to = nrepeat, by = jump)

####  ===========================  Perform MCMC  ===============================
mcmc.result <- MCMC_function_SEIR_pois(data_R, S_init, E_init, I_init, R_init, 
                                       initial_phase, estim_infect,
                                       latent_mean, latent_var, infect_mean, infect_var,
                                       prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                       nrepeat, tun_beta, dist_type = dist_type_infer, filename = "test.RData")
mcmc.result.proc = mcmc.result[selrow,]


### ===================== CODE for reproducing FIGURE 7 =================== ###

set.seed(12) # set the random seed

args1 = 1 # From 1 to 2, each number represent one specific setting in Figure 7.
flg1 = 0
set_comb = matrix(NA, nrow = 2, ncol = 2)

for(dt_data in c(2)){
    for(dt_infer in c(1,2)){
        flg1 = flg1 +1
        set_comb[flg1, ] = c(dt_data, dt_infer)            
    }
}

estim_infect = 1
R0_set = 2
maxT = 43
dt_data = set_comb[args1, 1]
dt_infer = set_comb[args1, 2]

R0_data = R0_set # vary this value in c(1,2,4,8)
mean_tauL = 5.52
var_tauL = 5.24
mean_tauLh = 4.43
var_tauLh = 4.80
mean_tauI = 6
var_tauI = 10
mean_tauIA = 7.25
var_tauIA = 24.64

dt_list = c("exp", "gamma")
dist_type_data = dt_list[dt_data]
dist_type_infer = dt_list[dt_infer]

beta = R0_data/mean_tauI
h_var = 0.6
e_asym = 0.5
e_vacc = 0.1
vacc = 0.005
rec = 0.1
death = 0.01
p_A = 0.3

timespan0 = 0:maxT
h = timespan0[2] - timespan0[1]
tmp = mean_var_to_param(mean_tauL, var_tauL, dist_type = dist_type_data)
k_tauL_data = tmp$par1
s_tauL_data = tmp$par2
tmp = mean_var_to_param(mean_tauLh, var_tauLh, dist_type = dist_type_data)
k_tauLh_data = tmp$par1
s_tauLh_data = tmp$par2
tmp = mean_var_to_param(mean_tauI, var_tauI, dist_type = dist_type_data)
k_tauI_data = tmp$par1
s_tauI_data = tmp$par2
tmp = mean_var_to_param(mean_tauIA, var_tauIA, dist_type = dist_type_data)
k_tauIA_data = tmp$par1
s_tauIA_data = tmp$par2

param_set_data = list(beta = beta, h_var = h_var, e_asym = e_asym, e_vacc = e_asym, 
                      vacc = vacc, rec = rec, death = death, p_A = p_A,
                      k_tauL = k_tauL_data, s_tauL = s_tauL_data, k_tauLh = k_tauLh_data, s_tauLh = s_tauLh_data,
                      k_tauI = k_tauI_data, s_tauI = s_tauI_data, k_tauIA = k_tauIA_data, s_tauIA = s_tauIA_data)

tmp = mean_var_to_param(mean_tauL, var_tauL, dist_type = dist_type_infer)
k_tauL_infer = tmp$par1
s_tauL_infer = tmp$par2
tmp = mean_var_to_param(mean_tauLh, var_tauLh, dist_type = dist_type_infer)
k_tauLh_infer = tmp$par1
s_tauLh_infer = tmp$par2
tmp = mean_var_to_param(mean_tauI, var_tauI, dist_type = dist_type_infer)
k_tauI_infer = tmp$par1
s_tauI_infer = tmp$par2
tmp = mean_var_to_param(mean_tauIA, var_tauIA, dist_type = dist_type_infer)
k_tauIA_infer = tmp$par1
s_tauIA_infer = tmp$par2

param_set_infer = list(beta = beta, h_var = h_var, e_asym = e_asym, e_vacc = e_vacc, 
                       vacc = vacc, rec = rec, death = death, p_A = p_A,
                       k_tauL = k_tauL_infer, s_tauL = s_tauL_infer, k_tauLh = k_tauLh_infer, s_tauLh = s_tauLh_infer,
                       k_tauI = k_tauI_infer, s_tauI = s_tauI_infer, k_tauIA = k_tauIA_infer, s_tauIA = s_tauIA_infer)

y_init0 = list(S_init = 10^8, V_init = 0, E_init = 0, Eh_init = 0, I_init = 1, Ih_init = 1, 
               A_init = 0, Ah_init = 0, C_init = 0, R_init = 0, D_init = 0)

dd_lower = get_divisor_dist(c(param_set_data$k_tauLh, param_set_data$s_tauLh, param_set_data$k_tauIA, param_set_data$s_tauIA), dist_type = dist_type_data)
dd = max(dd_lower, get_divisor_dist(c(param_set_data$k_tauL, param_set_data$s_tauL, param_set_data$k_tauI, param_set_data$s_tauI), dist_type = dist_type_data))

ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
result0 = mean_trajectory_LARGE_dist_input(timespan = ext_tspan, param_set = param_set_data, y_init = y_init0, dist_type = dist_type_data)

cumul_true = result0$cumul_Ct[subseq0]
data_daily = cbind(rep(NA, length(diff(cumul_true))), diff(cumul_true))

# Specify the mean and variance for prior distributions of beta and latent period distribution parameters
prior_mean_beta = 1
prior_mean_infect_shape = 1
prior_mean_infect_scale = 1
prior_var_beta = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_shape = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_scale = 10^6 # large prior variance indicates the non-informative prior

# Specify parameters for MCMC sampling
# the expected run time highly depends on the number of iterations for the MCMC algorithm, the number of data points, and other settings.
burn = 10000 #The length of the burn-in period, which determines the number of initial MCMC samples discarded to ensure that the remaining MCMC samples are independent from their initial value. 
thinning = 1/100 # The thinning rate of the MCMC iteration. 
jump = round(1/thinning) # The reciprocal of the thinning rate.
effnum = 1000 #The number of effective MCMC iterations. Its default value is 1000. 
nrepeat = burn + effnum * jump # The total number of MCMC iterations is burn + effnum * jump.
tun_beta = 0.1 # tuning parameter; proposal variance for the MH algorithm sampling the transmission rate, beta.
selrow = seq(from = burn + jump, to = nrepeat, by = jump)

####  ===========================  Perform MCMC  ===============================

mcmc.result <- MCMC_function_extended_pois(data_daily = data_daily, y_init = y_init0, param_set = param_set_infer, estim_infect,
                                           prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                                           nrepeat, tun_beta, dist_type = dist_type_infer, filename = "test.RData")

mcmc.result.proc = mcmc.result[selrow,]


