setwd("/Users/hyukpyohong/Library/CloudStorage/Dropbox/Delayed_Epidemic_model/delay_epidemic_code/IONISE_v3") # set the working directory contains the code file named "IONISE_functions.R"
source("IONISE_functions.R") # load the library of source codes needed to perform the package.

set.seed(1234) # set the random seed

args0=(commandArgs(TRUE))
args1 = as.numeric(args0[1])
flg1 = 0
set_comb = matrix(NA, nrow = 192, ncol = 5)

for(estim_infect in c(0, 1)){ 
    for(lik_type in c(1,2)){
        for(dt_data in c(1,2,4,5)){
            for(dt_infer in c(1,2,4,5)){
                for(R0_set in c(1,2,4)){
                    flg1 = flg1 +1
                    set_comb[flg1, ] = c(estim_infect, lik_type, dt_data, dt_infer, R0_set)
                }
            } 
        }
    } 
}

estim_infect = set_comb[args1, 1]
lik_type = set_comb[args1, 2]
dt_data = set_comb[args1, 3]
dt_infer = set_comb[args1, 4]
R0_set = set_comb[args1, 5]

dist_type_list = c("exp", "gamma", "invgamma", "lognormal", "weibull")
dist_type_data = dist_type_list[dt_data]
dist_type_infer = dist_type_list[dt_infer]

# estimate beta, tau1, or tau2. if estim_TF = c(1,1,0) then the beta and tau1 are estimated while the tau2 is fixed as its true value.

# data_R <- read.csv("input_data_example.csv")
# data_R <- read.csv("input_data.csv")

####  ========================= Initialize variables ===========================
####  =========== Users need to specify variables in this section ==============

initial_phase = TRUE # TRUE if the data comes from the initial phase of the disease spread

S_init = 100000
E_init = 0
I_init = 1
R_init = 0

# The default values(COVID-19) below are from Lauer et al, 2020, Annals of Internal Medicine.
latent_mean = 5.51
latent_var = 5.22

# If estim_infect is TRUE then infect_mean and infect_var are used as the initial values for MCMC sampling.
# estim_infect = TRUE # TRUE if the user wants to estimate the infectious period distribution 
infect_mean = 6
infect_var = 10

latent_par1 = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type_data)$par1
latent_par2 = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type_data)$par2
infect_par1 = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type_data)$par1
infect_par2 = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type_data)$par2

theta_init = c(R0_set / infect_mean, latent_par1, latent_par2, infect_par1, infect_par2)

maxT = 43

timespan0 = 0:maxT
h = timespan0[2] - timespan0[1]
dd = get_divisor_dist(theta = c(latent_par1,latent_par2,infect_par1,infect_par2), dist_type = dist_type_data)
ext_tspan = seq(from = 0, to = max(timespan0), by = h/dd)
subseq0 = seq(from = 1, to = dd * max(timespan0)/h + 1, by = dd)
if (dist_type_data == "exp"){
    init = c(S=S_init, E1 = E_init, I1 = I_init, R = R_init)
    theta_ode <- c(beta = theta_init[1], s1 = theta_init[3]^(-1), s2 = theta_init[5]^(-1), N = sum(init))
    sol = ode(y = init, times=ext_tspan, func=SEIR11, parms = theta_ode)
    S_true = sol[,2]
    E_true = sol[,3]
    I_true = sol[,4]
    R_true = sol[,5]
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

data_R = cbind(rep(NA, length(R_true)),R_true)

# Specify the mean and variance for prior distributions of beta and latent period distribution parameters
prior_mean_beta = 1
prior_mean_latent_shape = 1
prior_mean_latent_scale = 1
prior_var_beta = 1000 # large prior variance indicates the non-informative prior
prior_var_latent_shape = 1000 # large prior variance indicates the non-informative prior
prior_var_latent_scale = 1000 # large prior variance indicates the non-informative prior

# Specify parameters for MCMC sampling
nrepeat = 100
tun_beta = 0.1 # tuning parameter; proposal variance for the MH algorithm sampling the transmission rate, beta. 

filename = paste("Revision_v1_maxt",maxT,"est_infect",estim_infect,"lik_type", lik_type, "dt_data",dt_data, "dt_infer",dt_infer, "R0_",R0_set,".RData", sep = "")

####  ===========================  Perform MCMC  ===============================
if(lik_type == 1){
    mcmc.result <- MCMC_function(data_R, S_init, E_init, I_init, R_init, 
                                      initial_phase, estim_infect,
                                      latent_mean, latent_var, infect_mean, infect_var,
                                      prior_mean_beta, prior_mean_latent_shape, prior_mean_latent_scale, prior_var_beta, prior_var_latent_shape, prior_var_latent_scale, 
                                      nrepeat, tun_beta, dist_type = dist_type_infer, filename = filename)
}else{
    mcmc.result <- MCMC_function_pois(data_R, S_init, E_init, I_init, R_init, 
                                      initial_phase, estim_infect,
                                      latent_mean, latent_var, infect_mean, infect_var,
                                      prior_mean_beta, prior_mean_latent_shape, prior_mean_latent_scale, prior_var_beta, prior_var_latent_shape, prior_var_latent_scale, 
                                      nrepeat, tun_beta, dist_type = dist_type_infer, filename = filename)
}


# write.table(x = cbind(mcmc.result[,c(1,4,5)], mcmc.result[,4]*mcmc.result[,5], mcmc.result[,1] * mcmc.result[,4] * mcmc.result[,5]),
#             file="post_samples.csv",sep=",", 
#             col.names=c("beta", "infect_shape", "infect_scale", "mean_infect", "reproduction"), row.names=F)

save.image(file = filename)