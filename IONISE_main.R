setwd("/Users/hyukpyohong/Dropbox/Delay_SEIR_paper/Release code") # set the working directory contains the code file named "IONISE_source.R"
# load the module of source codes needed to perform the package.
source("IONISE_source.R")

####  ======== PART 1: Common setting ==========
set.seed(1234) # set the random seed
# load the input data for estimation
data_R <- read.csv("input_data_example.csv")
# data_R <- read.csv("input_data.csv")

initial_phase = TRUE # TRUE if the data comes from the initial phase of the disease spread
estim_infect = TRUE # TRUE if the user wants to estimate the infectious period distribution 

# Specify parameters for MCMC sampling
# the expected run time highly depends on the number of iterations for the MCMC algorithm, the number of data points, and other settings.
nrepeat = 1000 # for the example data, the expected run time is approximately 1 hours per 10,000 iterations.
tun_beta = 0.1 # tuning parameter; proposal variance for the MH algorithm sampling the transmission rate, beta. 

# specify the distribution type for the delay distributions in the model (i.e., latent and infectious distributions)
# dist_type is one of c("exp", "gamma", "invgamma", "lognormal", "weibull")
dist_type = "gamma"

# specify the type of the likelihood function for the inference.
# lik_type is one of c("poisson", "gaussian")
lik_type = "gaussian"

# Specify the mean and variance for prior distributions of beta and infectious period distribution parameters
prior_mean_beta = 1
prior_mean_infect_shape = 1
prior_mean_infect_scale = 1
prior_var_beta = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_shape = 10^6 # large prior variance indicates the non-informative prior
prior_var_infect_scale = 10^6 # large prior variance indicates the non-informative prior

model_type = "SEIR" # this code is for the SEIR model in Fig.1.

####  ======== PART 2-1: If you choose to run the SEIR model, adjust PART 2-1 ==========
if(model_type == "SEIR"){
    # specify the initial conditions
    # The SEIR model requires to specify the initial conditions of the four variables.
    y_init = list(S_init = 10^8, E_init = 0, I_init = 1, R_init = 0) 
    
    # The default values(COVID-19) below are from Lauer et al, 2020, Annals of Internal Medicine.
    latent_mean = 5.51
    latent_var = 5.22
    
    # If estim_infect is TRUE then infect_mean and infect_var are used as the initial values for MCMC sampling.
    infect_mean = 6
    infect_var = 10
}

####  ======== PART 2-2: If you choose to run the extended SEIR model, adjust PART 2-2 ==========
if(model_type == "Extended_SEIR"){
    
    y_init = list(S_init = 10^8, V_init = 0, E_init = 0, Eh_init = 0, I_init = 1, Ih_init = 1, 
                  A_init = 0, Ah_init = 0, C_init = 0, R_init = 0, D_init = 0)
    
    # The mean and variances for the delay distributions in the model
    # tauL: latent period 
    # tauLh: latent period of the highly contagious variant
    # tauI: infectious period
    # tauIA: infectious period for asymptomatic individuals
    mean_tauL = 5.52
    var_tauL = 5.24
    mean_tauLh = 4.43
    var_tauLh = 4.80
    mean_tauI = 6
    var_tauI = 10
    mean_tauIA = 7.25
    var_tauIA = 24.64
    
    # specify the parameter values
    beta = 0.17 # the transmission rate beta is estimated. Thus, this value is used as the initial value for MCMC sampling.
    h_var = 0.6
    e_asym = 0.5
    e_vacc = 0.1
    vacc = 0.005
    rec = 0.1
    death = 0.01
    p_A = 0.3
}

####  ======== PART 3: Perform MCMC. Users do not need to adjust this part ==============

if(model_type == "SEIR"){
    tmp = mean_var_to_param(latent_mean, latent_var, dist_type = dist_type)
    k_tauL = tmp$par1
    s_tauL= tmp$par2
    tmp = mean_var_to_param(infect_mean, infect_var, dist_type = dist_type)
    k_tauI = tmp$par1
    s_tauI = tmp$par2
    
    param_set = list(k_tauL = k_tauL, s_tauL = s_tauL, k_tauI = k_tauI, s_tauI = s_tauI)
}
if(model_type == "Extended_SEIR"){
    tmp = mean_var_to_param(mean_tauL, var_tauL, dist_type = dist_type)
    k_tauL = tmp$par1
    s_tauL = tmp$par2
    tmp = mean_var_to_param(mean_tauLh, var_tauLh, dist_type = dist_type)
    k_tauLh = tmp$par1
    s_tauLh = tmp$par2
    tmp = mean_var_to_param(mean_tauI, var_tauI, dist_type = dist_type)
    k_tauI = tmp$par1
    s_tauI = tmp$par2
    tmp = mean_var_to_param(mean_tauIA, var_tauIA, dist_type = dist_type)
    k_tauIA = tmp$par1
    s_tauIA = tmp$par2
    
    param_set = list(beta = beta, h_var = h_var, e_asym = e_asym, e_vacc = e_vacc, 
                     vacc = vacc, rec = rec, death = death, p_A = p_A,
                     k_tauL = k_tauL, s_tauL = s_tauL, k_tauLh = k_tauLh, s_tauLh = s_tauLh,
                     k_tauI = k_tauI, s_tauI = s_tauI, k_tauIA = k_tauIA, s_tauIA = s_tauIA)
}

mcmc.result = MCMC_function(data_R, y_init, param_set, initial_phase, estim_infect,
                            prior_mean_beta, prior_mean_infect_shape, prior_mean_infect_scale, prior_var_beta, prior_var_infect_shape, prior_var_infect_scale, 
                            nrepeat, tun_beta, dist_type = dist_type, model_type = model_type, lik_type = lik_type)

write.table(x = cbind(mcmc.result[,c(1,4,5)], mcmc.result[,4]*mcmc.result[,5], mcmc.result[,1] * mcmc.result[,4] * mcmc.result[,5]),
            file="post_samples.csv",sep=",", 
            col.names=c("beta", "infect_shape", "infect_scale", "mean_infect", "reproduction"), row.names=F)
