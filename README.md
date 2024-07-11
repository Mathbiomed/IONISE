# IONISE - Inference Of Non-markovIan SEIR Model

## Overview
IONISE is a user-friendly computational package developed for the gamma model approach, specifically targeting the SEIR (Susceptible-Exposed-Infectious-Removed) model. The package is implemented in R software version 4.3.2 and aims to estimate epidemiological parameters from sole confirmed case data using gamma-distributed latent and infectious periods.

## Getting Started

### Prerequisites
- R software version 4.3.2 or higher
- Additional R packages:
  - `deSolve`
  - `MASS`
  - `mvtnorm`
  - `dplyr`
  - `invgamma`
  - `ramcmc`

### Modules
Source functions are modularized with four modules below:
- `IONISE_module_basic_computing.R`: This module contains basic functions for the packages, such as calculating probability density functions and convolutions between density functions.
- `IONISE_module_diff_eqs.R`: This module contains numerical solvers for differential equations describing the SEIR and extended SEIR models with delays. 
- `IONISE_module_MH.R`: This module contains functions for the Metropolis-Hastings (MH) algorithm used in Markov chain Monte Carlo (MCMC) sampling.
- `IONISE_module_MCMC.R`: This module contains functions performing MCMC sampling.

### Installation
1. Download the IONISE package and example CSV files from the [GitHub repository](#https://github.com/Mathbiomed/IONISE) (link will be activated after acceptance).
2. Save confirmed case data as a CSV file with the name `input_data.csv` in the folder containing the R code of the main function, `IONISE_main.R`.

## Usage

1. Open `IONISE_main.R` in your R environment.
2. Set the folder containing the R code as the working directory using the function `setwd()`.
3. Install and load the required packages:
   ```R
   install.packages(c("deSolve", "MASS", "mvtnorm", "dplyr", "invgamma", "ramcmc"))
   library(deSolve)
   library(MASS)
   library(mvtnorm)
   library(dplyr)
   library(invgamma)
   library(ramcmc)
   ```
4. Compile the necessary function libraries using the line `source("IONISE_source.R")`.
5. Import the input data CSV file using the line `read.csv("input_data.csv")`.

### Simple Manual
For more details, refer to a step-by-step manual provided in the relevant reference (Supplementary Note 2). 

1. Set the variable `initial_phase` to `TRUE` if the data comes from the initial phase of disease spread; otherwise, set it to `FALSE`.
2. Specify whether to estimate the infectious period distribution (`estim_infect`).
   - If `TRUE`, the code will estimate shape and scale parameters (`k_(τ_I)`, `s_(τ_I)`).
   - If `FALSE`, specify the mean and variance of the infectious period distribution (`infect_mean`, `infect_var`).
3. Specify the parameters for the MCMC procedure. The variable `nrepeat` represents the number of MCMC iterations, and the variable `tun_beta` represents the proposal variance for the Metropolis-Hastings algorithm for sampling the `β`. The default values for `nrepeat` and `tun_beta` are 10,000 and 0.1, respectively.
4. Specify the distribution type for the delays (i.e., latent and infectious periods) in the model using the variable `dist_type`. The variable is one of `exp`, `gamma`, `invgamma`, `lognormal`, and `weibull`. 
5. Specify the type of the likelihood function for the inference using the variable `lik_type`. The variable is one of `poisson` and `gaussian`.
6. Set up hyperparameters for prior distributions:
   - Transmission rate (`prior_mean_beta`, `prior_var_beta`).
   - Infectious period (`prior_mean_infect_shape`, `prior_mean_infect_scale`, `prior_var_infect_shape`, `prior_var_infect_scale`).
   - Default values are non-informative priors (mean=1, variance=1,000,000).
7. Specify the type of a model used for the inference using the variable “model_type.” If the variable is `SEIR`, the SEIR model is used for the inference. If the variable is `Extended_SEIR`, the extended model is used for the inference. See the relevant reference to see descriptions of the SEIR model and the extended model.
8. In the code section `PART 2-1` or `PART 2-2`, specify the initial values of the compartment and the values of the parameters in the chosen model. 
9. Run the code section named `PART 3: Perform MCMC`.

## Output

Once `IONISE_main.R` is properly executed, a CSV file named `post_samples.csv` is created. This file contains posterior samples of β, `k_(τ_I)`, `s_(τ_I)`, `μ_(τ_I)`, and `R` in the first, second, third, fourth, and fifth columns, respectively. Note that the burn-in period and thinning for the MCMC iterations are not automatically applied in the code. To obtain independent posterior samples from the MCMC iterations, users may need to apply the burn-in period and thinning.

## Expected Runtimes

- **Installation:** Less than 10 seconds.
- **Execution:** Approximately 1 hour per 10,000 MCMC iterations for the example data when it the code tested with MacBook Air (2023) running macOS Sonoma (14.4.1) with processor M2 8-Core with R version 4.3.2 (2023-10-31).

## Reference

1. To be updated.
