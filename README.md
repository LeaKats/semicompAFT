# semicompAFT
The script for running the estimation procedure is called "semicompAFT". Before running the estimation procedure, the file "rcpp_based_functions.R" in this directory should be placed in the same working directory.

Here we list the arguments of the estimation procedures:

## Estimation With Frailty - "estimation_with_frailty_f"
Function Arguments:

X01 - covariate matrix for the healthy-illness transition.

X02 - covariate matrix for the healthy-death transition.

X12 - covariate matrix for the illness-death transition.

V - minimum of disease time, death time (from healthy state) or censoring time.

W - minimum of death time (from illness state) and censoring time. If the person died or was censored directly from healthy state, W should be set to 0.

delta1 - 1 if V = illness time, 0 otherwise.

delta2 - 1 if V = death time, 0 otherwise.

delta3 - 1 if W = death time, 0 otherwise.

zeta_beta - zeta value in the bandwidth computation used for betas estimation. Defaults to 65.

zeta_h - zeta value in the bandwidth computation used for baseline hazard functions estimation. Defaults to 0.01.

initial_sigma - starting value for sigma. Defaults to 5.

conv_betas_bound - tolerance for betas' convergence. Defaults to 0.00001.

conv_Hs_bound - tolerance for cumulative hazard functions' convergence. Defaults to 0.0001.

conv_sigma_bound - tolerance for sigma's convergence. Defaults to 0.0001.

stop_iter_num - maximal number of iterations. Defaults to 1000.

B - number of required weighted bootstrap samples. Defaults to 100. Weights are randomly sampled from the standard exponential distribution.

print - T if the iterative process should be printed, F  otherwise.

## Estimation Without Frailty - "estimation_without_frailty_f"

Function Arguments:

X01 - covariate matrix for the healthy-illness transition.

X02 - covariate matrix for the healthy-death transition.

X12 - covariate matrix for the illness-death transition.

V - minimum of disease time, death time (from healthy state) or censoring time.

W - minimum of death time (from illness state) and censoring time. If the person died or was censored directly from healthy state, W should be set to 0.

delta1 - 1 if V = illness time, 0 otherwise.

delta2 - 1 if V = death time, 0 otherwise.

delta3 - 1 if W = death time, 0 otherwise.

zeta_beta - zeta value in the bandwidth computation used for betas estimation. Defaults to 65.

zeta_h - zeta value in the bandwidth computation used for baseline hazard functions estimation. Defaults to 0.01.

B - number of required weighted bootstrap samples. Defaults to 100. Weights are randomly sampled from the standard exponential distribution.

print - T if the iterative process should be printed, F  otherwise.


## Example for running the estimation procedure
Example for running the estimation procedure based on a simulated data found in the script named "example_with_simulated_data".
