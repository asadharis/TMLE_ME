# The main functions for calculating the
# Average Treatment effect and related quantities.
#
##### LOAD LIBRARIES #######
#
source("helper_functions.R")
source("integral_functions.R")
source("datagen_functions.R")
#

# Estimates psi_1 = EY_1 
#
# Args:
#   Y: Vector of response
#   X: Vector of exposure/treatment
#   L_err: Vector of covariate measurements.
#   MU: A function, defining estimate of outcome model,
#       first argument is covariate value second argument is 
#       indicator value, i.e MU(L,1) or MU(L,0)
#   TAU: A function, defining the propensity score model
#   lambda: defines the threshold value
est_psi1 <- function(Y, X, L_err, 
                     MU, TAU, lambda = 1) {
  # Define the functions which we need to transform
  # 1. 1/tau
  myf1 <- function(z) {
    1/TAU(z)
  }
  # 2. mu
  myf2 <- function(z) {
    MU(z,1)
  }
  # 3. mu/tau -> Note: For future we can use 
  #             convolution property of fourier transform.
  myf3 <- function(z) {
    MU(z,1)/TAU(z)
  }
  term1 <- sapply(dat$L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf1)
  })
  term2 <- sapply(dat$L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf2)
  })
  term3 <- sapply(dat$L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf3)
  })
  
  mean(term1 * dat$Y * dat$X) + mean(term2) - mean(term3 * dat$X)
  
}

# Estimates psi_0 = EY_0 
#
# Args:
#   SAME AS est_psi1

est_psi0 <- function(Y, X, L_err, 
                     MU, TAU, lambda = 1) {
  # Define the functions which we need to transform
  # 1. 1/tau
  myf1 <- function(z) {
    1/(1 - TAU(z))
  }
  # 2. mu
  myf2 <- function(z) {
    MU(z,0)
  }
  # 3. mu/tau -> Note: For future we can use 
  #             convolution property of fourier transform.
  myf3 <- function(z) {
    MU(z,0)/(1 - TAU(z))
  }
  term1 <- sapply(dat$L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf1)
  })
  term2 <- sapply(dat$L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf2)
  })
  term3 <- sapply(dat$L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf3)
  })
  
  mean(term1 * dat$Y * (1-dat$X)) + 
    mean(term2) - 
    mean(term3 * (1-dat$X))
  
}

# 
est_ate <- function(Y, X, L_err, 
                     MU, TAU, lambda = 1){
  est_psi1(Y,X,L_err, MU,TAU, lambda = lambda) - 
    est_psi0(Y,X,L_err, MU,TAU, lambda = lambda)
}

