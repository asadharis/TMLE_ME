# The main functions for calculating the
# Average Treatment effect and related quantities.
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
  term1 <- sapply(L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf1)
  })
  term2 <- sapply(L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf2)
  })
  term3 <- sapply(L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf3)
  })
  
  mean(term1 * Y * X) + 
    mean(term2) - 
    mean(term3 * X)
  
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
  term1 <- sapply(L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf1)
  })
  term2 <- sapply(L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf2)
  })
  term3 <- sapply(L_err, FUN = function(z_t){
    calE_transform(z_t = z_t, lim_t = lambda,
                   lim_z = lambda, char_noise = NULL, 
                   sd_e = 1, myf3)
  })
  
  mean(term1 * Y * (1-X)) + 
    mean(term2) - 
    mean(term3 * (1-X))
  
}

# 
est_ate <- function(Y, X, L_err, 
                     MU, TAU, lambda = 1){
  est_psi1(Y,X,L_err, MU,TAU, lambda = lambda) - 
    est_psi0(Y,X,L_err, MU,TAU, lambda = lambda)
}

# Estimate psi0, psi1 when there is no measurement error
# using the doubly robust estimator. 
est_ate_noError <- function(Y, X, L, 
                            MU, TAU) {
  psi1 <- mean(X * Y * 1/TAU(L)) +
    mean(MU(L,1)) - 
    mean(X * MU(L,1)/TAU(L))
  
  psi0 <- mean((1-X) * Y * 1/(1-TAU(L)) ) +
    mean(MU(L,0)) - 
    mean((1-X) * MU(L,0)/(1-TAU(L)))
  ate <- psi1- psi0
  return(c("psi0" = psi0, "psi1" = psi1, "ate" = ate))
}
