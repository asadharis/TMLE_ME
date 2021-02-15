

# A basic data generating mechanism for estimating
# the average treatment effect. 
#
# Varying: alpha, sigma, sigma_e
#
# Y|X,L = 1 + X - 1.5*L + Normal(0, sigma)
# Pr(X=1| L) = expit(alpha*L)
# L = Normal(0,1)
# L^* = L + Normal(0, sigma_e);
#
# EY(1) = 1 + 1 - E(1.5*L) = 2
# EY(0) = 1 + 0 - E(1.5*L) = 1
# ATE = EY(1) - EY(0) = 1.
#
gen.data <- function(n=100, alpha = 0.5, sigma = 0.2,
                     sigma_e = 1) {
  L <- rnorm(n)
  L_err <- L + rnorm(n, sd = sigma_e)
  
  # First define the TAU function.
  TAU <- function(L) {
    expit(alpha * L)
  }
  pr_x <- TAU(L)
  X <- rbinom(n, size = 1, prob = pr_x)
  
  # Then define the MU function
  MU <- function(L, X) {
    1 + X - 1.5*L
  }
  Y <- MU(L,X) + rnorm(n, sd = sigma)
  return(list("Y"  = Y, "X" = X, "L" = L, 
                    "L_err" = L_err, "PrX" = pr_x,
                    "MU" = MU, "TAU" = TAU))
}
