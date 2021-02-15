# Main functions for carrying out integration.
# We now have a function for the main transformations which 
# we will use. 

# First we have the transform which takes something from the 
# "no-measurement error domain" to the 
# "measurement error domain".

# The double integrand involved in calculating this function.
#   Args:
#       x: a 2xN matrix, first row corresponds variable t and 
#          second to varaible z (as defined in current version).
#       z_t: A scalar z_tilde value at which we will calculate the 
#             transformation.
#       char_noise: the characteristic function of the noise. 
#       FUN: the vectorized function we wish to transform
#       ...: Other arguments for function, FUN.
calE_integrand <- function(x, z_t = 0, 
                           char_noise = NULL, 
                           FUN, ...) {

  Re(matrix(FUN(x[2,], ...) *
           (1/char_noise(x[1,])) *
           exp(-1i * x[1,]*(x[2,] - z_t)) * 
           (2*pi)^(-1), ncol = ncol(x)))
}

# The function to calculate the cal_E transformation.
# 
#   Args:
#       z_t: Scalar value at which to calculate the transform.
#       lim_t: Limits of integration for the t variable. 
#       lim_z: Limits of integration for the z variable.
#       char_noise: the characteristic function of the noise. If NULL,
#                   noise is assumed to be Gaussian.
#       sd_e: If char_noise == NULL, SD of Gaussian noise.
#       FUN: the vectorized function we wish to transform
#       ...: Other arguments for function, FUN.
calE_transform <- function(z_t = 0, lim_t = 3,
                           lim_z = lim_t, 
                           char_noise = NULL, sd_e = 1, 
                           FUN, ...) {
  if(is.null(char_noise)) {
    char_noise <- function(t) char_gaus(t, sd = sd_e)
  }
  
  require(cubature)
  ans <- hcubature(f = calE_integrand,
                   lowerLimit = c(-lim_t, -lim_z),
                   upperLimit = c(lim_t,lim_z),
                   z_t=z_t, char_noise = char_noise, 
                   FUN = FUN, ...,
                   vectorInterface = TRUE)$integral
  ans
}


############### EXTRA TEST FUNCTIONS ####################

# integrand_cdf <- function(t, z_t = 0, 
#                            char_noise = NULL, 
#                            gamma = 1, ...) {
#   
#   Re((1i)/(t*char_noise(t)) * exp(-1i * t *(gamma - z_t)))
#   
# }
# 
# transform_cdf <- function(z_t = 0, lim_t = 3,
#                           lim_z = lim_t, 
#                           char_noise = NULL, sd_e = 0.1, 
#                           ...) {
#   if(is.null(char_noise)) {
#     char_noise <- function(t) char_gaus(t, sd = sd_e)
#   }
#   
#   require(cubature)
#   ans <- hcubature(f = calE_integrand,
#                    lowerLimit = c(-lim_t, -lim_z),
#                    upperLimit = c(lim_t,lim_z),
#                    z_t=z_t, char_noise = char_noise, 
#                    FUN = FUN, ...,
#                    vectorInterface = TRUE)$integral
#   ans
# }
# 


################# TESTING ###############################
# alpha <- 1.5
# calE_transform(.1, FUN =  function(x){1/expit(alpha*x)})
