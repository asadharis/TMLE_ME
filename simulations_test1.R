############################################################
############################################################
############################################################

# First pass of simulations for testing our methodology.

##### SIMULATION SETUP ####
# In this simulation we run a first test of our ate 
# estimation functions. Our goal is to estimate the 
# Average Treament Effect (ATE) and related mean values 
# of counter factuaals using the TMLE framework for data
# with measurement error. 
# 
# We generate a linear model with binary treatment X and
# continuous confounder L. We observe an error in our 
# confounder term. For details of data generating mechanism
# see 'datagen_functions.R'.
# 
# Our method relies on using the deconvolution problem 
# originally done for density/distribution estimation. 
# Thus we calculate (numerically) a fourier type integral
# however, calculation over the entire real line/hypercube is 
# unstable and thus we introduce a tuning parameter lambda. 
# 
# In this first simulation we would like to observe how our 
# estimator performs as a function of lambda. 
#
# For this simple setup we assume the true propensity
# score or outcome model is known. In practice, these would
# have to be estimated either from the data or from external
# data/sources. 
############################################################
############################################################
############################################################


##### LOAD LIBRARIES #######
#
source("helper_functions.R")
source("integral_functions.R")
source("datagen_functions.R")
source("ate_functions.R")
#

# Estimates ATE for sequence of lambda values.
#
# Args:
#   seed: Seed value for data generation
#   n: Sample size
#   sigma: The SD of the residual (not the error in data)
#   sigma_e: The SD of the measurement error
#   lambdas: A sequence of lambda values to run
#   ncores: The number of cores to use for the simulation
run_sim <- function(seed = 1, n = 100, sigma = 1, 
                    sigma_e = 0.5, lambdas = 1:3,
                    verbose = FALSE) {
  require(cubature)
  # seed = 1; n = 100; sigma = 1;
  # sigma_e = 0.5; lambdas = 1:3;ncores = 8

  # Generate data for fixed seed value. 
  set.seed(seed)
  dat <- gen.data(n = n, sigma = sigma, sigma_e = sigma_e)
  psi0 <- psi1 <- numeric(length(lambdas))
  
  exports_par <- c("expit", "est_psi0", "est_psi1", 
                   "char_gaus", "calE_integrand",
                   "calE_transform", "dat", 
                   "lambdas")
  
  nlam <- length(lambdas)
  for(i in 1:length(lambdas)) {
    if(verbose) {
      cat("Iteration number:", i, "\n")
    }

    psi0[i] <-   est_psi0(dat$Y, dat$X, dat$L_err, 
                          dat$MU, dat$TAU, 
                          lambda = lambdas[i])
    psi1[i] <-  est_psi1(dat$Y, dat$X, dat$L_err, 
                         dat$MU, dat$TAU, 
                         lambda = lambdas[i])
  }
  ate <- psi1 - psi0
  pars_noError <- est_ate_noError(dat$Y, dat$X, dat$L, 
                                  dat$MU, dat$TAU)
  fin <- data.frame("lam" = lambdas, 
                    "psi0" = psi0, "psi1" = psi1, 
                    "ate" = ate, 
                    "ate_noError" = pars_noError[3],
                    "ate_TRUE" = 1)
  
  dirname <- paste0("sim1/n=", n,
                    "/sigma=", sigma,
                    "/sigma_e=", sigma_e)
  filename <- paste0(dirname, "/",seed, ".RData")
  
  if(dir.exists(dirname)) {
    save(list = c("fin"), file = filename)
  } else {
    dir.create(dirname, recursive = TRUE)
    save(list = c("fin"), file = filename)
  }
  
}


args <-  commandArgs(T)
seed <- as.numeric(args[[1]])
n <- as.numeric(args[[2]])
sigma<- as.numeric(args[[3]])
sigma_e<- as.numeric(args[[4]])


# seed <- 1
# n <- 100
# sigma <- 1
# sigma_e <- 0.5

lambdas <- seq(0.5, 3, length = 100)

run_sim(seed, n, sigma, sigma_e,lambdas, verbose = TRUE)

q(save = "no")


