# Some helper functions.

# The expit function.
expit <- function(x) {
  1/(1 + exp(-x))
}

# Characteristic function of mean 0 Normal RV
char_gaus <- function(t, sd = 1) {
  exp(-1 * (sd)^2*t^2/2)
}

