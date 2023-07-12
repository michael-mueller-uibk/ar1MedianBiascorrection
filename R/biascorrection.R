# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# Helper functions --------------------------------------------------------

# Applying the probabilists' Hermite polynomials (K = 3) to the
# transformed estimates of phi
.h3_g <- function(phi) {
  # compute logit transformation
  phi <- qlogis((phi + 1) / 2)
  y <- c(1, phi, phi^2 - 1, phi * (phi^2 - 1) - 2 * phi)
  return(y)
}

# Main function -----------------------------------------------------------

ar1_bias_corr <- function(phi, n, method = "yw") {
  #' Title
  #'
  #' @param cutoff
  #'
  #' @return
  #' @export
  #'
  #' @examples

  # TODO Sanity checks
  # ...
  # Get optimal betas from sysdata.RDA
  switch(method,
         "yw"   = betas <- ar1MedianBiascorrection:::yule_walker_betas,
         "burg" = betas <- ar1MedianBiascorrection:::burg_betas,
  )
  # Get the betas for the specific period
  b <- betas[[n - 4]]
  # compute unbiased coefficient
  y <- b %*% matrix(.h3_g(phi), 4, 1)
  # compute inverse transformation (g1)
  y <- 1 - 2 / (1 + exp(y))
  # Return the unbiased estimate y
  return(as.numeric(y))
}
