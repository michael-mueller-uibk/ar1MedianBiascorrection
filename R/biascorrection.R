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


  require("Matrix")
  # TODO Sanity checks
  # ...
  # Read optimal b's depending on the method the ar parameter is estimated
  path_sim <- file.path("data")
  file <- switch(method,
                 "yw" = list.files(path = path_sim, pattern = "opt_b_yw*"),
                 "mle" = list.files(path = path_sim, pattern = "opt_b_mle*"),
                 "burg" = list.files(path = path_sim, pattern = "opt_b_burg*")
  )
  load(file.path(path_sim, file[[1]]))
  # Get the b's for the specific period
  b <- sims[[n - 4]]
  # compute unbiased coefficient
  y <- b %*% matrix(.h3_g(phi), 4, 1)
  # compute inverse transformation (g1)
  y <- 1 - 2 / (1 + exp(y))
  # Return the unbiased estimate y
  return(as.numeric(y))
}