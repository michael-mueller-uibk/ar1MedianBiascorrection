#--------------------------------------------------------
# Title: Ar1 biascorrection function
# Author: Michael Mueller
# Date: 2023-07-27
# Purpose: Corrects the bias of estimates
# Additional notes:
# Any additional information or notes about the script
#--------------------------------------------------------

# Helper functions --------------------------------------------------------

#' Applying the probabilists' Hermite polynomials (K = 3) to the
#' transformed estimates of phi
#'
#' @param phi A numeric value representing the estimate of phi.
#' @return A numeric vector of length 4 containing the Hermite polynomials.
#'
#' @examples
#' .hermite3_logit(0.5)
#'
#' @keywords internal
.hermite3_logit <- function(phi) {
  # compute logit transformation
  phi <- qlogis((phi + 1) / 2)
  y <- c(1, phi, phi^2 - 1, phi * (phi^2 - 1) - 2 * phi)
  return(y)
}

# Main function -----------------------------------------------------------

#' Corrects the bias of estimates based on the AR(1) process.
#'
#' This function corrects the bias of estimates using the AR(1) process for a
#' specific period 'n'. The 'method' parameter can be used to choose the method
#' to obtain optimal betas from sysdata.RDA.
#'
#' @param phi A numeric value representing the estimate of phi.
#' @param n An integer representing the length of the time series for which the
#'   value of phi was computed.
#' @param method A character string indicating the method that was used to
#'   obtain phi.
#'   Options: "yw" (default) - Yule-Walker method, "burg" - Burg method.
#' @return A numeric value representing the bias corrected estimate 'y'.
#' @export
#'
#' @examples
#' ar1_bias_corr(0.3, 5)
ar1_bias_corr <- function(phi, n, method = "yw") {

  # Sanity checks
  if (!is.numeric(phi) || length(phi) != 1) {
    stop("Argument 'phi' must be a numeric value.")
  }
  if (n <= 4 || n > 100) {
    stop("Argument 'n' must be an integer between 5 and 100")
  }
  if (!method %in% c("yw", "burg")) {
    stop("Invalid 'method' argument. Supported options are 'yw' and 'burg'.")
  }

  # Get optimal betas from sysdata.RDA
  switch(method,
         "yw"   = betas <- ar1MedianBiascorrection:::yule_walker_betas,
         "burg" = betas <- ar1MedianBiascorrection:::burg_betas,
  )
  # Get the betas for the specific period
  b <- betas[[n - 4]]
  # compute unbiased coefficient
  y <- b %*% matrix(.hermite3_logit(phi), 4, 1)
  # compute inverse transformation (g1)
  y <- 1 - 2 / (1 + exp(y))
  # Return the unbiased estimate y
  return(as.numeric(y))
}
