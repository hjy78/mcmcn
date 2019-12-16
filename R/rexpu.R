#' @param n Number of observations.
#' @param a,b Lower and upper limits of the uniform distribution.
#' @param rate Rate of exponential distribution.

rexpu <- function(n, a = 0, b = 0, rate = 1) {
  x.unif <- runif(n, a, b)
  x.exp <- b + rexp(n, rate)
  p <- rbinom(n, 1, 1 / (rate * (b - a) + 1))
  x <- ifelse(p, x.exp, x.unif)
  return(x)
}
