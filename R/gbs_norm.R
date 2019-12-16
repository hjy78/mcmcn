#' Gibbs Algorithm for Generating Multivariate Normal Distribution
#'
#' MCMC method for generating vectors with multivariate normal distribution
#' using Gibbs sampler.
#'
#' @param n The numbers of samples one wants to obtain.
#' @param mu,A Parameters for normal distribution, mu being the mean vecter while
#' A being the covariance matrix
#' @param init The initial value vector, which indicates the dimensions.
#' @param burn Times of iterations one wants to omit before recording.
#'
#' @return A "mcmcn" object `list("chain" = chain)` with chain storing samples by row.
#' @export
#'
#' @examples
#' # Generating Multivariate Normal Distribution Samples------------------------
#'
#' # provide some parameters
#' mu <- c(1,3)
#' A <- matrix(c(1, 0.1, 0.1, 1), nrow = 2)
#'
#' # generating random variates using function `gibbs_norm`
#' x.norm <- gbs_norm(mu, A, 10000, mu, burn = 0)
#'
#' # exploring the results
#' cat("The sample mean is", colMeans(x.norm), "\n")
#' cat("The sample covariance matrix is\n")
#' cov(x.norm)
#' par(mfrow = c(1, 2))
#' plot(x.norm[1:100, 1], x.norm[1:100, 2], type = "b", pch = 20,
#'      main = "first 100 iters", xlab = "1st variable", ylab = "2nd variable")
#' plot(x.norm[9901:10000, 1], x.norm[9901:10000, 2], type = "b",
#'      pch = 20, main = "last 100 iters", xlab = "1st variable", ylab = "2nd variable")
#' par(mfrow = c(2, 1))
#' plot(1:500, x.norm[1:500, 1], type = "l", ylab = "1st variable",
#'      main = "first 500 iters")
#' plot(9501:10000, x.norm[9501:10000, 1], type = "l", ylab = "1st variable",
#'      main = "last 500 iters")
#' par(mfrow = c(1, 2))
#' plot(x.norm, pch = 20, cex = 0.5, main = "Sample Distribution",
#'      xlab = "1st variable", ylab = "2nd variable")
#' qqnorm(x.norm[, 1], main = "QQ plot, 1st variable")
#' par(mfrow = c(1, 1))

gbs_norm <- function(n,
                     mu,
                     A,
                     init,
                     burn = 1000) {
  stopifnot(length(mu) == nrow(A) && nrow(A) == ncol(A))
  stopifnot(burn >= 0)

  # number of variates of the distribution
  nvar <- length(mu)

  # Calculate the inverse of A[-i, -i] for i in 1:nvar and store the results in A.inv
  # A.inv will be used for calculating the conditional distribution
  A.inv <- array(0, dim = c(nvar - 1, nvar - 1, nvar))
  for (i in 1:nvar) {
    A.inv[, , i] <- solve(A[-i, -i])
  }

  # number of iterations needed to operate
  iters <- n + burn

  # the variates generated in the ith iteration will be stored in chain[i, ],
  # with chain[1, ] storing the initial value
  chain <- matrix(0, iters, nvar)
  chain[1, ] <- init
  for (i in 2:iters) {
    j <- sample(1:nvar, 1)
    chain[i, -j] <- chain[i - 1, -j]
    chain[i, j] <- rnorm(1,
                         mu[j] + A[j, -j] %*% A.inv[, , j] %*% (chain[i, -j] - mu[-j]),
                         (A[j, j] - A[j, -j] %*% A.inv[, , j] %*% A[-j, j]) ^ (1 / 2))
  }

  ifelse(burn, return(structure(list("chain" = chain[-(1:burn), ]), class = "mcmcn")),
         return(structure(list("chain" = chain), class = "mcmcn")))
}
