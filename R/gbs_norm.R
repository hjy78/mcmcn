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
#' # mu <- c(1,3)
#' # A <- matrix(c(1, 0.1, 0.1, 1), nrow = 2)
#'
#' # generating random variates using function `gibbs_norm`
#' x.norm <- gbs_norm(10000, c(1,3), matrix(c(1, 0.1, 0.1, 1), nrow = 2), c(1,3), burn = 0)
#'
#' # exploring the results
#' summary(x.norm)
#' plot(x.norm)
#' qqnorm(x.norm$chain[, 1], main = "QQ plot, 1st variable")

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

  if (burn)
    structure(list(chain = chain[-(1:burn), ]), class = "mcmcn")
  else
    structure(list(chain = chain), class = "mcmcn")
}
