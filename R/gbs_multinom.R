#' Gibbs Algorithm for Generating Multinom Distribution
#'
#' MCMC method for generating vectors with multinom distribution using Gibbs sampler.
#'
#' @param n The numbers of samples one wants to obtain.
#' @param size,prob Parameters for multinom distribution. size being a integer,
#' say N, specifying the total number of objects that are put into d boxes
#' in the typical multinomial experiment. prob being numeric non-negative vector
#' of length d, specifying the probability for the d classes;
#' should be normalized to sum 1 manually before providing.
#' @param init The initial value vector, which indicates the dimensions.
#' @param burn Times of iterations one wants to omit before recording.
#'
#' @return A "mcmcn" object `list("chain" = chain)` with chain storing samples by row.
#' @export
#'
#' @examples
#' # Generating Multinom Distribution Samples-----------------------------------
#'
#' # setting some parameters
#' set.seed(123)
#' prob <- c(0.2, 0.1, 0.1, 0.1, 0.1, 0.1,
#'           0.05, 0.05, 0.05, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02)
#' p <- sample(prob)
#'
#' # generating random varaites using function `gbs_multinom`
#' # and exploring the results
#' par(mfrow = c(2, 1))
#' for (size in 10^(2:6)) {
#'   init <- size * p
#'   x.multinom <- gbs_multinom(10000, size, prob, init, burn = 0)
#'   plot(x.multinom[1:2000, 1], type = "l", ylab = "1st variable",
#'        ylim=c(0.95*0.2*size, 1.05*0.2*size),
#'        main = paste("size =", size, "first 2000 iters"))
#'   plot(density(x.multinom[-(1:1000), 2]), xlab = "occurs",
#'        xlim = c(0.08*size, 0.22*size),
#'        main = paste("size =", size))
#'   lines(density(x.multinom[-(1:1000), 1]))
#' }
#' par(mfrow = c(1, 1))

gbs_multinom <- function(n,
                      size,
                      prob,
                      init,
                      burn = 1000) {
  stopifnot(sum(prob) == 1)
  stopifnot(length(prob) == length(init))
  stopifnot(sum(init) == size)
  stopifnot(burn >= 0)

  # number of variates of the distribution
  nvar <- length(prob)

  # number of iterations needed to operate
  iters <- n + burn

  # the variates generated in the ith iteration will be stored in chain[i, ],
  # with chain[1, ] storing the initial value
  chain <- matrix(0, iters, nvar)
  chain[1, ] <- init
  for (i in 2:iters) {
    x <- sample(1:nvar, 2)
    s <- sum(chain[i - 1, c(x[1], x[2])])
    chain[i, -c(x[1], x[2])] <- chain[i - 1, -c(x[1], x[2])]
    chain[i, x[1]] <- rbinom(1, s, prob[x[1]] / (prob[x[1]] + prob[x[2]]))
    chain[i, x[2]] <- s - chain[i, x[1]]
  }

  ifelse(burn, return(structure(list("chain" = chain[-(1:burn), ]), class = "mcmcn")),
         return(structure(list("chain" = chain), class = "mcmcn")))
}
