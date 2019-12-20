#' Metropolis-Hastings Algorithm Using Normal Distribution
#'
#' MCMC method for continuous random vector with normal distribution as proposal
#' distribution using a random walk Metropolis algorithm.
#'
#' @param f Density function from which one wants to sample.
#' @param n The numbers of samples one wants to obtain.
#' @param init The initial value vector, which indicates the dimensions.
#' @param stepsize A vector with stepsize[i] being the standard deviation of
#' proposal for updating 'i'th variable.
#' @param burn Times of iterations one wants to omit before recording.
#'
#' @return A "mcmcn" object `list("chain" = chain, "reject" = k/iters, "acpt" = acpt)` with
#' chain storing samples by row, reject being the rejection rate,
#' acpt being whether to be accepted each iter.
#' @export
#'
#' @examples
#' # Generating Multivariate Normal Distribution Samples------------------------
#'
#' # multivariate normal distribution pdf(mu, A)
#' # mu <- c(1, 3)
#' # A <- matrix(c(1, 0.1, 0.1, 1), nrow = 2)
#' f <- pdff("norm", c(1, 3), matrix(c(1, 0.1, 0.1, 1), nrow = 2))
#'
#' # generating random variates using function `mtrp`
#' x.norm <- mtrp(f, 10000, c(3, 3), burn = 0)
#'
#' # exploring the results
#' summary(x.norm)
#' plot(x.norm)
#' qqnorm(x.norm$chain[, 1], main = "QQ plot, 1st variable")
#'
#' # Exploring the Results of Different Initial Values--------------------------
#'
#' # exploring the results of different initial values
#' set.seed(1234)
#' par(mfrow = c(2, 1))
#' x.norm1 <- mtrp(f, 10000, c(-4, 3), burn = 0)
#' x.norm2 <- mtrp(f, 10000, c(1, 3), burn = 0)
#' x.norm3 <- mtrp(f, 10000, c(6, 3), burn = 0)
#' plot(1:500, x.norm1$chain[1:500, 1], type = "l", col = "red",
#'      ylab = "1st variable", ylim = c(-4, 6), main = "first 500 iters")
#' lines(1:500, x.norm2$chain[1:500, 1], col = "green")
#' lines(1:500, x.norm3$chain[1:500, 1], col = "blue")
#' points(rep(1, 3), c(-4, 1, 6), col = c("red", "green", "blue"), pch = 20)
#' plot(9501:10000, x.norm1$chain[9501:10000, 1], type = "l", col = "red",
#'      ylab = "1st variable", ylim = c(-4, 6), main = "last 500 iters")
#' lines(9501:10000, x.norm2$chain[9501:10000, 1], col = "green")
#' lines(9501:10000, x.norm3$chain[9501:10000, 1], col = "blue")
#' par(mfrow = c(1, 1))
#'
#' # Exploring the Results of Different Stepsizes-------------------------------
#'
#' # exploring the results of different stepsizes
#' set.seed(1234)
#' par(mfrow = c(2, 1))
#' for (stepsize in exp(seq(-2, 1, length.out = 4))) {
#'   x.norm <- mtrp(f, 1000, c(5, 3), stepsize = stepsize, burn = 0)
#'   plot(
#'     x.norm$chain[, 1],
#'     type = "l", ylab = "1st variable",
#'     main = paste("stepsize", round(stepsize,2), "rejection rate", round(x.norm$reject,2))
#'   )
#' }
#' par(mfrow = c(1, 1))

mtrp <- function(f,
                 n,
                 init,
                 stepsize = 1,
                 burn = 1000) {
  stopifnot(is.finite(f(init)))
  finit <- f(init)
  stopifnot(burn >= 0)

  # number of variates of the distribution
  nvar <- length(init)

  # number of iterations needed to operate
  iters <- n + burn

  # the variates generated in the ith iteration will be stored in chain[i, ],
  # with chain[1, ] storing the initial value
  acpt <- logical(iters)
  acpt[1] <- TRUE
  chain <- matrix(0, iters, nvar)
  chain[1, ] <- init

  # k: counter of rejections
  k <- 0
  for (i in 2:iters) {
    y <- chain[i - 1, ] + rnorm(nvar) * stepsize

    # garantee y is in S
    while (f(y) / finit < 1e-16) {
      y <- chain[i - 1, ] + rnorm(nvar) * stepsize
    }

    # check whether to update
    if (f(y) >= runif(1) * f(chain[i - 1, ])) {
      acpt[i] <- TRUE
      chain[i, ] <- y
    } else {
      k <- k + 1
      acpt[i] <- FALSE
      chain[i, ] <- chain[i - 1, ]
    }
  }

  if (burn) {
    structure(list(
      chain = chain[-(1:burn), ],
      reject = k / iters,
      acpt = acpt[-(1:burn)]
    ), class = "mcmcn")
  } else {
    structure(list(
      chain = chain,
      reject = k / iters,
      acpt = acpt
    ), class = "mcmcn")
  }
}
