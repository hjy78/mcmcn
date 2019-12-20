#' Metropolis-Hastings Algorithm Using Exponential Distribution
#'
#' MCMC method for continuous random vector with exponential distribution as
#' proposal distribution using an independence sampler algorithm.
#'
#' @param f Density function from which one wants to sample.
#' @param n The numbers of samples one wants to obtain.
#' @param init The initial value vector, which indicates the dimensions.
#' @param a The left support boundary vector. if a numeric number is provided,
#' it will be broadcast to a vector with the same length as init.
#' @param rate A vector with rate[i] being the rate of
#' proposal for updating 'i'th variable. if a numeric number is provided,
#' it will be broadcast to a vector with the same length as init.
#' @param burn Times of iterations one wants to omit before recording.
#'
#' @return A "mcmcn" object `list("chain" = chain, "reject" = k/iters, "acpt" = acpt)` with
#' chain storing samples by row, reject being the rejection rate,
#' acpt being whether to be accepted each iter.
#' @export
#'
#' @examples
#' # the Exp(1) Distribution Left-Truncated at 1--------------------------------
#' (i.e., with support [1,+Infinity))
#'
#' # pdf f(x) = exp(1 - x)
#' f <- function(x) ifelse(x >= 1, exp(-x), 0)
#'
#' # generating random variates using function `mtrp_exp`
#' x.exp1 <- mtrp_exp(f, 10000, 1, a = 1, burn = 0)
#'
#' # exploring the results
#' summary(x.exp1)
#' plot(x.exp1, "c")
#' hist(x.exp1$chain, freq = FALSE, breaks = 50,
#'      main = "Histogram of Samples",
#'      xlab = "X", xlim = c(0,4), ylim = c(0,  exp(1)))
#' curve(exp(1 - x), col = "red", add = TRUE)

mtrp_exp <- function(f,
                     n,
                     init,
                     a = 0,
                     rate = 1,
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
    y <- a + rexp(nvar) / rate

    # garantee y is in S
    while (f(y) / finit < 1e-16) {
      y <- a + rexp(nvar) / rate
    }

    # check whether to update
    if (f(y) * exp(-sum(rate * (chain[i - 1, ] - a))) >= runif(1) *
        f(chain[i - 1, ]) * exp(-sum(rate * (y - a)))) {
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
