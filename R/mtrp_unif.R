#' Metropolis-Hastings Algorithm Using Uniform Distribution
#'
#' MCMC method for continuous random vector with uniform distribution as proposal
#' distribution using an independence sampler algorithm.
#'
#' @param f Density function from which one wants to sample.
#' @param n The numbers of samples one wants to obtain.
#' @param init The initial value vector, which indicates the dimensions.
#' @param a,b The left and right support boundary vector. if a numeric number
#' is provided, it will be broadcast to a vector with the same length as init.
#' @param burn Times of iterations one wants to omit before recording.
#'
#' @return A "mcmcn" object `list("chain" = chain, "reject" = k/iters, "acpt" = acpt)` with
#' chain storing samples by row, reject being the rejection rate,
#' acpt being whether to be accepted each iter.
#' @export
#'
#' @examples
#' # the Exp(1) distribution right-truncated at 1-------------------------------
#' # i.e., with support [0,1]
#'
#' # pdf f(x) = exp(1 - x) / (e - 1)
#' f <- function(x) ifelse(x >= 0 && x <=1, exp(-x), 0)
#'
#' # generating random variates using function `mtrp_unif`
#' x.exp2 <- mtrp_unif(f, 10000, 1, a = 0, b = 1, burn = 0)
#'
#' # exploring the results
#' summary(x.exp2)
#' plot(x.exp2, "c")
#' hist(x.exp2$chain, freq = FALSE, main = "Histogram of Samples", xlab = "X")
#' curve(exp(1 - x) / (exp(1) - 1), from = 0, to = 1, col = "red", add = TRUE)

mtrp_unif <- function(f,
                      n,
                      init,
                      a,
                      b,
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
  sc <- b - a
  for (i in 2:iters) {
    y <- a + sc * runif(nvar)

    # garantee y is in S
    while (f(y) / finit < 1e-16) {
      y <- a + sc * runif(nvar)
    }

    # check whether to update
    if (f(y) >= runif(1) * f(chain[i - 1, ])) {
      acpt[i] <- TRUE
      chain[i, ] <- y
    } else {
      acpt[i] <- FALSE
      k <- k + 1
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
