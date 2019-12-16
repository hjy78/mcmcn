#' Metropolis-Hastings Algorithm Using Exponential-Uniform Mixture Distribution
#'
#' MCMC method for continuous random vector with exponential-uniform mixture
#' distribution as proposal distribution using Metropolis-Hastings algorithm.
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
#' @return A "mcmcn" object `list("chain" = chain, "alpha" = k/iters, "acpt" = acpt)` with
#' chain storing samples by row, alpha being the rejection rate,
#' acpt being whether to be accepted each iter.
#' @export
#'
#' @examples
#' # the Exponential-Uniform Mixture--------------------------------------------
#'
#' # suppose X ~ exponential-uniform mixture with parameters(a = 0, b = 1, rate = 1)
#' # the plot below shows the different pdfs of X and Exp(1)
#' curve(0.5 * exp(1-x), from = 1, xlim = c(0, 3), ylim = c(0, 1), ylab = "y")
#' lines(c(0, 1), c(0.5, 0.5))
#' curve(exp(-x), col = "red", add = TRUE)
#' legend("topright", c("Exp(1)", "X"), lty = c(1, 1), col = c("red", "black"))
#'
#' # the Exponential-Gamma mixture----------------------------------------------
#' # (Suppose that the rate parameter Λ has Gamma ( r, beta ) distribution and X has
#' # Exp(Λ) distribution. That is, X|Λ = λ ∼ f_X (x|λ) = λe^{−λy}.)
#'
#' # pdf f(x) = r * beta^r / (x + beta)^{r+1}
#' r <- 1
#' beta <- 2
#' f <- function(x){
#'   ifelse(x >= 0, (x + beta)^(-(r+1)), 0)
#' }
#'
#' # generating random variates using function `mtrp_expu`
#' x.pareto <- mtrp_expu(f, 10000, 1, rate = 0.5, burn = 0)
#'
#' # exploring the results
#' cat("The sample rejection rate is", x.pareto$alpha)
#' par(mfrow = c(2, 1))
#' plot(1:500, x.pareto$chain[1:500], type = "l",
#'      main = "first 500 iters", ylab = "X")
#' plot(9501:10000, x.pareto$chain[9501:10000], type = "l",
#'      main = "last 500 iters", ylab = "X")
#' par(mfrow = c(1, 1))
#' hist(x.pareto$chain, freq = FALSE, main = "Histogram of Samples", xlab = "X")
#' curve(r * beta ^ r * (x + beta)^(-(r+1)), from = 0, col = "red", add = TRUE)

mtrp_expu <- function(f,
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

  # check a and rate
  stopifnot(length(a) == 1 || length(a) == nvar)
  if (length(a) == 1)
    a <- rep(a, nvar)
  stopifnot(length(rate) == 1 || length(rate) == nvar)
  if (length(rate) == 1)
    rate <- rep(rate, nvar)

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

  # generating chain[i, ]
  for (i in 2:iters) {
    x <- chain[i - 1, ]
    y <- numeric(nvar)
    for (j in 1:nvar) {
      y[j] <- rexpu(1, a[j], x[j], rate[j])
    }

    # garantee y is in S
    while (f(y) / f(init) < 1e-16) {
      for (j in 1:nvar) {
        y[j] <- rexpu(1, a[j], x[j], rate[j])
      }
    }

    # proposal distribution g(y|x)
    g <- function(y, x) {
      return(prod(ifelse(
        y <= x, 1 / (x - a + 1 / rate), exp(-rate * (y - x)) / (x - a + 1 / rate)
      )))
    }

    # check whether to update
    if (f(y) * g(x, y) >= runif(1) * f(x) * g(y, x)) {
      acpt[i] <- TRUE
      chain[i, ] <- y
    } else {
      k <- k + 1
      acpt[i] <- FALSE
      chain[i, ] <- x
    }
  }

  ifelse(burn, return(structure(list(
    "chain" = chain[-(1:burn), ], "alpha" = k / iters, "acpt" = acpt[-(1:burn)]
  ), class = "mcmcn")), return(structure(
    list("chain" = chain, "alpha" = k / iters, "acpt" = acpt), class = "mcmcn"
  )))
}
