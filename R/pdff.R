#' Getting the pdf of some distribution
#'
#' Getting the pdf of some distribution, where the Const is ignored for MCMC method.
#'
#' @param type Type of the distribution one wants to get.
#' c("norm", "mixnorm") is provided.
#' @param ... Additional arguments for type.
#'
#' @return The pdf of `type`.
#' @export
#'
#' @examples
#' # norm --------------------------------------------------------------------
#' f <- pdff("norm")
#' # mixnorm -----------------------------------------------------------------
#' p <- 0.5
#' mu1 <- c(0, 0)
#' mu2 <- c(6, 6)
#' A1 <- matrix(c(1, 0.1, 0.1, 1), nrow = 2)
#' A2 <- matrix(c(1, -0.1, -0.1, 1), nrow = 2)
#' g <- pdff("mixnorm", p, mu1, mu2, A1, A2)
pdff <- function(type = c("norm", "mixnorm"), ...) {
  stopifnot(!is.null(type))
  type <- match.arg(type)
  UseMethod("pdff", structure(list(), class = type))
}


#' @rdname pdff
#'
#' @param mu,A Parameters for normal distribution,
#' mu being the mean vecter while A being the covariance matrix
#' @export
pdff.norm <- function(type, mu=0, A=1) {
  if ((nvar <- length(mu)) == 1L) {
    stopifnot(length(A) == 1L)
    function(x) {
      stopifnot(length(x) == nvar)
      exp(-(x - mu) ^ 2 / (2 * A))
    }
  } else {
    stopifnot(nvar == nrow(A))
    A.inv <- solve(A)
    function(x)
      as.vector(exp(-0.5 * t(x - mu) %*% A.inv %*% (x - mu)))
  }
}

#' @rdname pdff
#'
#' @param p,m1,mu2,A1,A2 Parameters for mixture normal distribution.
#' The pdf can be written as p*pdf1 + (1-p)*pdf2
#' @export
pdff.mixnorm <- function(type, p, mu1, mu2, A1, A2) {
  stopifnot(length(p) == 1L)
  stopifnot((nvar <- length(mu1)) == length(mu2))
  if (nvar == 1L) {
    stopifnot(length(A1) == 1L && length(A2) == 1L)
    function(x)
      p * A1 ^ -0.5 * exp(-(x - mu1) ^ 2 / (2 * A)) +
      (1 - p) * A2 ^ -0.5 * exp(-(x - mu2) ^ 2 / (2 * A2))
  } else {
    stopifnot(nvar == nrow(A1) && nvar == nrow(A2))
    A1.det <- det(A2)
    A2.det <- det(A2)
    A1.inv <- solve(A1)
    A2.inv <- solve(A2)
    function(x)
      as.vector(
        p * A1.det ^ (-1 / 2) * exp(-0.5 * t(x - mu1) %*% A1.inv %*% (x - mu1)) +
          (1 - p) * A2.det ^ (-1 / 2) * exp(-0.5 * t(x - mu2) %*% A2.inv %*% (x - mu2))
      )
  }
}
