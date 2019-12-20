#' Summary for an mcmcn Object
#'
#' Summary function for displaying the information of a "mcmcn" class,
#' including rejection rate, univariate mean quantile, multivariate covariance matrix, etc.
#'
#' @param x An mcmcn object to summary.
#'
#' @export
#'
#' @examples
#' f <- pdff("norm", c(1, 3), matrix(c(1, 0.1, 0.1, 1), nrow = 2))
#' x.norm <- mtrp(f, 10000, c(3, 3), burn = 0)
#' summary(x.norm)

summary.mcmcn <- function(x) {
  nvar = ncol(x$chain)
  if(!is.null(x$reject)) cat("Rejection rate:", x$reject, "\n\n")
  print(summary(x$chain))
  if (nvar == 1L) {
    cat("\nVarience:", var(x$chain), "\n")
  } else {
    cat("\nCovarience Matrix:\n")
    print(var(x$chain))
  }
  invisible()
}
