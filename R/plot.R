#' Plot for an mcmcn Object
#'
#' Plot to display Markov chain and sample distribution for "mcmcn" class (sample).
#'
#' @param x An mcmcn object to plot.
#' @param type Charactor vector choosen from c("c", "d").
#' "c" short for chain, "d" short for distribution.
#' @param index Integer vector specifying dimensions to be displayed.
#'
#' @section Details:
#' When the sample is 1-dimensional, the parameter index does not work.
#' Draw the sample by passing "c" for type (the first 500 times and the last
#' 500 times) ~ index line chart, and draw the sample histogram by passing "d"
#' for type. And kernel density estimates.
#'
#' When the sample is 2 dimensions and above, the parameter index defaults to
#' NULL. By passing "c" for type, draw a line chart of the first two dimensions
#' (first 500 times and last 500 times) of the sample, and pass "d" for type.
#' To plot a scatter plot of the first two dimensions of the sample (used to
#' show the sample). If the parameter index is set as an index array,
#' for example, index = c (2, 3), then pass "c" to draw the samples of the
#' second dimension and the third dimension (the first 500 times and the last
#' 500 times) ~ index respectively Line chart, pass in "d"
#' To plot sample histograms and kernel density estimates
#' for the 2nd and 3rd dimensions, respectively.
#'
#' @export
#' @examples
#' f <- pdff("norm", c(1, 3), matrix(c(1, 0.1, 0.1, 1), nrow = 2))
#' x.norm <- mtrp(f, 10000, c(3, 3), burn = 0)
#' plot(x.norm)

plot.mcmcn <- function(x, type=c("c", "d"), index = NULL) {
  type <- match.arg(type, several.ok = TRUE)
  n = nrow(x$chain)
  nvar = ncol(x$chain)
  if ("c" %in% type) {
    if (nvar == 1L) {
      par(mfrow = c(2, 1))
      plot(1:500, x$chain[1:500, 1], type = "l",
           xlab = "iters", ylab = "X",
           main = "first 500 iters")
      plot((n-499):n, x$chain[(n-499):n, 1], type = "l",
           xlab = "iters", ylab = "X",
           main = "last 500 iters")
    } else if (is.null(index)) {
      par(mfrow = c(1, 2))
      plot(x$chain[1:100, 1], x$chain[1:100, 2], type = "b", pch = 20,
           main = "first 100 iters", xlab = "1st variable", ylab = "2nd variable")
      plot(x$chain[(n-99):n, 1], x$chain[(n-99):n, 2], type = "b",
           pch = 20, main = "last 100 iters", xlab = "1st variable", ylab = "2nd variable")
    } else {
      par(mfrow = c(2, 1))
      for (i in index) {
        plot(1:500, x$chain[1:500, i], type = "l",
             xlab = "iters", ylab = paste0(i, "th variable"),
             main = "first 500 iters")
        plot((n-499):n, x$chain[(n-499):n, i], type = "l",
             xlab = "iters", ylab = paste0(i, "th variable"),
             main = "last 500 iters")
      }
    }
    par(mfrow = c(1, 1))
  }
  if ("d" %in% type) {
    if(nvar == 1L) {
      hist(x$chain, freq = FALSE, main = "Histogram of Samples", xlab = "X")
      lines(density(x$chain), col = "red")
    } else if (is.null(index)) {
      plot(x$chain, pch = 20, cex = 0.5, main = "Sample Distribution",
           xlab = "1st variable", ylab = "2nd variable")
    } else {
      for (i in index) {
        hist(x$chain[, i], freq = FALSE, xlab = paste0(i, "th variable"),
           main = paste0("Histogram of ", i, "th variable"))
        lines(density(x$chain[, i]), col = "red")
      }
    }
  }
  invisible()
}
