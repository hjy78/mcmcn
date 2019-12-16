#' @export

plot.mcmcn <- function(x) {
  n = nrow(x$chain)
  nvar = ncol(x$chain)
  if(nvar >= 2) {
    par(mfrow = c(1, 2))
    plot(x$chain[1:100, 1], x$chain[1:100, 2], type = "b", pch = 20,
         main = "first 100 iters", xlab = "1st variable", ylab = "2nd variable")
    plot(x$chain[(n-99):n, 1], x$chain[(n-99):n, 2], type = "b",
         pch = 20, main = "last 100 iters", xlab = "1st variable", ylab = "2nd variable")
  }
  par(mfrow = c(2, 1))
  plot(1:500, x$chain[1:500, 1], type = "l",
       xlab = "iters", ylab = "1st variable",
       main = "first 500 iters")
  plot((n-499):n, x$chain[(n-499):n, 1], type = "l",
       xlab = "iters", ylab = "1st variable",
       main = "last 500 iters")
  par(mfrow = c(1, 1))
  return(invisible())
}
