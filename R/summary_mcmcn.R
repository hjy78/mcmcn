#' @export

summary.mcmcn <- function(x) {
  nvar = ncol(x$chain)
  if(!identical(x$alpha, NULL)) cat("Rejection rate:", x$alpha, "\n\n")
  print(summary(x$chain))
  cat("\nCovarience Matrix:\n")
  print(var(x$chain))
  return(invisible())
}
