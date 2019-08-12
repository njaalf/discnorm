# compute thresholds
pc_th <- function(Y, freq=NULL, prop=NULL) {
  # my line:
  Y <- Y-min(Y)+1# start at 1, required by tabulate
  if(is.null(prop)) {
    if(is.null(freq)) freq <- tabulate(Y)
    prop <- freq / sum(freq)
  }
  prop1 <- prop[-length(prop)]
  qnorm(cumsum(prop1))
}