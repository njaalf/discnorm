# compute thresholds
pc_th <- function(Y, freq=NULL, prop=NULL) {
  # my line:
  Y <- Y-min(Y, na.rm=T)+1# start at 1, required by tabulate
  if(is.null(prop)) {
    if(is.null(freq)) freq <- base::tabulate(Y)
    prop <- freq / sum(freq)
  }
  prop1 <- prop[-length(prop)]
  stats::qnorm(base::cumsum(prop1))
}