
getdisc <- function(X,  thresholds){# must start at zero for sirt::polychorc
  z <- X
  my.data <- sapply(1:ncol(X), function(i)  arules::discretize(z[,i], method = "fixed", labels=F, breaks=thresholds[[i]]))
  colnames(my.data ) <- paste("d", 1:ncol(X), sep="")
  sapply(data.frame(my.data), function(col) col-min(col))
}