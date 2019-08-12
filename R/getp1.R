getp1sam <- function(my.data){
  unlist(apply(my.data,2,  tabulate))/nrow(my.data)
}

#the dot means that we concatenate many variables
getp1pop <- function(thresholds){
  unlist(lapply(thresholds, function(x) diff(pnorm(x))))
}

