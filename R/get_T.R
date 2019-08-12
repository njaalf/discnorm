
get_T <- function(thresholds){
  T.mat <- NULL
  n.var <- length(thresholds)
  #calculate column dimension of T.mat
  coldim <- 0
  for(j in 1:(n.var-1)){
    th1 <- thresholds[[j]]
    for(i in (j+1):n.var){
      th2 <- thresholds[[i]]
      coldim <- coldim + (length(th1)-sum(is.infinite(th1))+1)*(length(th2)-sum(is.infinite(th2))+1)
    }
  }
  
  # probabilities for var 1 extracted from pi_21 (first)
  th <- thresholds[[1]]
  K1 <- length(th)-sum(is.infinite(th))+1 # there are K categories
  th <- thresholds[[2]]
  K2 <- length(th)-sum(is.infinite(th))+1 # there are K categories
  for(probnr in 1:K1){
    my.row <- rep(0, coldim)
    select <- which(probnr == rep(1:K1, K2))
    my.row[select] <- 1
    T.mat <- rbind(T.mat, my.row)
  }
  
  # for the others
  counter <- 0
  for(var in 2:n.var){#extract from pi_var,1
    th <- thresholds[[var]]
    K <- length(th)-sum(is.infinite(th))+1 # there are K categories
    for(probnr in 1:K){
      #task is to find row in T that gives this probability
      my.row <-rep(0, coldim)
      for(s in (counter+1):(counter+K1))
        my.row[s] <- 1
      counter <- counter + K1
      T.mat <- rbind(T.mat, my.row)
    }
  }
  T.mat
}
