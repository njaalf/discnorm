#' Bootstrap test for underlying non-normality in ordinal data.
#' 
#'@param my.data a dataset containing ordinal data.
#'@param B number of bootstrap samples.
#'@return p-value associated with the underlying normality hypothesis.
#'@export
bootTest <- function(my.data, B=1000){
  #sirt needs minimum value to be zero
  my.data <- sapply(data.frame(my.data), function(col) col-min(col,na.rm=T))
  P.hat <- sirt::polychoric2(my.data, cor.smooth=TRUE, use_pbv=FALSE)$rho
  thresholds.hat <- lapply(data.frame(my.data), function(x) unique(c(-Inf, pc_th(x), Inf)))
  d <- ncol(P.hat)
  indices <- NULL #has nrow*(nrow-1)/2 elements
  for(i2 in 1:(d-1)){
    for( i1 in (i2+1):d){
      indices <- rbind(indices, c(i1,i2))
    }
  }
  
  #bootnorm returns NA if original data can  not use computeT, and a negative number if there are many NA's among the bootstrap samples
  Tstat <- tryCatch(computeT(my.data, indices),error=function(w) { NULL})
  if(is.null(Tstat)){
    stop("Error: could not compute test statistic in original sample.\n")
    return(NA) 
  }
  cat(" B = ")
  TstatBoot <- rep(0, B)
  for(i in (1:B)) {
    if(!(i %% 100))
      cat( i, "\n")
    norm.sample <- MASS::mvrnorm(n=nrow(my.data), mu=rep(0, ncol(P.hat)), Sigma=P.hat)
    boot.sample <- data.frame(getdisc(norm.sample, thresholds.hat))
    TstatBoot[i] <- tryCatch(computeT(boot.sample, indices), error=function(w) { NA})
  }
  naprop <- sum( is.na(TstatBoot))/length(TstatBoot)
  
  if(naprop  > 0.5){
    stop("Error: could not compute test statistic in more than 50% of the bootstrap samples.\n")
    return(NA)
  }
  
  mean(Tstat < TstatBoot, na.rm=T)
}
