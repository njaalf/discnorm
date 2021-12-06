#' Bootstrap test for discretized normality
#' 
#' \code{bootTest} is a bootstrap test for whether an ordinal dataset is consistent with being
#' a discretization of a multivariate normal dataset. 
#' 
#'@param my.data A dataset containing ordinal data. Must contain only integer values. 
#'@param B Number of bootstrap samples.
#'@param verbose If true, bootstrap progress is printed to the console.
#'@return p-value associated with the underlying normality hypothesis.
#'@references Njål Foldnes & Steffen Grønneberg (2019) Pernicious Polychorics: The Impact and Detection of Underlying Non-normality, Structural Equation Modeling: A Multidisciplinary Journal, DOI: 10.1080/10705511.2019.1673168

#
#'@examples 
#'set.seed(1)
#'norm.data <- MASS::mvrnorm(300, m=rep(0,3), 
#'Sigma=cov(MASS::mvrnorm(15, mu=rep(0,3), Sigma=diag(3))))
#'disc.data <- apply(norm.data,2,  cut, 
#'breaks = c(-Inf, 0,1, Inf), labels=FALSE)# normal data discretized
#'pvalue <- bootTest(disc.data, B=500)
#'#no support for underlying non-normality
#'@export
bootTest <- function(my.data, B=1000, verbose=TRUE){
  #sirt needs minimum value to be zero
  my.data <- sapply(data.frame(my.data), function(col) col-min(col,na.rm=T))
  #P.hat <- sirt::polychoric2(my.data, cor.smooth=TRUE, use_pbv=FALSE)$rho
  P.hat <- lavaan::lavCor(data.frame(my.data), ordered=colnames(my.data), cor.smooth = TRUE)
  
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
  if(verbose)
    cat("Progress 0% ")
  TstatBoot <- rep(0, B)
  for(i in (1:B)) {
    if(verbose & !(i %% ceiling(B/10)))
      cat( 100*i/B, "% ", sep="")
    norm.sample <- MASS::mvrnorm(n=nrow(my.data), mu=rep(0, ncol(P.hat)), Sigma=P.hat)
    boot.sample <- data.frame(getdisc(norm.sample, thresholds.hat))
    TstatBoot[i] <- tryCatch(computeT(boot.sample, indices), error=function(w) { NA})
  }
  if(verbose)
    cat("\n")
  naprop <- sum( is.na(TstatBoot))/length(TstatBoot)
  
  if(naprop  > 0.5){
    stop("Error: could not compute test statistic in more than 50% of the bootstrap samples.\n")
    return(NA)
  }
  
  mean(Tstat < TstatBoot, na.rm=T)
}
