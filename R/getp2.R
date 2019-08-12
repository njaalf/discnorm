
getp2sam <- function(my.data){
  pvector <- NULL
  nrow <- ncol(my.data)
  
  for(i2 in 1:(nrow-1)){
    for( i1 in (i2+1):nrow){
      pvector <- c(pvector, lav_matrix_vec(t(table(my.data[,i1], my.data[,i2]))))
    }}
  pvector/nrow(my.data)
}

getp2pop <- function(P,  thresholds){
  pvector <- NULL
  nrow <- nrow(P)
  for(i2 in 1:(nrow-1)){
    for( i1 in (i2+1):nrow){
      rho <- P[i1, i2]
      th1 <- thresholds[[i1]]; th1 <- th1[!is.infinite(th1)]
      th2 <- thresholds[[i2]]; th2 <- th2[!is.infinite(th2)]
      pvector <- c(pvector, lav_matrix_vec(t(pc_PI(rho, th1, th2))))
    }}
  pvector 
}
