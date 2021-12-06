
computeT <- function(my.data, indices) {
  suppressWarnings(Phat <- sirt::polychoric2(my.data, cor.smooth=F, use_pbv = 2)$rho)
  #suppressWarnings(Phat <- lavaan::lavCor(data.frame(my.data), 
              #                            cor.smooth=FALSE)
  thresholds.hat <- lapply(data.frame(my.data), function(x) unique(c(-Inf, pc_th(x), Inf)))
    
  eStat <- apply(indices, MARGIN=1, function(X) 
  {
    th1 <- thresholds.hat[[X[1]]]; th1 <- th1[!is.infinite(th1)]
    th2 <- thresholds.hat[[X[2]]]; th2 <- th2[!is.infinite(th2)] #-c(1,thLength)
    residual <- lavaan::lav_matrix_vec(t(table(my.data[,X[1]], my.data[,X[2]]))/nrow(my.data) - t(pc_PI(Phat[X[1], X[2]], th1, th2)) )
    sum(residual^2)
  }
  )
  #T is now the sum of sqares of eStat.
  return(sum(eStat))
}
