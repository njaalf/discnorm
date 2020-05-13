#this is the implied theorethical probabilities based on normality
pc_PI <- function(rho, th.y1, th.y2) {
  nth.y1 <- length(th.y1); nth.y2 <- length(th.y2)
  pth.y1 <- stats::pnorm(th.y1);  pth.y2 <- stats::pnorm(th.y2)
  
  # catch special case: rho = 0.0
  if(rho == 0.0) {
    rowPI <- diff(c(0,pth.y1,1))
    colPI <- diff(c(0,pth.y2,1))
    PI.ij <- outer(rowPI, colPI)
    return(PI.ij)
  }
  
  # prepare for a single call to pbinorm
  upper.y <- rep(th.y2, times=rep.int(nth.y1, nth.y2))
  upper.x <- rep(th.y1, times=ceiling(length(upper.y))/nth.y1)
  #rho <- rep(rho, length(upper.x)) # only one rho here
  
  BI <- pbivnorm::pbivnorm(x=upper.x, y=upper.y, rho=rho)
  #BI <- pbinorm1(upper.x=upper.x, upper.y=upper.y, rho=rho)
  dim(BI) <- c(nth.y1, nth.y2)
  BI <- rbind(0, BI, pth.y2, deparse.level = 0)
  BI <- cbind(0, BI, c(0, pth.y1, 1), deparse.level = 0)
  
  
  # get probabilities
  nr <- nrow(BI); nc <- ncol(BI)
  PI <- BI[-1L,-1L] - BI[-1L,-nc] - BI[-nr,-1L] + BI[-nr,-nc]
  
  # all elements should be strictly positive
  PI[PI < .Machine$double.eps] <- .Machine$double.eps
  
  PI
}

