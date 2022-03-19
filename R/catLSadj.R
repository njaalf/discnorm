#' Adjusted polychoric correlation
#' 
#' \code{catLSadj} estimates the underlying correlations assuming bivariate normal copulas, and marginal underlying distributions as provided by the user.  
#' 
#'@param data.df A dataset containing ordinal data.
#'@param marginslist A list of length equal to the number of columns in data.df. Each element in the list
#'specifies a univariate marginal distribution. 
#'Each element must contain a function F, which is the univariate
#'CDF of the marginal. It must accept vectorial input.
#'
#'F is assumed to be continuous, and have support on a interval, which
#'may be all numbers.
#'That is, a random variable X which is F distributed can take
#'any values in an interval (which could be the set of all real numbers).
#
#'In addition, elements "qF", "sd" and "support" may be included.
#'"qF" is to be the quantile function of F.
#'"sd" is to be the standard deviation of F.
#'"support" is the support of the distribution of F.
#
#'If "support" is not included, it is assumed to be all numbers.
#'If qF or sd is not included, they are numerically approximated.
#'
#'For optimal performance, both qF and sd should be passed.
#'If they are not provided, they will be approximated numerically, sometimes
#'at great cost of both precision and execution speed. 
#'For all well-known univariate distributions, both qF and sd are well-known.
#'Implementations for most quantiles are available in R, and standard deviation
#'formulas for most distributions are available on Wikipedia.
#'@param verbose If true, additional information is printed to screen.
#'@return A list of two elements: The adjusted polychoric correlation matrix and its associated asymptotic covariance matrix.
#'@references Steffen Grønneberg & Njål Foldnes (2022) Factor Analyzing Ordinal Items Requires Substantive Knowledge of Response Marginals, Psychological Methods, DOI: 10.1037/met0000495

#'@examples 
#'shape= 2
#'scale = 1/sqrt(shape)
#'m1 <- list(F=function(x) pchisq(x, df=1), qF=function(x) qchisq(x, df=1), sd=sqrt(2))
#'G3 <- function(x) pgamma(x+shape*scale, shape=shape, scale=scale)
#'G3flip <- function(x) 1- G3(-x)
#'qG3 <- function(x) qgamma(x, shape=shape, scale=scale)-shape*scale
#'qG3flip <- function(x) -qG3(1-x)
#'marginslist <- list(m1, list(F=G3, qF=qG3), list(F=G3flip, qF=qG3flip))
#'Sigma <- diag(3)
#'Sigma[Sigma==0] <- 0.5
#'Sigma
#'# [,1] [,2] [,3]
#'# [1,]  1.0  0.5  0.5
#'# [2,]  0.5  1.0  0.5
#'# [3,]  0.5  0.5  1.0
#'set.seed(1)
#'norm.data <- MASS::mvrnorm(10^5, rep(0,3), Sigma)
#'colnames(norm.data) <- c("x1", "x2", "x3")
#'#With normal marginals, the correlation matrix is (approximately)
#'#Sigma.
#'#Transform the marginals to follow the elements in marginslist:
#'nonnorm.data <- data.frame(x1=marginslist[[1]]$qF(pnorm(norm.data[, 1])), 
#'                           x2=marginslist[[2]]$qF(pnorm(norm.data[, 2])),
#'                           x3=marginslist[[3]]$qF(pnorm(norm.data[, 3])))
#'#the transformed data has the following  correlation
#'cor(nonnorm.data)
#'#           x1        x2        x3
#'# x1 1.0000000 0.4389354 0.3549677
#'# x2 0.4389354 1.0000000 0.4256547
#'# x3 0.3549677 0.4256547 1.0000000
#'#The original normal dataset ia fitted perfectly to a factor model with
#'#the following parameters:
#'head(standardizedsolution(cfa("F=~ x1+x2+x3", norm.data)),3)
#'#Extract from output:
#'#   lhs op rhs est.std    se       z pvalue ci.lower ci.upper
#'# 1   F =~  x1   0.710 0.002 286.201      0    0.705    0.715
#'# 2   F =~  x2   0.707 0.002 284.859      0    0.702    0.712
#'# 3   F =~  x3   0.710 0.002 286.078      0    0.705    0.715
#'
#'#With non-normal marginals, the factor loadings (and remaining parameters)
#'#change:
#'head(standardizedsolution(cfa("F=~ x1+x2+x3", nonnorm.data)),3)
#'#Extract from output:
#'#   lhs op rhs est.std    se       z pvalue ci.lower ci.upper
#'# 1   F =~  x1   0.605 0.003 191.004      0    0.599    0.611
#'# 2   F =~  x2   0.725 0.003 219.829      0    0.719    0.732
#'# 3   F =~  x3   0.587 0.003 185.940      0    0.581    0.593
#'#Discretize the non-normal dataset:
#'disc.data <- data.frame(x1=cut(nonnorm.data[, 1], breaks= c(-Inf, 0.1, 1, Inf), labels=FALSE), 
#'                        x2= cut(nonnorm.data[, 2], breaks= c(-Inf, -.7, 0,1, Inf), labels=FALSE),
#'                        x3=cut(nonnorm.data[, 3], breaks= c(-Inf, -1, 0,1, Inf), labels=FALSE))
#'#The polychoric correlation is not close to the correlation matrix
#'#of the non-normal data, but is close to the correlation matrix
#'#of the -normal- dataset:
#'lavaan::lavCor(disc.data, ordered=colnames(disc.data))
#'#    x1    x2    x3   
#'# x1 1.000            
#'# x2 0.503 1.000      
#'# x3 0.506 0.503 1.000
#'#Compute the adjustments:
#'\dontrun{adjusted <- catLSadj(disc.data, marginslist, verbose=T )
#'
#'#The estimated adjusted polychoric correlation matrix is close to the
#'#correlation matrix of the non-normal data:
#'adjusted[[1]]
#'#       x1    x2    x3   
#'# x1 1.000            
#'# x2 0.440 1.000      
#'# x3 0.357 0.427 1.000
#'# run cat LS without adjustment
#'head(standardizedsolution(fcat <- cfa("F=~ x1+x2+x3", disc.data, ordered=colnames(disc.data))),3)
#'#Extract from output:
#'#   lhs op rhs est.std    se       z pvalue ci.lower ci.upper
#'# 1   F =~  x1   0.712 0.003 221.553      0    0.705    0.718
#'# 2   F =~  x2   0.707 0.003 225.835      0    0.701    0.713
#'# 3   F =~  x3   0.711 0.003 226.119      0    0.705    0.717
#'#These parameter estimates are close to the parameters of the continuous model for
#'#normal data, and not to the model parameters obtained from the discretized non-normal dataset
#'#To get consistent estimates of these parameters
#'#we need to use the adjusted polychoric correlation.
#'#To get lavaan to compute the adjusted factor estimates, we need:
#'sample.th   <- lavInspect(fcat, "sampstat")$th
#'attr(sample.th, "th.idx") <- lavInspect(fcat, "th.idx")
#'#the asymptotic covariance matrix of the adjusted polychorics: 
#'gamma.adj <- adjusted[[2]]
#'WLS.V.new <- diag(1/diag(gamma.adj))
#'fcat.adj  <- cfa("F=~ x1+x2+x3", sample.cov=adjusted[[1]],
#'                 sample.nobs=nrow(disc.data),  sample.th=sample.th,
#'                 NACOV = gamma.adj, WLS.V=WLS.V.new)
#'head(standardizedsolution(fcat.adj), 3)
#'#Extract from output:
#'#   lhs op rhs est.std    se       z pvalue ci.lower ci.upper
#'# 1   F =~  x1   0.607 0.003 224.011      0    0.602    0.612
#'# 2   F =~  x2   0.725 0.003 224.485      0    0.719    0.731
#'# 3   F =~  x3   0.589 0.002 237.887      0    0.584    0.593
#'#Closely matches the model parameters obtained with the non-normal dataset}





#'@export
catLSadj <- function(data.df, marginslist, verbose=FALSE) {
  
  
  #Step 1: Standardize the marginals
  
  marginslist <- lapply(marginslist, function(margin) {
    hasF <- hasName(margin, "F")
    hasQf <- hasName(margin, "qF")
    hasSd <- hasName(margin, "sd")
    
    if(!hasF) {
      stop("marginlist has an element that lacks F.")
    }
    F <- margin$F
    
    #Basic sanity check for whether F is a valid CDF:
    if( (F(-Inf) != 0) | (F(Inf) != 1)) {
      stop("A CDF provided for a marginal does not evaluate to zero and one at -Inf and Inf.")
    }
    
    qF <- margin$qF #will be NULL if not provided
    if(!hasQf) {
      if(hasName(margin, "support")) {
        if(!is.numeric(margin$support) | length(margin$support) != 2) {
          stop("support provided for a marginal, but in an incorrect format.")
        }
        qF <- GoFKernel::inverse(F, lower=margin$support[1], upper=margin$support[2])
      } else {
        qF <- GoFKernel::inverse(F)  
      }
    }
    sd <- margin$sd
    if(!hasSd) {
      sd <- calcSd(F, qF)
      if(verbose) cat("Standard deviation approximated to be:", sd, "\n")
    }
    if(sd < 2*.Machine$double.eps) { #we will divide by sd
      stop("sd of a marginal is too close to zero.")
    }
    
    #We transform the CDF F and the quantile function qF in such a way that 
    #they have standard deviation equal to 1:
    #
    #1) CDF
    #The transformed versions are indicated by a tilde.
    #\tilde F(x) = P(X/sigma \leq x) = P ( X \leq x*sigma) = F(x*sigma)
    #
    #2) Quantile function:
    #We calculate
    #\tilde F(x) = y
    #F(x*sigma) = y
    #x*sigma = Q(y)
    #x = Q(y)/sigma.
    #Therefore, tildeQ(y) = Q(y)/sigma.
    
    return(list(F = function(x) F(x*sd), qF = function(x) qF(x)/sd))
  })
  
  #Step 2: Compute adjusted polychoric matrix
  polcorr <- lavaan::lavCor(data.df, ordered=names(data.df), cor.smooth = T)
  
  d <- dim(data.df)[2]
  polcorrAdj <- polcorr
  
  if(verbose) {
    message("Original polychoric correlation matrix:")
    print(polcorr)
  }
  
  for(i in (2:d)) {
    for (j in (1:(i-1))) {
      if(verbose) {
        message("Calculating adjusted polychorics. Computing marginal pair:", i,"-",j)
      }
      polcorrAdj[i,j] <- Psi(polcorr[i,j],
                             F1=marginslist[[i]][[1]], 
                             F2=marginslist[[j]][[1]], 
                             qF1=marginslist[[i]][[2]], 
                             qF2=marginslist[[j]][[2]])
      polcorrAdj[j,i] <- polcorrAdj[i,j]
    }
  }
  
  if(verbose) {
    message("Adjusted polychoric correlation matrix:")
    print(polcorrAdj)
  }
  
  #Step 3: Extract original gamma, and adjust it.
  
  #We fit an arbitrary factor model, just to get lavaan to compute Gamma,
  #which is the same for all models. 
  model <- paste("factor =~", paste(colnames(data.df), collapse = " + "))
  #lavaan will always give a warning, correctly warning against a model that has
  #not converged. We do not care about the model parameters,
  #but use its interface to later extract the (original, non-adjusted) gamma-matrix:
  suppressWarnings(fitMot <- lavaan::cfa(model, data = data.df, ordered = names(data.df), std.lv=T, 
                                         optim.force.converged=TRUE, 
                                         optim.method = "GN", optim.gn.iter.max=1))
  gamma.orig <- lavInspect(fitMot, "gamma")
  gamma.new <- get_newgamma(gamma.orig, polcorr, marginslist)
  
  return(list(polcorrAdj, gamma.new))  
}






#Calculates the standard deviation of a general CDF using
#numerical integraion:
calcSd <- function(F, qF) {
  
  delta <- 10^{-8}
  lowerF <- qF(0)
  upperF <- qF(1)
  
  if(!is.finite(qF(0))) {
    lowerF <- qF(delta)
  }
  
  if(!is.finite(qF(1))) {
    upperF <- qF(1-delta)
  }
  
  successInt <- 0 #variable to keep track of possible failures from 
  #numerical integrations.
  
  #step 1: Expectation:
  
  #Positive part:
  res1Int <- 0
  if(upperF > 0) {
    lowerLimit <- max(0,lowerF)
    upperLimit <- upperF
    res1 <- cubature::hcubature(f=function(x) {
      return(1-F(x[1])) #P(X > x), assuming X continuous
    }, lowerLimit=lowerLimit, upperLimit=upperLimit)
    successInt <- successInt + res1$returnCode
    res1Int <- res1$integral
  }
  #Negative part:
  res2Int <- 0
  if(lowerF < 0) {
    lowerLimit <- lowerF
    upperLimit <- min(upperF,0)
    res2 <- cubature::hcubature(f=function(x) {
      return(-F(x[1]))
    }, lowerLimit=lowerLimit, upperLimit=upperLimit)
    successInt <- successInt + res2$returnCode
    res2Int <- res2$integral
  }
  expInt <- res1Int + res2Int
  
  #step 2: E(X^2)
  #Positive part
  res3Int <- 0
  if(upperF > 0) {
    lowerLimit <- 0
    #If X ~ F, then X^2 has potentially a greater 
    #support than F. If we exceed the support here,
    #1-F(x[1]) is zero. And in the next case
    #(the negative part), F(x[1]) is zero.
    upperLimit <- max(upperF^2,upperF)
    res3 <- cubature::hcubature(f=function(x) {
      return(2*x[1]*(1-F(x[1])))
    }, lowerLimit=lowerLimit, upperLimit=upperLimit)
    successInt <- successInt + res3$returnCode
    res3Int <- res3$integral
  }
  #negative part:
  res4Int <- 0
  if(lowerF < 0) {
    lowerLimit <- min(-lowerF^2, -lowerF)
    upperLimit <- 0
    res4 <- cubature::hcubature(f=function(x) {
      return(-2*x[1]*F(x[1]))
    }, lowerLimit=lowerLimit, upperLimit=upperLimit)
    successInt <- successInt + res4$returnCode
    res4Int <- res4$integral
  }
  expSqInt <- res3Int + res4Int
  varInt <- expSqInt - (expInt^2)
  
  if(varInt < 0) {
    stop("Numerical integration of a standard deviation failed.")
  }
  
  sdInt <- sqrt(varInt)
  if(successInt != 0) { 
    warning("hcubature gave returnCode != 0. Numerical integration may have failed.")
  }
  return(sdInt)
}






#NB: Assumes that the marginal functions are standardized
Psi <- function(rho, F1, F2, qF1, qF2) {
  cop <- copula::normalCopula(rho)
  C <- function(u,v) {
    #cat(u,v,"\n")
    
    #pCopula gives a warning if u or v are 1, and appears
    #to numerically integrate the CDF also there (though
    #gives the correct answer)
    #For these cases, we know by marginal normality
    #the value of the CDF:
    if(u == 1) {
      return(v)
    }
    if(v == 1) {
      return(u)
    }
    
    return(copula::pCopula(copula = cop, u=c(u,v)))
    
    
  }
  CovInt <- function(C,F1,F2, lowerLimit,upperLimit, ...) {
    intFunct <- function(x) {
      return(C(F1(x[1]), F2(x[2])) - F1(x[1])*F2(x[2]))
    }
    return(cubature::hcubature(f=intFunct, lowerLimit=lowerLimit, upperLimit=upperLimit, ...))
  }
  delta <- 10^{-7}
  lowerY1 <- qF1(delta)
  upperY1 <- qF1(1-delta)
  lowerY2 <- qF2(delta)
  upperY2 <- qF2(1-delta)
  lowerLimit <- c(lowerY1,lowerY2)
  upperLimit<- c(upperY1, upperY2)
  
  res <- CovInt(C,F1,F2,lowerLimit, upperLimit)
  if(res$returnCode != 0) {
    warning("hcubature gave returnCode != 0. Numerical integration may have failed.")
    print(res)
  }
  res <- res$integral
  
  #basic sanity check
  if(res > 1 | res < -1) {
    stop("Numerical integration failed.")
  }
  
  return(res)
}

Psi_vect <-Vectorize(Psi, "rho")

#NB: Assumes that the marginal functions are standardized
Der_Psi <- function(rho, F1, F2, qF1, qF2, output = FALSE) {
  phi2rho <- function(u,v) {
    #cat("u,v=",u,v,"\n")
    res <- mnormt::dmnorm(c(u,v), varcov=matrix(c(1,rho,rho,1),2,2)) #default zero mean
    #cat("res=", res, "\n")
    
    #mnormt::dmnorm gives out NaN if u or v are infinite. The density is then equal to zero:
    if(is.nan(res)) { 
      res <- 0
    }
    return(res)
  }
  derInt <- function(phi2rho,F1,F2, lowerLimit,upperLimit, ...) {
    intFunct <- function(x) {
      return(phi2rho(stats::qnorm(F1(x[1])), stats::qnorm(F2(x[2]))) )
    }
    return(cubature::hcubature(f=intFunct, lowerLimit=lowerLimit, upperLimit=upperLimit, ...))
  }
  delta <- 10^{-7}
  lowerY1 <- qF1(delta)
  upperY1 <- qF1(1-delta)
  lowerY2 <- qF2(delta)
  upperY2 <- qF2(1-delta)
  lowerLimit <- c(lowerY1,lowerY2)
  upperLimit<- c(upperY1, upperY2)
  res <- derInt(phi2rho,F1,F2,lowerLimit, upperLimit)
  if(res$returnCode != 0) {
    warning("hcubature gave returnCode != 0. Numerical integration may have failed.")
    print(res)
  }
  if(output == TRUE) {
    cat("H(",rho,") = ",res$integral,"\n")
  }
  return(res$integral)
}

### ASYMPTOTIC COV MATRIX 
getLambdagrad <- function(marginlist, polychoric.mat){
  ## column by column under the diagonal
  derivs <- NULL
  for(j in 1:(length(marginlist)-1)){
    for(i in (j+1):length(marginlist)){
      derivs <- c(derivs, Der_Psi(polychoric.mat[i,j], F1=marginlist[[i]][[1]],
                                  F2=marginlist[[j]][[1]],
                                  qF1=marginlist[[i]][[2]],
                                  qF2=marginlist[[j]][[2]]))
    }
  }
  diag(derivs)
}


get_newgamma <- function(gamma.orig, pol.cor, marginlist){
  dim <- length(marginlist)
  dimgamma <- ncol(gamma.orig)
  subgamma <- gamma.orig[(dimgamma-dim*(dim-1)/2+1):dimgamma, 
                         (dimgamma-dim*(dim-1)/2+1):dimgamma]
  lambda <- getLambdagrad(marginlist, pol.cor)
  ngamma <- gamma.orig
  ngamma[(dimgamma-dim*(dim-1)/2+1):dimgamma, 
         (dimgamma-dim*(dim-1)/2+1):dimgamma] <- lambda %*% subgamma %*% t(lambda)
  ngamma
}


