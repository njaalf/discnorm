rm(list=ls())
library(discnorm)

require(MPsychoR)
require(lavaan)
data("Rmotivation")
Rmot <- na.omit(Rmotivation[, 1:7])
#Some example distributions that are all standardized:
shape= c(1.5, 2.5, 10)
scale = 1/sqrt(shape)

G3 <- function(x) pgamma(x+shape[1]*scale[1], shape=shape[1], scale=scale[1])
G3flip <- function(x) 1- G3(-x)
G2 <- function(x) pgamma(x+shape[2]*scale[2], shape=shape[2], scale=scale[2])
G2flip <- function(x) 1- G2(-x)
G1 <- function(x) pgamma(x+shape[3]*scale[3], shape=shape[3], scale=scale[3])
G1flip <- function(x) 1- G1(-x)
#QUANTILES
qG3 <- function(x) qgamma(x, shape=shape[1], scale=scale[1])-shape[1]*scale[1]
qG3flip <- function(x) -qG3(1-x)
qG2 <- function(x) qgamma(x, shape=shape[2], scale=scale[2])-shape[2]*scale[2]
qG2flip <- function(x) -qG2(1-x)
qG1 <- function(x) qgamma(x, shape=shape[3], scale=scale[3])-shape[3]*scale[3]
qG1flip <- function(x) -qG1(1-x)


# a non-standardized example, with bounded support:
Fbeta <- function(x) pbeta(x, shape1=2, shape2=3) #example with bounded support

k <- shape[1]
theta <- 2

#List with marginals with various configurations
marginslist <- list(list(F=pnorm, qF=qnorm, sd=1), 
                    list(F=G1), 
                    list(F=Fbeta, support=c(0,1)), 
                    list(F=G3, qF=qG3), 
                    list(F=G1flip, qF=qG1flip), 
                    list(F=G2flip),
                    list(F=function(x) pgamma(x,shape=k, scale=theta) , qF=function(x) qgamma(x,shape=k, scale=theta)))


length(marginslist)

res <- catLSadj(Rmot, marginslist, verbose=TRUE)

f = sem("ext1~~ext2+ext3+ext4+ext5+ext6+ext7", Rmot, ordered=names(Rmot))
gamma.orig <- lavInspect(f, "gamma")
gamma.new <- res[[2]]
library(ggplot2)
qplot(lav_matrix_vec(gamma.orig[8:28, 8:28]),lav_matrix_vec(gamma.new[8:28, 8:28]))+geom_abline()

## EXAMPLE
shape= 2
scale = 1/sqrt(shape)
m1 <- list(F=function(x) pchisq(x, df=1), qF=function(x) qchisq(x, df=1), sd=sqrt(2))
G3 <- function(x) pgamma(x+shape*scale, shape=shape, scale=scale)
G3flip <- function(x) 1- G3(-x)
qG3 <- function(x) qgamma(x, shape=shape, scale=scale)-shape*scale
qG3flip <- function(x) -qG3(1-x)
marginslist <- list(m1, list(F=G3, qF=qG3), list(F=G3flip, qF=qG3flip))
                    
Sigma <- diag(3)
Sigma[Sigma==0] <- 0.5
Sigma
# [,1] [,2] [,3]
# [1,]  1.0  0.5  0.5
# [2,]  0.5  1.0  0.5
# [3,]  0.5  0.5  1.0
set.seed(1)
norm.data <- MASS::mvrnorm(10^5, rep(0,3), Sigma)
colnames(norm.data) <- c("x1", "x2", "x3")
#With normal marginals, the correlation matrix is (approximately)
#Sigma.

#Transform the marginals to follow the elements in marginslist:
nonnorm.data <- data.frame(x1=marginslist[[1]]$qF(pnorm(norm.data[, 1])), 
                     x2=marginslist[[2]]$qF(pnorm(norm.data[, 2])),
                     x3=marginslist[[3]]$qF(pnorm(norm.data[, 3])))
#the transformed data has the following  correlation
cor(nonnorm.data)
#           x1        x2        x3
# x1 1.0000000 0.4389354 0.3549677
# x2 0.4389354 1.0000000 0.4256547
# x3 0.3549677 0.4256547 1.0000000


#The original normal dataset ia fitted perfectly to a factor model with
#the following parameters:
head(standardizedsolution(cfa("F=~ x1+x2+x3", norm.data)),3)
#Extract from output:
#   lhs op rhs est.std    se       z pvalue ci.lower ci.upper
# 1   F =~  x1   0.710 0.002 286.201      0    0.705    0.715
# 2   F =~  x2   0.707 0.002 284.859      0    0.702    0.712
# 3   F =~  x3   0.710 0.002 286.078      0    0.705    0.715

#With non-normal marginals, the factor loadings (and remaining parameters)
#change:
head(standardizedsolution(cfa("F=~ x1+x2+x3", nonnorm.data)),3)
#Extract from output:
#   lhs op rhs est.std    se       z pvalue ci.lower ci.upper
# 1   F =~  x1   0.605 0.003 191.004      0    0.599    0.611
# 2   F =~  x2   0.725 0.003 219.829      0    0.719    0.732
# 3   F =~  x3   0.587 0.003 185.940      0    0.581    0.593

#Discretize the non-normal dataset:
disc.data <- data.frame(x1=cut(nonnorm.data[, 1], breaks= c(-Inf, 0.1, 1, Inf), labels=FALSE), 
                   x2= cut(nonnorm.data[, 2], breaks= c(-Inf, -.7, 0,1, Inf), labels=FALSE),
                   x3=cut(nonnorm.data[, 3], breaks= c(-Inf, -1, 0,1, Inf), labels=FALSE))

#The polychoric correlation is not close to the correlation matrix
#of the non-normal data, but is close to the correlation matrix
#of the -normal- dataset:
lavaan::lavCor(disc.data, ordered=colnames(disc.data))
#    x1    x2    x3   
# x1 1.000            
# x2 0.503 1.000      
# x3 0.506 0.503 1.000

#Compute the adjustments:
adjusted <- catLSadj(disc.data, marginslist, verbose=T )

#The estimated adjusted polychoric correlation matrix is close to the
#correlation matrix of the non-normal data:
adjusted[[1]]
#       x1    x2    x3   
# x1 1.000            
# x2 0.440 1.000      
# x3 0.357 0.427 1.000

# run cat LS without adjustment
head(standardizedsolution(fcat <- cfa("F=~ x1+x2+x3", disc.data, ordered=colnames(disc.data))),3)
#Extract from output:
#   lhs op rhs est.std    se       z pvalue ci.lower ci.upper
# 1   F =~  x1   0.712 0.003 221.553      0    0.705    0.718
# 2   F =~  x2   0.707 0.003 225.835      0    0.701    0.713
# 3   F =~  x3   0.711 0.003 226.119      0    0.705    0.717

#These parameter estimates are close to the parameters of the continuous model for
#normal data, and not to the model parameters obtained from the discretized non-normal dataset
#To get consistent estimates of these parameters
#we need to use the adjusted polychoric correlation.
#To get lavaan to compute the adjusted factor estimates, we need:
sample.th   <- lavInspect(fcat, "sampstat")$th
attr(sample.th, "th.idx") <- lavInspect(fcat, "th.idx")
#the asymptotic covariance matrix of the adjusted polychorics: 
gamma.adj <- adjusted[[2]]
WLS.V.new <- diag(1/diag(gamma.adj))

fcat.adj  <- cfa("F=~ x1+x2+x3", sample.cov=adjusted[[1]],
                  sample.nobs=nrow(disc.data),  sample.th=sample.th,
                  NACOV = gamma.adj, WLS.V=WLS.V.new)
head(standardizedsolution(fcat.adj), 3)
#Extract from output:
#   lhs op rhs est.std    se       z pvalue ci.lower ci.upper
# 1   F =~  x1   0.607 0.003 224.011      0    0.602    0.612
# 2   F =~  x2   0.725 0.003 224.485      0    0.719    0.731
# 3   F =~  x3   0.589 0.002 237.887      0    0.584    0.593

#Closely matches the model parameters obtained with the non-normal dataset
                   



