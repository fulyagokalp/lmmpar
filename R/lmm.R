
# https://github.com/fulyagokalp/lmmpar
# devtools::install_github("fulyagokalp/lmmpar")

# library(MASS)
# library(mnormt)
# library(lme4)
# library(nlme)
# library(LaplacesDemon)
# library(mvtnorm)
# library(plyr)
# # install.packages("doParallel")
# library(doParallel)
# library(psych)


#' Title
#'
#' description
#' @param y desc1
#' @param X desc1
#' @param Z desc1
#' @param beta desc1
#' @param R desc2
#' @param D desc3
#' @param cores desc4
#' @importFrom MASS ginv
#' @importFrom mnormt rmnorm
#' @importFrom stats median
#' @export
# ' @example
# ' #example code here
# ' lmm.ep.em()

#Reference paper: Schafer, J.L.S., 1998, some improved procedures for linear mixed models
#Ui, Wi and beta is calculated in parallel form. Then,sigma and D are calculated with final beta
lmm.ep.em <- function(y, X, Z, beta, R, D, sigma, cores = 3, nu=3, nrm=100, maxiter=500, first_parallel=FALSE, second_parallel=FALSE){
  #first_parallel <- TRUE

  doParallel::registerDoParallel(cores)

  a <- 0
  p = nrow(beta)
  n = ncol(y)
  m = nrow(y)
  N <- n*m
  q = dim(D)[1]

  Dinv <- ginv(D)
  Rinv <- ginv(R)

  update.beta.first <- array(0,c(p,p,n))
  update.beta.second <- array(0,c(p,1,n))
  update.U <- array(0, c(q, q, n))
  update.W <- array(0, c(m, m, n))
  update.sigma <- matrix(0, n, 1)
  update.D <- array(0,c(q,q,n))


  repeat {
    if (a>maxiter||nrm<0.0005) {break}

    #It is the parallel version of ith_ans = lapply(1:n, function(i) {} )
    ith_ans = plyr::llply(1:n, function(i) {
      #E Step

      U.i <- ginv(Dinv+(t(Z[,,i])%*%Rinv%*%Z[,,i]))
      W.i <- Rinv-(Rinv%*%Z[,,i]%*%U.i%*%t(Z[,,i])%*%Rinv)

      #Update beta

      update.beta.first.i <- t(X[,,i])%*%W.i%*%X[,,i]
      update.beta.second.i <- t(X[,,i])%*%W.i%*%y[,i]

      return(list(
        uU  = U.i,
        uW  = W.i,
        ub1 = update.beta.first.i,
        ub2 = update.beta.second.i
      ))
    }, .parallel = first_parallel)

    for(i in 1:n) {
      ith_ans_i = ith_ans[[i]]
      update.U[,,i] = ith_ans_i$uU
      update.W[,,i] = ith_ans_i$uW
    }

    final.beta <- ginv(Reduce('+',lapply(ith_ans, function(x) x$ub1)))%*%Reduce('+',lapply(ith_ans, function(x) x$ub2))
    #CM Step 2 Fix beta=betaHat
    #Update sigma
    ith_ans2 = plyr::llply(1:n, function(i){
      update.sigma.i <- t(y[,i]-X[,,i]%*%final.beta)%*%update.W[,,i]%*%(y[,i]-X[,,i]%*%final.beta)

      #Update u (random effect)
      u.i <- update.U[,,i]%*%t(Z[,,i])%*%Rinv%*%(y[,i]-X[,,i]%*%final.beta)
      #CM Step 3
      #Update D

      update.D.i <- (1/sigma)*(u.i%*%t(u.i)+update.U[,,i])
      #browser()
      return(list(
        usigma  = update.sigma.i,
        uD  = update.D.i

      ))
    }, .parallel = second_parallel)

    #browser()

    #Final calculations

    final.D <- (1/n)*Reduce('+', lapply(ith_ans2, function(x) x$uD))
    final.sigma <- as.numeric((1/N)*Reduce('+', (lapply(ith_ans2, function(x) x$usigma))))

    ratio <- final.sigma/sigma
    nrm <- norm(final.beta-beta)

    beta = final.beta
    if (median(svd(final.D)$d)<10) D=final.D
    else D=D
    if (final.sigma > 0) sigma = final.sigma
    else sigma = sigma
    Dinv <- ginv(D)

    a <- a+1
  }

  return(list(
    iterations = a,
    beta = beta,
    D = D,
    sigma = sigma,
    norm = nrm
  ))

}
