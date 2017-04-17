
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
lmm.ep.em <- function(
  y, X, Z,
  beta, R, D, sigma,
  nu = 3, nrm = 100, maxiter = 500,
  first_parallel = FALSE, second_parallel = FALSE,
  cores = 8
){

  library("doParallel")
  registerDoParallel(cores)
  # library("doMC")
  # registerDoMC(cores)

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
    # ith_ans = plyr::llply(1:n, function(i) {
    ith_ans = plyr::llply(1:8, function(core) {
      if (core > n) stop("make n bigger!")
      positions <- seq(from = core, to = n, by = 8)

      ret <- lapply(positions, function(i) {
        U.i <- ginv(Dinv+(t(Z[,,i])%*%Rinv%*%Z[,,i]))
        U.i
      })
      return(ret)

      #E Step
      # print(paste("E", i))

      # W.i <- Rinv-(Rinv%*%Z[,,i]%*%U.i%*%t(Z[,,i])%*%Rinv)

      #Update beta

      # update.beta.first.i <- t(X[,,i])%*%W.i%*%X[,,i]
      # update.beta.second.i <- t(X[,,i])%*%W.i%*%y[,i]

      # print(paste("E", i, "end"))
      # return(U.i)
      # return(list(
      #   uU  = U.i,
      #   uW  = W.i,
      #   ub1 = update.beta.first.i,
      #   ub2 = update.beta.second.i
      # ))
    }, .parallel = first_parallel)

    ith_ans <- unlist(ith_ans, recursive = FALSE)
    # browser()

    ub1 <- numeric()
    ubfi <- array(0, c(p,p))
    ubsi <- array(0, c(p,1))
    for(i in 1:n) {
      # ith_ans_i = ith_ans[[i]]
      # update.U[,,i] = ith_ans_i$uU
      # update.W[,,i] = ith_ans_i$uW
      U.i <- ith_ans[[i]]
      W.i <- Rinv - (Rinv %*% Z[,,i] %*% U.i %*% t(Z[,,i]) %*% Rinv)
      # browser()
      ubfi <- ubfi + t(X[,,i])%*%W.i%*%X[,,i]
      ubsi <- ubsi + t(X[,,i])%*%W.i%*%y[,i]
      update.U[,,i] = U.i
      update.W[,,i] = W.i
    }
    final.beta <- ginv(ubfi) %*% ubsi

    # final.beta <- ginv(Reduce('+',lapply(ith_ans, function(x) x$ub1)))%*%Reduce('+',lapply(ith_ans, function(x) x$ub2))

    #CM Step 2 Fix beta=betaHat
    #Update sigma
    # ith_ans2 = plyr::llply(1:n, function(i){
    ith_ans2 = plyr::llply(1:8, function(core){
      if (core > n) stop("make n bigger!")
      positions <- seq(from = core, to = n, by = 8)

      sigma_total <- 0
      update.D <- array(0, c(q,q))

      for (i in positions) {
        update.sigma.i <- t(y[,i]-X[,,i]%*%final.beta)%*%update.W[,,i]%*%(y[,i]-X[,,i]%*%final.beta)

        #Update u (random effect)
        u.i <- update.U[,,i]%*%t(Z[,,i])%*%Rinv%*%(y[,i]-X[,,i]%*%final.beta)
        update.D.i <- (1/sigma)*(u.i%*%t(u.i)+update.U[,,i])
        sigma_total <- sigma_total + as.numeric(update.sigma.i)
        update.D <- update.D + update.D.i
      }

      #CM Step 3
      #Update D

      #browser()
      return(list(
        usigma  = sigma,
        uD  = update.D
      ))
    }, .parallel = second_parallel)

    #browser()

    #Final calculations

    final.D <- (1/n)*Reduce('+', lapply(ith_ans2, function(x) x$uD))
    final.sigma <- as.numeric((1/N)*Reduce('+', (lapply(ith_ans2, function(x) x$usigma))))

    ratio <- final.sigma/sigma
    nrm <- norm(final.beta-beta)

    beta = final.beta
    if (median(svd(final.D)$d)<10) {
      D=final.D
    } else {
      D=D
    }
    if (final.sigma > 0) {
      sigma = final.sigma
    } else {
      sigma = sigma
    }
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
