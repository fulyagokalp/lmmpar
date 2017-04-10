
# library(MASS)
# library(mnormt)
# library(lme4)
# library(nlme)
# library(LaplacesDemon)
# library(mvtnorm)
# library(plyr)
# # install.packages("doParallel")
# library(doParallel)


#' Title
#' (blank)
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
lmm.ep.em <- function(y, X, Z, beta,R,D, cores = 3){

  doParallel::registerDoParallel(cores)

  a <- 0
  aa <- 1
  nrm=100
  maxiter <- 500

  p = nrow(beta)
  n = ncol(y)
  m = nrow(y)
  N <- n*m


  V <- Z%*%D%*%t(Z)+R
  Dinv <- ginv(D)
  Rinv <- ginv(R)
  Vinv <- ginv(V)



  update.beta.first <- array(0,c(p,p,n))
  update.beta.second <- array(0,c(p,1,n))
  update.R <- array(0,c(m,m,n))
  update.D <- array(0,c(p,p,n))
  ETauSquare <- matrix(0,c(n,1))
  ETauSquareMinus <- matrix(0,c(n,1))
  update.lnL1 <- matrix(0,c(n,1))
  update.lnL2 <- matrix(0,c(n,1))




  repeat {
    if (a>maxiter||nrm<0.005) {break}   #nrm<0.0005 idi degistirdim

    #Parallel olmasa, asagidaki kod ith_ans = lapply(1:n, function(i) {} ) olarak yaziliyor.
    ith_ans = plyr::llply(1:n, function(i) {
      #E Step

      Temp1 <- ginv(Dinv+t(Z)%*%Rinv%*%Z)

      u.i <- Temp1%*%t(Z)%*%Rinv%*%(y[,i]-X%*%beta)

      Temp2 <- t(u.i)%*%Dinv%*%u.i

      ETauSquare.i <- 1
      ETauSquareMinus.i <- 1

      Omega <- as.numeric(ETauSquare.i)*Temp1

      #CM Step 1 Fix Ri=RiHat
      #Update beta

      update.beta.first.i <- as.numeric(ETauSquareMinus.i)*(t(X)%*%Rinv%*%X)
      update.beta.second.i <- as.numeric(ETauSquareMinus.i)*(t(X)%*%Rinv%*%(y[,i]-Z%*%u.i))

      #CM Step 2 Fix beta=betaHat

      update.R.i <- as.numeric(ETauSquareMinus.i)*((y[,i]-Z%*%u.i)%*%t(y[,i]-Z%*%u.i)+Z%*%Omega%*%t(Z)+X%*%beta%*%t(beta)%*%t(X)-(X%*%beta%*%t(y[,i]-Z%*%u.i)))

      #CM Step 3

      update.D.i <- Temp1+as.numeric(ETauSquareMinus.i)*(u.i%*%t(u.i))
      update.lnL1.i <- (t(y[,i]-X%*%beta-Z%*%u.i)*as.numeric(ETauSquareMinus.i))%*%(Rinv%*%(y[,i]-X%*%beta-Z%*%u.i))
      update.lnL2.i <- as.numeric(ETauSquareMinus.i)*Temp2

      return(list(
        ub1 = update.beta.first.i,
        ub2 = update.beta.second.i,
        uR  = update.R.i,
        uD  = update.D.i,
        uL1 = update.lnL1.i,
        uL2 = update.lnL1.i
      ))
    }, .parallel = FALSE)

    for(i in 1:n) {
      ith_ans_i = ith_ans[[i]]
      update.beta.first[,,i] = ith_ans_i$ub1
      update.beta.second[,,i] = ith_ans_i$ub2
      update.R[,,i] = ith_ans_i$uR
      update.D[,,i] = ith_ans_i$uD
      update.lnL1[i,] = ith_ans_i$uL1
      update.lnL2[i,] = ith_ans_i$uL2
    }

    #browser()

    #Final calculations
    final.beta <- ginv(apply(update.beta.first, c(1,2), sum))%*%apply(update.beta.second, c(1,2), sum)
    final.R <- (1/N)*(apply(update.R, c(1,2), sum))
    final.D <- (1/n)*(apply(update.D, c(1,2), sum))
    #final.sigma <- (1/n)*sum(sigma)

    nrm <- norm(final.beta-beta)
    beta=final.beta
    if (median(svd(final.R)$d)<10) R=final.R
    else R=R
    if (median(svd(final.D)$d)<10) D=final.D
    else D=D
    V <- Z%*%D%*%t(Z)+R
    Dinv <- ginv(D)
    Rinv <- ginv(R)
    Vinv <- ginv(V)
    final.lnL <- -2*(n/2*log2(2*pi))-(1/2)*(sum(update.lnL1)+log2(det(R))-(n*log2(det(D)))-sum(update.lnL2))
    AICLaplace <- 2*76-2*final.lnL
    #names(out) <- c("iteration","norm","beta","R","D","MeanETauSquare","MeanETauSquareMinus")
    #print(list(final.beta,final.R,final.D))
    a <- a+1
  }

  return(list(
    iterations = a,
    beta = final.beta,
    R = R,
    D = D,
    median_eigen = median(svd(R)$d),
    log_lik = final.lnL,
    AIC = AICLaplace
  ))

}

#
