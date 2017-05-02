
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

  a <- 0
  p = nrow(beta)
  n = ncol(y)
  m = nrow(y)
  N <- n*m
  q = dim(D)[1]

  Dinv <- ginv(D)
  Rinv <- ginv(R)

  cores <- floor(cores)
  if (cores > 1) {
    c2 <- parallel::makeCluster(cores)
    on.exit({
      parallel::stopCluster(c2)
    })
    env <- base::environment()
    parallel::clusterExport(c2, ls(envir = env), envir = env)
  } else {
    c2 <- NULL
  }


  repeat {
    if (a>maxiter||nrm<0.0005) {break}

    #It is the parallel version of ith_ans = lapply(1:n, function(i) {} )
    ith_thread_fn <- function(core_i) {
      positions <- parallel::splitIndices(n, cores)[[core_i]]

      ubfi_sum <- array(0, c(p,p))
      ubsi_sum <- array(0, c(p,1))
      sigma_sum <- 0
      D_sum <- array(0, c(q,q))

      for (i in positions) {
        U_i <- MASS::ginv(
          Dinv + (t(Z[,,i]) %*% Rinv %*% Z[,,i])
        )
        W_i <- Rinv - (Rinv %*% Z[,,i] %*% U_i %*% t(Z[,,i]) %*% Rinv)

        ubfi_sum <- ubfi_sum + t(X[,,i]) %*% W_i %*% X[,,i]
        ubsi_sum <- ubsi_sum + t(X[,,i]) %*% W_i %*% y[,i]

        #Update u (random effect)
        b_i <- U_i %*% t(Z[,,i]) %*% Rinv %*% (y[,i] - X[,,i] %*% beta)
        sigma_sum <- sigma_sum + as.numeric(
          t(y[,i] - X[,,i] %*% beta) %*% W_i %*% (y[,i] - X[,,i] %*% beta)
        )
        D_sum <- D_sum + (1 / sigma) * (b_i %*% t(b_i) + U_i)
      }

      list(
        ubfi_sum = ubfi_sum,
        ubsi_sum = ubsi_sum,
        sigma_sum = sigma_sum,
        D_sum = D_sum
      )
    }

    if (cores == 1) {
      answers <- lapply(1, ith_thread_fn)
    } else {
      answers <- parallel::clusterApply(c2, 1:cores, ith_thread_fn)
    }

    ubfi_total <- array(0,c(p,p)) + Reduce('+', lapply(answers, `[[`, "ubfi_sum"))
    ubsi_total <- array(0,c(p,1)) + Reduce('+', lapply(answers, `[[`, "ubsi_sum"))
    sigma_total <- 0 + Reduce('+', lapply(answers, `[[`, "sigma_sum"))
    D_total <- array(0,c(q,q)) + Reduce('+', lapply(answers, `[[`, "D_sum"))

    # str(list(
    #   ubfi_total = ubfi_total,
    #   ubsi_total = ubsi_total,
    #   sigma_total = sigma_total,
    #   D_total = D_total
    # ))

    final.beta <- ginv(ubfi_total) %*% ubsi_total
    #Final calculations

    final.D <- (1/n) * as.matrix(D_total)
    final.sigma <- as.numeric(
      (1 / N) * as.matrix(sigma_total)[1,1]
    )

    ratio <- final.sigma / sigma
    nrm <- norm(final.beta - beta)

    # str(list(
    #   final.beta = final.beta,
    #   final.D = final.D,
    #   final.sigma = final.sigma
    # ))

    beta = final.beta
    if (median(svd(final.D)$d)<10) {
      D = final.D
    } else {
      D = D
    }
    if (final.sigma > 0) {
      sigma = final.sigma
    } else {
      sigma = sigma
    }
    Dinv <- ginv(D)

    a <- a + 1
  }

  return(list(
    iterations = a,
    beta = beta,
    D = D,
    sigma = sigma,
    norm = nrm
  ))

}
