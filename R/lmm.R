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
# library(matrixcalc) #Is positive definite package
# library(magrittr)


#' Title
#'
#' description
#' @param y matrix of responses with observations/subjects on column and repeats for each observation/subject on rows. It is (m $\times$ n) dimensional.
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
#Function is updated for stacked vector and matrices.
lmm.ep.em <- function(
  y, X, Z, subject,
  beta, R, D, sigma,
  nrm = 100, maxiter = 500,
  cores = 8
){

  a <- 0
  p = nrow(beta)
  n = length(unique(subject))
  m = nrow(y) / n
  N <- n*m
  q = dim(D)[1]

  Dinv <- ginv(D)
  Rinv <- ginv(R)

  cores <- floor(cores)
  if (cores > 1) {
    doParallel::registerDoParallel(cores)
    #c2 <- parallel::makeCluster(cores)
    # c2 <- parallel::makeCluster(cores, type="FORK")
    # on.exit({
    #   parallel::stopCluster(c2)
    # })
    # env <- base::environment()
    # parallel::clusterExport(c2, ls(envir = env), envir = env)
  } else {
    c2 <- NULL
  }

  # colnames(y) <- "Y"
  # colnames(X) <- paste("X", 1:ncol(X), sep = "")
  # colnames(Z) <- paste("Z", 1:ncol(Z), sep = "")
  # y <- cbind(y, subject = subject)
  # X <- cbind(X, subject = subject)
  # Z <- cbind(Z, subject = subject)
  #
  # y_mem <- bigmemory::as.big.matrix(y)
  # X_mem <- bigmemory::as.big.matrix(X)
  # Z_mem <- bigmemory::as.big.matrix(Z)
  #
  # subject_list <- split(seq_along(subject), subject)

  repeat {
    if (a>maxiter||nrm<0.0005) {break}
    cat("iter: ", a, "\n")

    #It is the parallel version of ith_ans = lapply(1:n, function(i) {} )
    ith_thread_fn <- function(core_i) {
      positions <- parallel::splitIndices(n, cores)[[core_i]]

      ubfi_sum <- array(0, c(p,p))
      ubsi_sum <- array(0, c(p,1))
      sigma_sum <- 0
      D_sum <- array(0, c(q,q))

      for (i in positions) {
        # cat("pos: ", i, "\n")
        subject_rows <- subject == i
        X_i = subset(X, subject_rows)
        Z_i = subset(Z, subject_rows)
        y_i = subset(y, subject_rows)

        U_i <- MASS::ginv(
          Dinv + (t(Z_i) %*% Rinv %*% Z_i)
        )
        W_i <- Rinv - (Rinv %*% Z_i %*% U_i %*% t(Z_i) %*% Rinv)

        ubfi_sum <- ubfi_sum + t(X_i) %*% W_i %*% X_i
        ubsi_sum <- ubsi_sum + t(X_i) %*% W_i %*% y_i

        #Update u (random effect)
        b_i <- U_i %*% t(Z_i) %*% Rinv %*% (y_i - X_i %*% beta)
        sigma_sum <- sigma_sum + as.numeric(
          t(y_i - X_i %*% beta) %*% W_i %*% (y_i - X_i %*% beta)
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
      answers <- plyr::llply(1:cores, ith_thread_fn, .parallel = TRUE)
      # answers <- parallel::clusterApply(c2, 1:cores, ith_thread_fn)
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

    final.D=round(final.D,10)
    final.D[!matrixcalc::is.positive.definite(final.D)] = D    #If the matrix is not positive definite, use the initial
    D = final.D

    # if (median(svd(final.D)$d)<10) {
    #   D = final.D
    # } else {
    #   D = D
    # }

    final.sigma[final.sigma < 0 ] = 1
    sigma = final.sigma

    # if (final.sigma > 0) {
    #   sigma = final.sigma
    # } else {
    #   sigma = sigma
    # }
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
