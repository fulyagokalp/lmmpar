# https://github.com/fulyagokalp/lmmpar # nolint
# devtools::install_github("fulyagokalp/lmmpar") # nolint


#' Parallel Linear Mixed Model
#'
#' Embarrassingly Parallel Linear Mixed Model calculations spread across local cores which repeat until convergence.  All calculations are currently done locally, but theoretically, the calculations could be extended to multiple machines.
#'
#' @param Y matrix of responses with observations/subjects on column and repeats for each observation/subject on rows. It is (m x n) dimensional.
#' @param X observed design matrices for fixed effects. It is (m*n x p) dimensional.
#' @param Z observed design matrices for random effects. It is (m*n x q) dimensional.
#' @param subject vector of positions that belong to each subject.
#' @param beta fixed effect estimation vector with length p.
#' @param R variance-covariance matrix of residuals.
#' @param D variance-covariance matrix of random effects.
#' @param sigma initial sigma value.
#' @param maxiter the maximum number of iterations that should be calculated.
#' @param cores the number of cores. Why not to use maximum?!
#' @param verbose boolean that defaults to print iteration context
#' @importFrom MASS ginv
#' @importFrom stats median
#' @export
#' @examples
#'
#' # Set up fake data
#' n <- 1000  # number of subjects
#' m <- 4      # number of repeats
#' N <- n * m  # true size of data
#' p <- 15     # number of betas
#' q <- 2      # width of random effects
#'
#' # Initial parameters
#' # beta has a 1 for the first value.  all other values are ~N(10, 1)
#' beta <- rbind(1, matrix(rnorm(p, 10), p, 1))
#' R <- diag(m)
#' D <- matrix(c(16, 0, 0, 0.025), nrow = q)
#' sigma <- 1
#'
#' # Set up data
#' subject <- rep(1:n, each = m)
#' repeats <- rep(1:m, n)
#'
#' subj_x <- lapply(1:n, function(i) cbind(1, matrix(rnorm(m * p), nrow = m)))
#' X <- do.call(rbind, subj_x)
#' Z <- X[, 1:q]
#' subj_beta <- lapply(1:n, function(i) mnormt::rmnorm(1, rep(0, q), D))
#' subj_err <- lapply(1:n, function(i) mnormt::rmnorm(1, rep(0, m), sigma * R))
#'
#' # create a known response
#' subj_y <- lapply(
#'   seq_len(n),
#'   function(i) {
#'     (subj_x[[i]] %*% beta) +
#'       (subj_x[[i]][, 1:q] %*% subj_beta[[i]]) +
#'       subj_err[[i]]
#'   }
#' )
#' Y <- do.call(rbind, subj_y)
#'
#' # run the algorithm in parallel to recover the known betas
#' ans <- lmmpar(
#'   Y,
#'   X,
#'   Z,
#'   subject,
#'   beta = beta,
#'   R = R,
#'   D = D,
#'   cores = 2,
#'   sigma = sigma,
#'   verbose = TRUE
#' )
#' str(ans)
# Reference paper: Schafer, J.L.S., 1998, some improved procedures for linear mixed models
#Ui, Wi and beta is calculated in parallel form. Then,sigma and D are calculated with final beta
#Function is updated for stacked vector and matrices.
lmmpar <- function(
  Y, X, Z, subject,
  beta, R, D, sigma,
  maxiter = 500,
  cores = 8,
  verbose = TRUE
){

  iter <- 0
  p <- nrow(beta)
  n <- length(unique(subject))
  m <- nrow(Y) / n
  N <- n * m
  q <- dim(D)[1]

  Dinv <- ginv(D)
  Rinv <- ginv(R)

  cores <- floor(cores)
  if (cores > 1) {
    doParallel::registerDoParallel(cores)
  }

  Y_mem <- bigmemory::as.big.matrix(Y)
  X_mem <- bigmemory::as.big.matrix(X)
  Z_mem <- bigmemory::as.big.matrix(Z)

  subject_list <- split(seq_along(subject), subject)

  verbose <- isTRUE(verbose)
  nrm <- 1

  repeat {
    if (iter > maxiter || nrm < 0.0005) break

    if (verbose) cat("iter: ", iter, "\n")

    #It is the parallel version of ith_ans = lapply(1:n, function(i) {} )
    core_fn <- function(core_i) {
      positions <- parallel::splitIndices(n, cores)[[core_i]]

      ubfi_sum <- array(0, c(p, p))
      ubsi_sum <- array(0, c(p, 1))
      sigma_sum <- 0
      D_sum <- array(0, c(q, q))

      for (i in positions) {
        subject_rows <- subject_list[[i]]
        Y_i <- Y_mem[subject_rows, , drop = FALSE] # nolint
        X_i <- X_mem[subject_rows, , drop = FALSE] # nolint
        Z_i <- Z_mem[subject_rows, , drop = FALSE] # nolint

        U_i <- MASS::ginv(
          Dinv + (t(Z_i) %*% Rinv %*% Z_i)
        )
        W_i <- Rinv - (Rinv %*% Z_i %*% U_i %*% t(Z_i) %*% Rinv)

        ubfi_sum <- ubfi_sum + t(X_i) %*% W_i %*% X_i
        ubsi_sum <- ubsi_sum + t(X_i) %*% W_i %*% Y_i

        #Update u (random effect)
        b_i <- U_i %*% t(Z_i) %*% Rinv %*% (Y_i - X_i %*% beta)
        sigma_sum <- sigma_sum + as.numeric(
          t(Y_i - X_i %*% beta) %*% W_i %*% (Y_i - X_i %*% beta)
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
      answers <- lapply(1, core_fn)
    } else {
      answers <- plyr::llply(1:cores, core_fn, .parallel = TRUE)
    }

    ubfi_total <- array(0, c(p, p)) + Reduce("+", lapply(answers, `[[`, "ubfi_sum"))
    ubsi_total <- array(0, c(p, 1)) + Reduce("+", lapply(answers, `[[`, "ubsi_sum"))
    sigma_total <- 0 + Reduce("+", lapply(answers, `[[`, "sigma_sum"))
    D_total <- array(0, c(q, q)) + Reduce("+", lapply(answers, `[[`, "D_sum"))

    #Final calculations
    final.beta <- ginv(ubfi_total) %*% ubsi_total

    final.D <- (1 / n) * as.matrix(D_total)
    final.sigma <- as.numeric(
      (1 / N) * as.matrix(sigma_total)[1, 1]
    )

    nrm <- norm(final.beta - beta)

    beta <- final.beta

    final.D <- round(final.D, 10)
    #If the matrix is not positive definite, use the previous parameter
    final.D[!matrixcalc::is.positive.definite(final.D)] <- D
    D <- final.D

    ## nolint start
    # if (median(svd(final.D)$d)<10) {
    #   D = final.D
    # } else {
    #   D = D
    # }
    ## nolint end

    if (any(final.sigma < 0)) {
      final.sigma[final.sigma < 0] <- sigma[final.sigma < 0]
    }
    sigma <- final.sigma

    Dinv <- ginv(D)

    iter <- iter + 1
  }

  return(list(
    iterations = iter,
    beta = beta,
    D = D,
    sigma = sigma,
    norm = nrm
  ))

}
