context("single_example")


test_that("single example", {

  is_personal <- !is.na(file.info(".lintr")$size)

  if (!is_personal) {
    expect_true(TRUE)
  } else {

    n <- 1000
    m <- 4
    N <- n*m
    p <- 50
    q <- 2
    nu <- 3 #degrees of freedom of t-distribution

    #Initial parameters
    beta = rbind(1,matrix(rmnorm(p, 10, 1),p,1))
    #beta = rbind(1,matrix(rbinom(p, 1, 0.1),p,1))
    R = diag(m)
    D = matrix(c(16, 0, 0, 0.025), nrow=q)
    sigma = 1

    beta_start = beta
    R_start = R
    D_start = D
    sigma_start = sigma

    # set up data
    subject <- rep(1:n, each = m)
    repeats <- rep(1:m, n)

    myresultX <- lapply(1:n, function(i) cbind(1,matrix(rnorm(m*p),nrow=m)))
    X <- do.call(rbind, myresultX)
    Z <- X[,1:q]
    myresultb <- lapply(1:n, function(i) rmnorm(1, rep(0, q), D))
    myresulte <- lapply(1:n, function(i) rmnorm(1, rep(0, m), sigma*R))

    # found data
    myresulty <- lapply(
      seq_len(n),
      function(i) {
        (myresultX[[i]] %*% beta) +
          (myresultX[[i]][,1:q] %*% myresultb[[i]]) +
          myresulte[[i]]
      }
    )
    y <- do.call(rbind, myresulty)

    cat("\nstarting timings...\n\n")
    cores_vals <- c(4, 1)

    timings <- list()
    pars <- list()
    for (i in seq_along(cores_vals)) {
      timing <- system.time({
        ans.gauss <- lmm_ep_em(
          y,
          X,
          Z,
          subject,
          beta = beta,
          R = R,
          D = D,
          cores = cores_vals[i],
          sigma = sigma,
          verbose = TRUE
        )
        print(ans.gauss)
      })

      print(timing)
      timing <- as.list(timing)
      timing$cores <- cores_vals[i]
      timing$n <- n
      timing$repeats <- m
      timing$p <- p
      timings[[i]] <- timing

      pars$Beta_start <- beta_start
      pars$R_start <- R_start
      pars$D_start <- D_start
      pars$sigma_start <- sigma_start

      expect_true(TRUE)
    }

    cat("\n")
    print(as.data.frame(do.call(rbind, timings)))
    print(pars)
  }

})


if (FALSE) {
  # nolint_start
  # Y <- as.vector(y)
  # m1 <- lme4::lmer(Y~X2+(1|subject)+(Z2-1|subject),REML=FALSE)
  # AICGauss <- 2*61-2*logLik(m1,REML=FALSE)
  # summary(m1)$coefficients[,1][1]
  # summary(m1)$coefficients[,1][2]
  # unlist(summary(m1))$varcor.subject
  # unlist(summary(m1))$varcor.subject.1
  # sigma(m1)
  #
  #
  # expect_true(abs(AICGauss - ans$AIC) < 10)
  #
  # expect_equal(length(ans), 7)
  # nolint_end
}
