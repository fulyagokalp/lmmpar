context("single_example")


test_that("single example", {


  n <- 10000
  m <- 4
  N <- n*m
  p <- 4
  q <- 2
  # c <- 2  #contamination factor
  # pb <- 0
  # pe <- 0
  nu <- 3 #degrees of freedom of t-distribution

  #Initial parameters
  beta = matrix(rmnorm(p+1, 10, 1),p+1,1)
  R = diag(m)
  D = matrix(c(16, 0, 0, 0.025), nrow=q)
  sigma = 1

  beta_start = beta
  R_start = R
  D_start = D
  sigma_start = sigma

  #Empty arrays for x and z; matrices for b, e, and y
  #
  #   X <- matrix(NA, (m*n), p+1)
  #   Z <- matrix(NA, (m*n), q)
  #   b <- matrix(NA, n, q)
  #   e <- matrix(NA, n, m)
  #   y <- matrix(NA, m*n, 1)
  subject <- rep(1:n, each = m)
  repeats <- rep(1:m, n)


  # X <- array(NA,c(m,p+1,n))
  # Z <- array(NA,c(m,q,n))
  # b <- matrix(NA, q, n)
  # e <- matrix(NA, m, n)
  # y <- matrix(NA, m, n)

  # myresultX <- list()
  # myresultb <- list()
  # myresulte <- list()

  myresultX <- lapply(1:n, function(i) cbind(1,matrix(rnorm(m*p),nrow=m)))
  X <- do.call(rbind, myresultX)
  Z <- X[,1:q]
  myresultb <- lapply(1:n, function(i) rmnorm(1, rep(0, q), D))
  myresulte <- lapply(1:n, function(i) rmnorm(1, rep(0, m), sigma*R))

  myresulty <- lapply(
    1:n,
    function(i) {
      (myresultX[[i]] %*% beta) +
        (myresultX[[i]][,1:q] %*% myresultb[[i]]) +
        myresulte[[i]]
    }
  )
  y <- do.call(rbind, myresulty)

  # for (i in 1:n) {
  #   X[,,i] <- cbind(1,matrix(rnorm(m*p),nrow=m))
  #   Z[,,i] <- X[,1:q,i,drop=FALSE]
  #   b[,i] <- rmnorm(1, rep(0, q), D)
  #   e[,i] <- rmnorm(1, rep(0, m), sigma*R)
  #   y[,i] <- X[,,i]%*%beta+Z[,,i]%*% b[,i]+e[,i]
  # }

  #Z <- matrix(c(1,1,1,1,8,10,12,14),m,p)
  #X <- matrix(c(1,1,1,1,8,10,12,14),m,p)
  #
  #   X2 <- rep(c(8,10,12,14) ,n)
  #   Z2 <- rep(c(8,10,12,14) ,n)
  #   subject <- rep(c(1:n),each=m)

  #VarBi <- (1+(c^2-1)*pb)*D
  #VarEi <- (1+(c^2-1)*pe)*sigma*R
  #V <- Z%*%D%*%t(Z)+sigma*R
  #Dinv <- ginv(D)
  #Rinv <- ginv(R)
  #Vinv <- ginv(V)

  #y <- matrix(0,m,n)
  #u <- array(0,dim=c(p,1,n))

  #for (i in 1:n){
  # b <- rmnorm(1, rep(0, 2), VarBi)
  # e <- rmnorm(1, rep(0, 4), VarEi)

  # set.seed(i)
  # y[,i] <- X%*%beta+Z%*%b+e
  # }

  cat("\nstarting timings...\n\n")


  cores_vals <- c(4, 1)

  timings <- list()
  pars <- list()
  for (i in seq_along(cores_vals)) {
    timing <- system.time({
      # profvis::profvis({
      ans.gauss <- lmm.ep.em(
        y,
        X,
        Z,
        subject,
        beta = beta,
        R = R,
        D = D,
        cores = cores_vals[i],
        sigma = sigma
      )
      print(ans.gauss)
    }) %>% print()

    timing <- as.list(timing)
    timing$cores <- cores_vals[i]
    timing$n <- n
    timing$repeats <- m
    timing$p <- p
    timings[[i]] <- timing
    pars$Beta_start = beta_start
    pars$R_start = R_start
    pars$D_start = D_start
    pars$sigma_start = sigma_start
    #expect_true(TRUE)
  }

  cat("\n")
  print(as.data.frame(do.call(rbind, timings)))
  print(pars)

  #   Y <- as.vector(y)
  #   m1 <- lme4::lmer(Y~X2+(1|subject)+(Z2-1|subject),REML=FALSE)
  #   AICGauss <- 2*61-2*logLik(m1,REML=FALSE)
  #   summary(m1)$coefficients[,1][1]
  #   summary(m1)$coefficients[,1][2]
  #   unlist(summary(m1))$varcor.subject
  #   unlist(summary(m1))$varcor.subject.1
  #   sigma(m1)


  # expect_true(abs(AICGauss - ans$AIC) < 10)

  # expect_equal(length(ans), 7)

})


m1 <- lme4::lmer(Y~X[,2]+X[,3]+X[,4]+X[,5]+X[,6]+X[,7]+X[,8]+X[,9]+X[,10]+X[,11]+X[,12]+X[,13]+X[,14]+X[,15]+X[,16]+X[,17]+X[,18]+X[,19]+X[,20]+X[,21]+X[,22]+X[,23]+X[,24]+X[,25]+X[,26]+X[,27]+X[,28]+X[,29]+X[,30]+X[,31]+X[,32]+X[,33]+X[,34]+X[,35]+X[,36]+X[,37]+X[,38]+X[,39]+X[,40]+X[,41]+X[,42]+X[,43]+X[,44]+X[,45]+X[,46]+X[,47]+X[,48]+X[,49]+X[,50]+X[,51]+(1|subject)+(Z[,2]-1|subject),REML=FALSE)
