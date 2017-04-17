context("single_example")

# library(doParallel)
# registerDoParallel(cores = 8)

test_that("single example", {


  n <- 100000
  m <- 4
  N <- n*m
  p <- 2
  q <- 2
  c <- 2  #contamination factor
  pb <- 0
  pe <- 0
  nu <- 3 #degrees of freedom of t-distribution

  #Initial parameters
  #beta = matrix(c(1,2,1),p+1,1)  #It should have p+1 values
  beta = matrix(c(10,0,1),p+1,1)
  R = diag(m)
  D = matrix(c(16, 0, 0, 0.025), nrow=q)
  sigma=1

  #Empty arrays for x and z; matrices for b, e, and y
  X <- array(NA,c(m,p+1,n))
  Z <- array(NA,c(m,q,n))
  b <- matrix(NA, q, n)
  e <- matrix(NA, m, n)
  y <- matrix(NA, m, n)


  for (i in 1:n) {
    X[,,i] <- cbind(1,matrix(rnorm(m*p),nrow=m))
    Z[,,i] <- X[,1:q,i,drop=FALSE]
    b[,i] <- rmnorm(1, rep(0, q), D)
    e[,i] <- rmnorm(1, rep(0, m), sigma*R)
    y[,i] <- X[,,i]%*%beta+Z[,,i]%*% b[,i]+e[,i]
  }


  #Z <- matrix(c(1,1,1,1,8,10,12,14),m,p)
  #X <- matrix(c(1,1,1,1,8,10,12,14),m,p)

  X2 <- rep(c(8,10,12,14) ,n)
  Z2 <- rep(c(8,10,12,14) ,n)
  subject <- rep(c(1:n),each=m)

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


  # system.time(ans <- lmm.t.em.par(
  #   y,
  #   X,
  #   Z,
  #   beta = beta,
  #   R = R,
  #   D = D,
  #   first_parallel = TRUE,
  #   sigma = sigma
  # )
  # )

  # system.time(ans2 <- lmm.t.em.par2(
  #   y,
  #   X,
  #   Z,
  #   beta = beta,
  #   R = R,
  #   D = D,
  #   first_parallel = TRUE,
  #   second_parallel = TRUE,
  #   sigma = sigma
  # )
  # )

  first_second <- expand.grid(first = c(FALSE, TRUE), second = c(FALSE, TRUE))
  # browser()

  timings <- list()
  for (i in seq_len(nrow(first_second))) {
    first <- first_second$first[i]
    second <- first_second$second[i]
    timing <- system.time({
    # profvis::profvis({
      ans.gauss <- lmm.ep.em(
        y,
        X,
        Z,
        beta = beta,
        R = R,
        D = D,
        first_parallel = first,
        second_parallel = second,
        sigma = sigma
      )
      print(ans.gauss)
    }) %>% print()

    timing <- as.list(timing)
    timing$first <- first
    timing$second <- second
    timings[[i]] <- timing

    expect_true(TRUE)
  }

  cat("\n")
  print(as.data.frame(do.call(rbind, timings)))

  print(pryr::object_size(y))
  print(pryr::object_size(X))
  print(pryr::object_size(Z))
  print(pryr::object_size(R))
  print(pryr::object_size(D))

  # browser()



  Y <- as.vector(y)
  m1 <- lme4::lmer(Y~X2+(1|subject)+(Z2-1|subject),REML=FALSE)
  AICGauss <- 2*61-2*logLik(m1,REML=FALSE)
  summary(m1)$coefficients[,1][1]
  summary(m1)$coefficients[,1][2]
  unlist(summary(m1))$varcor.subject
  unlist(summary(m1))$varcor.subject.1
  sigma(m1)


  # expect_true(abs(AICGauss - ans$AIC) < 10)

  # expect_equal(length(ans), 7)

})
