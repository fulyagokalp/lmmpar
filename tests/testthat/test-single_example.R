context("single_example")

test_that("single example", {




  n <- 5000
  m <- 4
  N <- n*m
  p <- 2
  q <- 2
  c <- 2  #contamination factor
  pb <- 0.10
  pe <- 0
  Z <- matrix(c(1,1,1,1,8,10,12,14),m,p)
  X <- matrix(c(1,1,1,1,8,10,12,14),m,p)

  X2 <- rep(c(8,10,12,14) ,n)
  Z2 <- rep(c(8,10,12,14) ,n)
  subject <- rep(c(1:n),each=m)

  beta = matrix(c(17,0.8),p,1)
  R = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),m,m)
  D = matrix(c(4,0,0,0.0225),q,q)

  VarBi <- (1+(c^2-1)*pb)*D
  VarEi <- (1+(c^2-1)*pe)*R
  V <- Z%*%D%*%t(Z)+R
  Dinv <- ginv(D)
  Rinv <- ginv(R)
  Vinv <- ginv(V)

  y <- matrix(0,m,n)
  u <- array(0,dim=c(p,1,n))

  for (i in 1:n){
    b <- rmnorm(1, rep(0, 2), VarBi)
    e <- rmnorm(1, rep(0, 4), VarEi)

    #set.seed(i)
    y[,i] <- X%*%beta+Z%*%b+e
  }




  ans <- lmm.ep.em(
    y,
    X,
    Z,
    beta = beta,
    R = R,
    D = D
  )


  Y <- as.vector(y)
  m1 <- lme4::lmer(Y~X2+(1|subject)+(Z2-1|subject),REML=FALSE)
  AICGauss <- 2*61-2*logLik(m1,REML=FALSE)
  summary(m1)$coefficients[,1][1]
  summary(m1)$coefficients[,1][2]
  unlist(summary(m1))$varcor.subject
  unlist(summary(m1))$varcor.subject.1
  sigma(m1)


  # expect_true(abs(AICGauss - ans$AIC) < 10)

  expect_equal(length(ans), 7)

})
