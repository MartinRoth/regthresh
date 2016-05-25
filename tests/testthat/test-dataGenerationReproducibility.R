context("Data Generation")

library(doParallel)
registerDoParallel(1)

test_that("correlated distorted weibull generation works", {
  set.seed(123)
  xpar <- list()
  xpar$S <- 3L
  xpar$T <- 10L
  tau     <- 0.95
  epsilon <- 0.25
  gamma   <- 0.5
  xi      <- 0.15
  kappasV <- rbeta(xpar$S, 2, 5) + 0.5
  betas   <- runif(xpar$S, 2, 4)

  params <- setMarginalParameters(tau = tau, epsilon = epsilon, betas = betas,
                                  kappas = kappasV, gamma = gamma, xi = xi,
                                  xpar = xpar)

  expect_equal_to_reference(params, "./marginalParameters.rds")

  copulaObj <- copula::normalCopula(0, xpar$S)

  expect_equal_to_reference(copulaObj, "./copulaObject.rds")

  simData <- generateDWData(params, copulaObj, xpar)

  expect_equal_to_reference(simData, "./simulatedData.rds")
})
