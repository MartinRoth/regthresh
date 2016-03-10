## Set parameters
hazDWParUpdate <- function(parameters) {
    # set threshold dependent on tau
    # and compute also where the GPD is pure
    parameters$loc <- hazDWSetThreshold(parameters)
    parameters$loc2 <- parameters$loc + parameters$epsilon
    return(parameters)
}

hazDWSetThreshold <- function(parameters) {
    # help function for hazDWParUpdate
    tau   <- parameters$tau
    beta  <- parameters$beta
    kappa <- parameters$kappa
    return((-log(1- tau) * beta^kappa)^(1/kappa))
}

## transition function
linearTransition <- function(x) {
    y <- rep(1, length(x))
    isBetweenZeroAndOne <- (x > 0 & x < 1)
    y[x >= 1] <- 0
    y[isBetweenZeroAndOne] <- 1 - x[isBetweenZeroAndOne]
    return(y)
}

polynomialTransition <- function(x) {
    y <- rep(1, length(x))
    isBetweenZeroAndOne <- (x > 0 & x < 1)
    y[x >= 1] <- 0
    y[isBetweenZeroAndOne] <- 2 * x[isBetweenZeroAndOne]^3 - 3 * x[isBetweenZeroAndOne]^2 + 1
    return(y)
}

#smoothTransition <- function(x) {
#    # the scaling has to be done before
#    a <- smoothTransitionHelpFunction(x)
#    b <- smoothTransitionHelpFunction(1 - x)
#    return(1 - a/(a + b))
#}
#
#smoothTransitionHelpFunction <- function(x) {
#    y <- rep(0, length(x))
#    index <- x > 0
#    y[index] <- exp(-1/x[index])
#    return(y)
#}

hazDW <- function(x, parameters) {
    a <- hazW(x, parameters)
    b <- hazGpd(x, parameters)
    c <- polynomialTransition((x-parameters$loc)/parameters$epsilon)
    return(c * a + (1 - c) * b)
}

hazW <- function(x, parameters) {
    kappa <- parameters$kappa
    beta  <- parameters$beta
    return(kappa/beta * (x / beta)^(kappa - 1))
    #return(hazardRateWeibull(x, parameters))
}

hazGpd <- function(x, parameters) {
    xi    <- parameters$xi
    gamma <- parameters$gamma
    loc   <- parameters$loc

    isGreaterLoc  <- x >= loc

    y               <- rep(0, length(x))
    y[isGreaterLoc] <- (gamma * loc + xi * (x[isGreaterLoc] - loc))^(-1)
    return(y)
}

## Cumulative hazard rate
HazDW <- function(x, parameters) {
    xi    <- parameters$xi
    gamma <- parameters$gamma
    loc   <- parameters$loc
    loc2  <- parameters$loc2

    y <- rep(NA, length(x))
    isSmaller <- x <= loc
    isGreater <- x >+ loc2
    isBetween <- !(isSmaller | isGreater)

    y[isSmaller] <- HazW(x[isSmaller], parameters)
    y[isBetween] <- HazW(loc, parameters) + HazTransitionWGpd(x[isBetween], parameters)
    y[isGreater] <- HazW(loc, parameters) + HazTransitionWGpd(loc2, parameters) +
                    HazGpd(x[isGreater], parameters) - HazGpd(loc2, parameters)
    return(y)
}

HazW <- function(x, parameters) {
    kappa <- parameters$kappa
    beta  <- parameters$beta
    return((x / beta)^kappa)
}

HazGpd <- function(x, parameters) {
    xi    <- parameters$xi
    gamma <- parameters$gamma
    loc   <- parameters$loc
    scale <- gamma * loc

    y <- 1/xi * log(1 + xi * (x - loc)/scale)
    return(y)
}

HazTransitionWGpd <- function(x, parameters) {
    a <- HazTransitionW(x, parameters)
    b <- HazTransitionGpd(x, parameters)
    return(a + b)
}

HazTransitionW <- function(x, parameters) {
    y <- HazTransitionWHelp(x, parameters) - HazTransitionWHelp(parameters$loc, parameters)
    return(y)
}

HazTransitionGpd <- function(x, parameters) {
    y <- HazTransitionGpdHelp(x - parameters$loc, parameters) - HazTransitionGpdHelp(0, parameters)
    return(y)
}

HazTransitionWHelp <- function(x, parameters) {
    kappa   <- parameters$kappa
    beta    <- parameters$beta
    loc     <- parameters$loc
    epsilon <- parameters$epsilon

    b <- loc / epsilon

    c <- 2 * kappa / (epsilon^3 * (kappa + 3)) * x^(kappa + 3)
    d <- -3 * kappa * (2 * b + 1) / (epsilon^2 * (kappa+2)) * x^(kappa + 2)
    e <- 6 * loc * kappa * (b + 1) / (epsilon^2 * (kappa + 1)) * x^(kappa + 1)
    f <- -(b^2 * (2 * b + 3) - 1) * x^kappa

    y <- (c + d + e + f) / beta^kappa

    return(y)
}

HazTransitionGpdHelp <- function(x, parameters) {
    gamma   <- parameters$gamma
    xi      <- parameters$xi
    loc     <- parameters$loc
    epsilon <- parameters$epsilon
    sigma   <- gamma * loc

    stopifnot(xi > -sigma/epsilon)

    a <- (epsilon * xi)^(-1)
    X <- xi * x + sigma

    b <- -2/3 * a * X^3
    c <- (3 * a * sigma + 3/2) * X^2
    d <- -6 * sigma * (a * sigma + 1) * X
    e <- sigma^2 * ( 2 * a * sigma + 3) * log(X)

    y <- a^2 / xi * (b + c + d + e)
    return(y)
}

#' Get marginal parameters of the distorted Weibull distribution
#' @param tau Quantile where GPD behaviors sets in
#' @param epsilon Length of transition interval
#' @param betas Weibull scale parameter
#' @param kappas Weibull shape parameter
#' @param gamma GPD dispersion coefficient
#' @param xi GPD shape parameter
#' @export
setMarginalParameters <- function(tau, epsilon, betas, kappas, gamma, xi) {
  mParams <- list()
  for (site in 1 : xpar$S) {
    mParams[[site]] <- list(tau = tau, epsilon = epsilon
                            ,beta = betas[site], kappa = kappas[site]
                            ,gamma = gamma, xi = xi)
  }
  mParams <- lapply(mParams, hazDWParUpdate)
  return(mParams)
}

#' Calculates the Weibull shape parameter
#'
#' that corresponds to a GPD shape parameter xi at quantile tau
#' @inheritParams setMarginalParameters
#' @export
calculateCorrespondingWeibullShape <- function(xi, tau) {
  kappa <- (xi * (- log(1 - tau)) + 1)^(-1)
  return(kappa)
}
calculatePenUltimateShapeParameter <- function(kappa, tau) {
  xi <- (1 - kappa)/kappa * (-log(1 - tau))^(-1)
  return(xi)
}

generateDWData <- function(marginalParameters, copula, xpar) {
  uniform <- rCopula(xpar$T, copula)
  data    <- foreach(s = 1 : xpar$S, .combine = "cbind") %dopar% {
    return(qDW(uniform[, s], parameters = marginalParameters[[s]]))
  }
  return(data)
}

getTrueReturnLevels <- function(params, xpar, returnPeriods, blockLength) {
  tmp <- foreach(i = 1 : xpar$S, .combine = "cbind") %do% {
    loc    = params[[i]]$loc
    scale  = params[[i]]$gamma * loc
    shape  = params[[i]]$xi
    lambda = params[[i]]$tau
    p      = 1 - 1/returnPeriods/blockLength
    return(qgpd(p = p, loc = loc, scale = scale, shape = shape
                , lambda = lambda))
  }
  colnames(tmp) = 1 : xpar$S
  rownames(tmp) = retPeriods
  return(tmp)
}

## Distorted Weibull
rDW <- function(n, parameters) {
    f <- function(x, parameters, tau) {
        return(1 - exp(-HazDW(x, parameters)) - tau)
    }
    unif <- runif(n)
    result <- numeric(n)
    for (i in 1 : n) {
       result[i] <- uniroot(f, c(0, 1000), parameters, tau = unif[i])$root
    }
    return(result)
}

qDW <- function(tau, parameters) {
    f <- function(x, parameters, tau) {
        return(1 - exp(-HazDW(x, parameters)) - tau)
    }
    result <- numeric()
    for (i in 1 : length(tau)) {
       result[i] <- uniroot(f, c(0, 1000), parameters, tau = tau[i])$root
    }
    return(result)
}

pDW <- function(q, parameters) {
    return(1 - exp(-HazDW(q, parameters)))
}

dDW <- function(x, parameters) {
    return(hazDW(x, parameters) * exp(-HazDW(x, parameters)))
}


