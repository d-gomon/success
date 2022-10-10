#' Exponential hazard, cumulative hazard and inverse cumulative hazard
#'
#' @description Functions which return the hazard, cumulative
#' hazard and inverse cumulative hazard at time t for an exponential distribution
#' with parameter \eqn{\lambda}{\lambda} and true hazard ratio \eqn{\mu}{\mu}.
#'
#' @details The hazard function of an exponential distribution is given by:
#' \deqn{h(t| \lambda, \mu) = \lambda e^\mu}{h(t|\lambda, \mu) = \lambda e^\mu}
#' The cumulative hazard (with true hazard ratio \eqn{\mu}{\mu}) is given by:
#' \deqn{H(t| \lambda, \mu) = \lambda t e^\mu}{H(t|\lambda, \mu) = \lambda * t * exp(\mu)}
#' The inverse cumulative hazard (with true hazard ratio \eqn{\mu}{\mu}) by:
#' \deqn{H^{-1}(t| \lambda, \mu) = \frac{t}{\lambda e^\mu}}{H^(-1)(t|\lambda, \mu) = t/(\lambda exp(\mu))}
#'
#' @return Value of specified function at time \eqn{t}{t}.
#'
#'
#' @param t time of evaluation.
#' @param lambda parameter of the exponential distribution.
#' @param mu (optional) true excess hazard rate \eqn{\mu}{\mu}.
#' @name exp_hazards
NULL
#> NULL


#' @rdname exp_hazards
#' @export
haz_exp <- function(t, lambda, mu = log(1)){
  lambda * exp(mu)
}

#' @rdname exp_hazards
#' @export
chaz_exp <- function(t, lambda, mu = log(1)){
  return(lambda * t * exp(mu))
}

#' @rdname exp_hazards
#' @export
inv_chaz_exp <- function(t, lambda, mu = log(1)){
  return(t/(lambda * exp(mu)))
}
