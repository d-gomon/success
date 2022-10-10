#' Weibull hazard, cumulative hazard and inverse cumulative hazard
#'
#' @description Functions which return the hazard, cumulative
#' hazard and inverse cumulative hazard at time t for a Weibull distribution
#' with shape parameter \eqn{\lambda}{\lambda}, scale parameter \eqn{\theta}{\theta}
#'  and true hazard ratio \eqn{\mu}{\mu}.
#'
#' @details The hazard function of a Weibull distribution is given by:
#' \deqn{h(t| \lambda, \theta, \mu) = \frac{\lambda}{\theta} \left(\frac{t}{\theta}  \right)^{\lambda -1} e^\mu}{h(t|\lambda, \theta, \mu) = \lambda/\theta (t/\theta)^(\lambda -1) exp(\mu)}
#' The cumulative hazard (with true hazard ratio \eqn{\mu}{\mu}) is given by:
#' \deqn{H(t| \lambda, \theta, \mu) = \left( \frac{t}{\theta} \right)^{\lambda} e^\mu}{H(\lambda, \theta, \mu) = ( t/\theta )^{\lambda} exp(\mu)}
#' The inverse cumulative hazard (with true hazard ratio \eqn{\mu}{\mu}) by:
#' \deqn{H^{-1}(t| \lambda, \theta, \mu) = \theta \left( \frac{t}{e^\mu}  \right)^{1/\lambda}}{H^(-1)(t| \lambda, \theta, \mu) = \theta (t/e^\mu)^{1/\lambda}}
#'
#' @return Value of specified function at time \eqn{t}{t}.
#'
#'
#' @param t time of evaluation.
#' @param lambda parameter of the exponential distribution.
#' @param mu (optional) true excess hazard rate \eqn{\mu}{\mu}.
#' @name weib_hazards
NULL
#> NULL


#' @rdname weib_hazards
#' @export
haz_weib <- function(t, lambda, theta, mu = log(1)){
  exp(mu) * (lambda/theta) * (t/theta)^(lambda -1)
}

#' @rdname weib_hazards
#' @export
chaz_weib <- function(t, lambda, theta, mu = log(1)){
  lambda * t * exp(mu)
}

#' @rdname weib_hazards
#' @export
inv_chaz_weib <- function(t, lambda, theta, mu = log(1)){
  theta * (t/exp(mu))^(1/lambda)
}
