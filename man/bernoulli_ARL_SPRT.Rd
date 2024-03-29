% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bernoulli_ARL.R
\name{bernoulli_ARL_SPRT}
\alias{bernoulli_ARL_SPRT}
\title{Average run length for Bernoulli CUSUM using Integral Equation methodology}
\usage{
bernoulli_ARL_SPRT(h, n_grid, Wncdf, glmmod, theta, theta_true, p0,
  tol = 1e-06)
}
\arguments{
\item{h}{Control limit for the Bernoulli CUSUM}

\item{n_grid}{Number of state spaces used to discretize the outcome space (when \code{method = "MC"})
or number of grid points used for trapezoidal integration (when \code{method = "SPRT"}).
Increasing this number improves accuracy, but can also significantly increase computation time.}

\item{Wncdf}{A function returning the values of the (risk-adjusted) cumulative
distribution function (cdf) for the singletons Wn.}

\item{glmmod}{Generalized linear regression model used for risk-adjustment as produced by
the function \code{\link[stats:glm]{glm()}}. Suggested: \cr
\code{glm(as.formula("(survtime <= followup) & (censorid == 1) ~ covariates"), data = data)}. \cr
Alternatively, a list containing the following elements:
\describe{
\item{\code{formula}:}{a \code{\link[stats:formula]{formula()}} in the form \code{~ covariates};}
\item{\code{coefficients}:}{a named vector specifying risk adjustment coefficients
for covariates. Names must be the same as in \code{formula} and colnames of \code{data}.}
}}

\item{theta}{The \eqn{\theta}{\theta} value used to specify the odds ratio
\eqn{e^\theta}{e^\theta} under the alternative hypothesis.
If \eqn{\theta >= 0}{\theta >= 0}, the average run length for the upper one-sided
Bernoulli CUSUM will be determined. If \eqn{\theta < 0}{\theta < 0},
the average run length for the lower one-sided CUSUM will be determined.
Note that \deqn{p_1 = \frac{p_0 e^\theta}{1-p_0 +p_0 e^\theta}.}{p1 = (p0 * e^\theta)/(1-p0+p0 * e^\theta).}}

\item{theta_true}{The true log odds ratio \eqn{\theta}{\theta}, describing the
true increase in failure rate from the null-hypothesis. Default = log(1), indicating
no increase in failure rate.}

\item{p0}{The baseline failure probability at \code{entrytime + followup} for individuals.}
}
\description{
Internal function that calculates the ARL using the connection
between the ARL of a Wald SPRT and a CUSUM.
}
\keyword{internal}
