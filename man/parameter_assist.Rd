% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_assist.R
\name{parameter_assist}
\alias{parameter_assist}
\title{Assist users in parameter selection}
\usage{
parameter_assist(baseline_data, data, formula, followup, theta = log(2),
  time, alpha = 0.05, maxtheta = log(6))
}
\arguments{
\item{baseline_data}{A \code{data.frame} for determining a baseline performance
metric. Rows should represent subjects and the
following named columns should be present: \describe{
\item{\code{entrytime}:}{time of entry into study (numeric);}
\item{\code{survtime}:}{time from entry until event (numeric);}
\item{\code{censorid}:}{censoring indicator (0 = right censored, 1 = observed),
(integer).}
} and optionally additional covariates used for risk-adjustment.}

\item{data}{A \code{data.frame} with data on which the user wants to construct
quality control charts.
Rows should represent subjects and the
following named columns should be present: \describe{
\item{\code{entrytime}:}{time of entry into study (numeric);}
\item{\code{survtime}:}{time from entry until event (numeric);}
\item{\code{censorid}:}{censoring indicator (0 = right censored, 1 = observed),
(integer).}
} and optionally additional covariates used for risk-adjustment.}

\item{formula}{A formula with right-hand side (RHS) indicating the form in which
the covariates should be used for the Cox and GLM
regression models. LHS of the formula will be ignored, and can be left empty.}

\item{followup}{(optional): The value of the follow-up time to be used to determine event time.
Event time will be equal to \code{entrytime + followup} for each subject.}

\item{theta}{The value of the expected log-hazard/odds ratio. In other words: the logarithm of the
expected increase in the odds/hazard ratio. Default is log(2) (detecting a
doubling of the odds/failure rate).}

\item{time}{Timeframe over which the type I error of the control chart should be
limited. Should be in the same unit as \code{survtime} in \code{data}. If left
unspecified, the maximum entrytime in \code{baseline_data} is taken. (numeric)}

\item{alpha}{Required maximal type I error (between 0 and 1) of the procedure
over the timeframe specified in \code{time}. Default is 0.05. (numeric)}

\item{maxtheta}{Maximum value the maximum likelihood estimate for
parameter \eqn{\theta}{\theta} can take. If \code{detection = "lower"}, \code{-abs(theta)}
will be the minimum value the maximum likelihood estimate for
parameter \eqn{\theta}{\theta} can take.  Default is \code{log(6)}, meaning that
at most a 6 times increase/decrease in the odds/hazard ratio is expected.}
}
\value{
A list of parameters to feed to quality control charts in this
package:
\itemize{
\item call: The call used to obtain output.
\item data: The data used in the call to the function.
\item baseline_data: The baseline_data used in the call to the function
\item glmmod: A \code{\link[stats:glm]{glm()}} model which can be fed to
the \code{\link[success:funnel_plot]{funnel_plot()}}
and \code{\link[success:bernoulli_cusum]{bernoulli_cusum()}} functions.
\item coxphmod: A \code{\link[survival:coxph]{coxph()}} model which can be
fed to the \code{\link[success:cgr_cusum]{cgr_cusum()}} and
\code{\link[success:cgr_cusum]{cgr_cusum()}} functions.
\item theta: Expected increase in the odds/hazard ratio.
\item psi: Estimated Poisson arrival rate in \code{data}.
\item time: Time frame over which to restrict type I error.
\item alpha: Desired level of type I error for control limit determination.
\item maxtheta: Maximum expected increase/decrease in the odds/hazard ratio.
}
}
\description{
This function can be used to determine some of the vital
parameters used to construct control charts in this package.
}
\details{
Depending on the specified arguments, the function will return
parameters. If \code{covariate_names} is not specified, the returned
risk-adjustment models will be trivial. If \code{formula} is not specified
but \code{covariate_names} are,
the function assumes the simplest form for the regression model
(cov1 + cov2 + ...). If \code{followup} is not specified, no \code{glmmod}
will be determined
}
\examples{
require(survival)

#Minimal example - no risk-adjustment
pars_min <- parameter_assist(baseline_data = surgerydat,
data = subset(surgerydat, unit == 1))

#Specifying all parameters
pars <- parameter_assist(baseline_data = surgerydat,
data = subset(surgerydat, unit == 1),
formula = formula("survtime ~ age + sex + BMI"), followup = 100)
}
\author{
Daniel Gomon
}
