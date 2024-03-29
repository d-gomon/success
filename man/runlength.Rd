% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runlength.R
\name{runlength}
\alias{runlength}
\alias{runlength.cgrcusum}
\alias{runlength.bkcusum}
\alias{runlength.bercusum}
\title{Determine run length of a CUSUM chart}
\usage{
runlength(chart, h)

\method{runlength}{cgrcusum}(chart, h, ...)

\method{runlength}{bkcusum}(chart, h, ...)

\method{runlength}{bercusum}(chart, h, ...)
}
\arguments{
\item{chart}{A \code{cgrcusum}, \code{bkcusum} or \code{bercusum} chart.}

\item{h}{Control limit h to be used when determining the run length}

\item{...}{Other parameters}
}
\value{
The run length of the chart with the given control limit.
}
\description{
This function can be used to calculate the run length of a 'cgrcusum', 'bkcusum'
or 'bercusum'
chart when using control limit h
}
\section{Methods (by class)}{
\itemize{
\item \code{runlength(cgrcusum)}: determines runlength of \code{cgrcusum} object

\item \code{runlength(bkcusum)}: determines runlength of \code{bkcusum} object

\item \code{runlength(bercusum)}: determines runlength of \code{bercusum} object

}}
\examples{
exprfitber <- as.formula("(survtime <= 100) & (censorid == 1) ~ age + sex + BMI")
glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 14), glmmod = glmmodber,
                   followup = 100, theta = log(2))
#Determine the run length of the above Bernoulli CUSUM when using a control limit
#of h = 1.
runlength(bercus, h = 1)
}
\author{
Daniel Gomon
}
