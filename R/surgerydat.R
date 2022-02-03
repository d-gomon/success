#' @title Simulated data set with data of surgery procedures
#' performed at multiple hospitals.
#'
#' @description Data about patients and their surgery procedure from 30 simulated hospitals
#' with patient arrivals in the first 400 days after the start of the study. \cr
#' Patient survival times were determined using a risk-adjusted Cox proportional hazards model
#' with coefficients age = 0.003, BMI = 0.02 and sexmale = 0.2 and exponential baseline hazard rate
#' \eqn{h_0(t, \lambda = 0.01) e^\mu}{h_0(t, \lambda = 0.01) exp(\mu)}.
#' Some hospitals have an increased failure rate:
#' \itemize{
#' \item Hospitals 1-15: \eqn{e^\mu = 1}{exp(\mu) = 1}
#' \item Hospitals 16-30: \eqn{e^\mu = 2}{exp(\mu) = 2}
#' } which means that the hazard rate at hospitals 16-30 is twice higher than exponential(\eqn{\lambda}). \cr
#' The arrival rate \eqn{\psi}{\psi} of patients at a hospital differs. The arrival rates are:
#' \itemize{
#' \item Hospitals 1-5 & 16-20: 0.5 patients per day (small hospitals)
#' \item Hospitals 6-10 & 21-25: 1 patient per day (medium sized hospitals)
#' \item Hospitals 11-15 & 26-30: 1.5 patients per day (large hospitals)
#' } These are then respectively small, medium and large hospitals.
#'
#'
#'
#' @format A \code{data.frame} with 12010 rows and 9 variables:
#' \describe{
#'   \item{entrytime}{Time of entry of patient into study (numeric)}
#'   \item{survtime}{Time from entry until failure of patient (numeric)}
#'   \item{censorid}{Censoring indicator (0 - right censored, 1 - observed) (integer)}
#'   \item{hosp_num}{Hospital number at which patient received treatment (integer)}
#'   \item{expmu}{True excess hazard used for generating patient survival (numeric)}
#'   \item{psival}{Poisson arrival rate at hospital which the patient was at (numeric)}
#'   \item{age}{Age of the patient (numeric)}
#'   \item{sex}{Sex of the patient (factor)}
#'   \item{BMI}{Body mass index of the patient (numeric)}
#' }
"surgerydat"
