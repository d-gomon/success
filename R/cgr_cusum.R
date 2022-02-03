#' @title Continuous time Generalized Rapid response CUSUM (CGR-CUSUM)
#'
#' @description This function performs the CGR-CUSUM procedure
#' described in Gomon et al. (in review). For detection purposes, it suffices
#' to determine the value of the chart at the times of failure. This can be
#'  achieved by leaving \code{ctimes} unspecified.
#' The function has two vital parameters, at least one of which must be specified:
#' \itemize{
#' \item{\code{coxphmod}: }{Cox proportional hazards model to be used for
#' risk-adjustment. If \code{cbaseh} is not specified, it will be determined
#' from \code{coxphmod} numerically.}
#' \item{\code{cbaseh}: }{The cumulative baseline hazard rate to use for chart
#' construction. If specified with \code{coxphmod} missing, no risk-adjustment
#' will be performed}
#' }
#'
#' @details The CGR-CUSUM can be used to test for a change of unknown positive fixed size \eqn{\theta}{\theta}
#'  in the subject-specific hazard rate from \eqn{h_i(t)}{h_i(t)} to \eqn{h_i(t) e^\theta}{h_i(t) exp(\theta)}
#'  starting from some unknown subject \eqn{\nu}{\nu}. The starting time of the first subject
#'  who had an increase in failure rate as well as the estimated increase in the
#'  hazard rate are shown in the output.
#'  The CGR-CUSUM is determined as
#' \deqn{\max_{1 \leq \nu \leq n} \left( \hat{\theta}_{\geq \nu}(t) N_{\geq \nu}(t) - \left( \exp\left( \hat{\theta}_{\geq \nu}(t) \right) - 1 \right) \Lambda_{\geq \nu}(t)\right),}{max{1<=\nu<=n} (\theta_{>=\nu}(t)N_{>=\nu}(t)) - (exp(\theta_{>=\nu}(t))-1) \Lambda_{>=\nu}(t)),}
#' where  \deqn{N(\geq \nu)(t) = \sum_{i \geq \nu} N_i(t),}{N_{>=\nu}(t) = \sum_{i>=\nu} N_i(t),}
#' with \eqn{N_i(t)}{N_i(t)} the counting process for the failure at time \eqn{t}{t} of subject \eqn{i}{i}
#' and \deqn{\Lambda_{\geq \nu}(t) = \sum_{i \geq \nu} \Lambda_i(t),}{\Lambda_{>=\nu}(t) = \sum_{i>=\nu}\Lambda_i(t),}
#' where \eqn{\Lambda_i(t)}{\Lambda_i(t)} is the cumulative intensity of subject \eqn{i}{i} at time \eqn{t}{t}.
#'
#' @param data A \code{data.frame} with rows representing subjects and the
#' following named columns: \describe{
#'   \item{\code{entrytime}:}{time of entry into study (numeric);}
#'   \item{\code{survtime}:}{time from entry until event (numeric);}
#'   \item{\code{censorid}:}{censoring indicator (0 = right censored, 1 = observed),
#'    (integer).}
#' } and optionally additional covariates used for risk-adjustment.
#' @param coxphmod A Cox proportional hazards regression model as
#' produced by
#' the function \code{\link[survival:coxph]{coxph()}}. Suggested: \cr
#' \code{coxph(Surv(survtime, censorid) ~ covariates, data = data)}. \cr
#' Alternatively, a list with the following elements:
#' \describe{
#'   \item{\code{formula}:}{a \code{\link[stats:formula]{formula()}} in the form \code{~ covariates};}
#'   \item{\code{coefficients}:}{a named vector specifying risk adjustment coefficients
#'   for covariates. Names must be the same as in \code{formula} and colnames of \code{data}.}
#' }
#' @param cbaseh A function that returns the unadjusted cumulative
#' baseline hazard \eqn{H_0(t)}{H_0(t)}. If \code{cbaseh} is missing but
#' \code{coxphmod} has been
#' specified as a survival object, this baseline hazard rate will be determined
#' using the provided \code{coxphmod}.
#' @param ctimes (optional): Vector of construction times at which the value of the chart should be
#' determined. When not specified, the chart is constructed at all failure times.
#' @param h (optional): Value of the control limit. The chart will only be
#' constructed until the value of the control limit has been reached or
#' surpassed.
#' @param stoptime (optional): Time after which the value of the chart should no
#' longer be determined. Default = max(failure time). Useful when ctimes
#' has not been specified.
#' @param C (optional): A numeric value indicating how long after entering the study
#' patients should no longer influence the value of the chart. This is
#' equivalent to right-censoring every observation at time \code{entrytime} + C.
# @param twosided (optional): A boolean indicating whether a two-sided CUSUM
#  should be constructed.
# If \code{TRUE}, 2 CUSUM charts are determined. One to check for an increase
# of \eqn{e^\theta}{exp(\theta)} and one for
# a decrease of \eqn{e^{-\theta}}{exp(-\theta)} in the hazard rate w.r.t.
#  the baseline hazard. Default is \code{FALSE}.
#' @param pb (optional): A boolean indicating whether a progress bar should
#' be shown. Default is \code{FALSE}.
#' @param ncores number of cores to use to parallelize the computation of the
#' CGR-CUSUM chart. If ncores = 1 (default), no parallelization is done. You
#' can use \code{\link[parallel:detectCores]{detectCores()}} to check how many
#' cores are available on your computer.
#' @param cmethod Method to calculate chart values. One of the following:
#' \itemize{
#' \item \code{"memory"} (default): matrix formulation of the problem
#' (faster for high volume/long time construction)
#' \item \code{"CPU"} calculates the value of the CGR-CUSUM for every
#' time point from scratch. Recommended for small data volume
#' (lower initialization time).
#' }
#' @param dependencies (optional): When \code{ncores > 1}, specify a list of
#' variables/functions/other dependencies to be exported to the core clusters
#' for parallel computation.
#' @param detection Should an \code{"upper"} or \code{"lower"} CGR-CUSUM be
#' constructed. Upper CUSUMs can be used to monitor for an increase in the
#' failure rate, while lower CUSUMs can be used to monitor for a decrease in the
#' failure rate.
#'
#' @return An object of class "cgrcusum" containing:
#' \itemize{
#' \item \code{CGR}: a \code{data.frame} with named columns:
#' \describe{
#'   \item{\code{time}:}{times at which chart is constructed;}
#'   \item{\code{value}:}{value of the chart at corresponding times;}
#'   \item{\code{exp_theta_t}:}{value of MLE \eqn{e^{\theta_t}}{e^(\theta_t)};}
#'   \item{\code{S_nu}}{time from which patients are considered for constructing the chart.}
#' }
#' \item \code{call}: the call used to obtain output;
#' \item \code{stopind}: indicator for whether the chart was stopped by
#' the control limit;
#' \item \code{h}: Specified value for the control limit.
#' }
# There are \code{\link[cgrcusum:plot.cgrcusum]{plot}} and
#  \code{\link[cgrcusum:runlength.cgrcusum]{runlength}} methods for "cgrcusum" objects.
#'
#' @importFrom stats loess
#' @importFrom stats predict
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#'
#' @author Daniel Gomon
#' @family quality control charts
#' @seealso \code{\link[cgrcusum]{plot.cgrcusum}}, \code{\link[cgrcusum]{runlength.cgrcusum}}
#'
#'
#' @examples
#' require(survival)
#' #Select only the data of the first hospital in the surgerydat data set
#' tdat <- subset(surgerydat, hosp_num == 1)
#'
#' #We know that the cumulative baseline hazard in the data set is
#' #Exponential(0.01). If you don't know the cumulative baseline, we suggest
#' #leaving the cbaseh argument empty and determining a coxphmod (see help)
#' tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
#'
#' #Determine a risk-adjustment model using a Cox proportional hazards model.
#' #Outcome (survival) is regressed on the available covariates:
#' exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
#' tcoxmod <- coxph(exprfit, data= surgerydat)
#'
#' #Determine the values of the chart
#' cgr <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE)
#' #Alternatively, cbaseh can be left empty when specifying coxphmod through coxph()
#' #cgr <- cgr_cusum(data = tdat, coxphmod = tcoxmod, pb = TRUE)
#'
#' #plot the CGR-CUSUM
#' plot(cgr)






cgr_cusum <- function(data, coxphmod, cbaseh, ctimes, h, stoptime,
                     C, pb = FALSE, ncores = 1, cmethod = "memory",
                     dependencies, detection = "upper"){
  call = match.call()
  #-------------------------------DATA CHECKS-----------------------------#
  #Basic data checks (global for BK, CGR and Bernoulli)
  if(missing(data)){
    stop("Please provide data to construct chart.")
  } else{
    data <- check_data(data)
  }
#  compriskcheck <- "cause" %in% colnames(data)
#  if(compriskcheck){
#    message("Competing risks specified.")
#  }
  # determine chronological failure times
  data$otime <- data$entrytime + data$survtime
  if(!missing(C)){
    tempidx <- which(data$otime < data$entrytime + C)
    data[tempidx,]$otime <- data$entrytime + C
    data[tempidx,]$censorid <- rep(length(tempidx), 0)
  }
  # determine the default construction times (all failtimes), if none specified
  if(missing(ctimes)){
    ctimes <- unique(data$otime)
  }
  if(missing(stoptime)){
    stoptime <- max(data$otime[is.finite(data$otime)])
  }
  checkcbase <- FALSE
  if(missing(coxphmod) | is.null(coxphmod)){
    coxphmod <- NULL
  } else if(inherits(coxphmod, "coxph") & missing(cbaseh)){
      checkcbase <- TRUE
      message("Missing cumulative baseline hazard. Determining using provided Cox PH model.")
      cbaseh <- extract_hazard(coxphmod)$cbaseh
      #Old Loess method:
      # cbaselo <- loess(cbase_temp$hazard~cbase_temp$time)
      # cbaseh <- function(x) predict(cbaselo, x)
  } else if(is.list(coxphmod)){
        if(all(c("formula", "coefficients") %in% names(coxphmod))){
          checkcoxlist <- TRUE
        } else{
          stop("coxphmod does not contain $formula and/or $coefficients.")
        }
  } else if(is.null(coxphmod)){
    coxphmod <- NULL
  } else{ stop("coxphmod is not a list or survival object.")}
  if(missing(cbaseh)){
    if(!checkcbase){
      stop("Please specify cbaseh (function) or coxphmod as Survival object.")
    }
  }
  if(!missing(h)){
    if(length(h) > 1){
      stop("Please specify only 1 control limit.")
    }
  }

  #----------------------------FUNCTION BODY--------------------------#

  #Only determine value of the chart at at relevant times
  ctimes <- sort(ctimes[which(ctimes <= stoptime)])
  #When determining CGR-CUSUM at only 1 time point, cmethod = "memory" doesnt work
  if(length(ctimes) == 1){
    cmethod = "CPU"
  }


  #If method = CPU, iteratively calculate the value at each time
  if(cmethod == "CPU"){
    Gt <- matrix(0, nrow = 0, ncol = 4)
    startval <- 0
    stopind <- FALSE
    #progress bar
    if(pb){
      pb2 <- txtProgressBar(min= 0, max = length(ctimes), style = 3)
    }
    #Calculate value of chart at each construction time separately
    for(j in seq_along(ctimes)){
      if(pb){
        setTxtProgressBar(pb2, value = j)
      }
      #Calculate value of CGR at time ctimes[j]
      temcgr <- cgr_helper(data = data, ctime = ctimes[j],
                           coxphmod = coxphmod, cbaseh = cbaseh)
      Gt <- rbind(Gt, c(ctimes[j], temcgr$val, exp(temcgr$theta), temcgr$starttime))
      if (!missing(h)){if(temcgr$val >= h) {stopind = TRUE; break}}
    }
    colnames(Gt) <- c("time", "value", "exp_theta_t", "S_nu")
    if(pb){
      close(pb2)
    }
    cgr <- list(CGR = Gt,
                stopind = stopind,
                call = call)
  } else if(cmethod == "memory" && detection == "upper"){ #Matrix form + no control limit + possibly multi core
    #Calculate using matrix formulation, without having to stop when control limit is surpassed
    #Can perform multi-core computations
    Gt <- cgr_helper_mat(data = data, ctimes = ctimes, h=h, coxphmod = coxphmod, cbaseh = cbaseh,
                           ncores = ncores, displaypb = pb, dependencies = dependencies)
    cgr <- list(CGR = Gt,
                call = call)
  } else if(cmethod == "memory" && detection == "lower"){
    Gt <- cgr_helper_mat_down(data = data, ctimes = ctimes, h=h, coxphmod = coxphmod, cbaseh = cbaseh,
                              ncores = ncores, displaypb = pb, dependencies = dependencies)
    cgr <- list(CGR = Gt,
                call = call)

  }



  if(!missing(h)){
    if(detection == "upper"){
      cgr$h <- abs(h)
    } else{
      cgr$h <- -abs(h)
    }
  }
  class(cgr) <- "cgrcusum"
  cgr
}


