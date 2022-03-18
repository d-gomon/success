#' @title Assist users in parameter selection
#'
#' @description  This function can be used to determine some of the vital
#' parameters used to construct control charts in this package.
#'
#' @details Depending on the specified arguments, the function will return
#' parameters. If \code{covariate_names} is not specified, the returned
#' risk-adjustment models will be trivial. If \code{formula} is not specified
#' but \code{covariate_names} are,
#' the function assumes the simplest form for the regression model
#' (cov1 + cov2 + ...). If \code{followup} is not specified, no \code{glmmod}
#' will be determined
#'
#'
# #' @param covariate_names A vector containing the covariate names to be used
# #' in risk-adjustment procedures. (character)
#' @param formula A formula with right-hand side (RHS) indicating the form in which
#' the covariates should be used for the Cox and GLM
#' regression models. LHS of the formula will be ignored, and can be left empty.
#' @param baseline_data A \code{data.frame} for determining a baseline performance
#' metric. Rows should represent subjects and the
#' following named columns should be present: \describe{
#'   \item{\code{entrytime}:}{time of entry into study (numeric);}
#'   \item{\code{survtime}:}{time from entry until event (numeric);}
#'   \item{\code{censorid}:}{censoring indicator (0 = right censored, 1 = observed),
#'    (integer).}
#' } and optionally additional covariates used for risk-adjustment.
#' @param data A \code{data.frame} with data on which the user wants to construct
#' quality control charts.
#' Rows should represent subjects and the
#' following named columns should be present: \describe{
#'   \item{\code{entrytime}:}{time of entry into study (numeric);}
#'   \item{\code{survtime}:}{time from entry until event (numeric);}
#'   \item{\code{censorid}:}{censoring indicator (0 = right censored, 1 = observed),
#'    (integer).}
#' } and optionally additional covariates used for risk-adjustment.
#' @param followup (optional): The value of the follow-up time to be used to determine event time.
#' Event time will be equal to \code{entrytime + followup} for each subject.
#' @param theta The value of the expected log-hazard ratio. In other words: the logarithm of the
#' expected increase in the odds of failure/failure rate. Default is log(2) (detecting a
#' doubling of the odds/failure rate).
#' @param time Timeframe over which the type I error of the control chart should be
#' limited. Should be in the same unit as \code{survtime} in \code{data}. If left
#' unspecified, the maximum entrytime in \code{baseline_data} is taken. (numeric)
#' @param alpha Required maximal type I error (between 0 and 1) of the procedure
#' over the timeframe specified in \code{time}. Default is 0.05. (numeric)
#'
#'
#'
#'
#' @return A list of parameters to feed to quality control charts in this
#' package
#' \itemize{
#' \item{call: }{The call used to obtain output.}
#' \item{data: }{The data used in the call to the function.}
#' \item{baseline_data: }{The baseline_data used in the call to the function}
#' \item{glmmod: }{A \code{\link[stats:glm]{glm()}} model which can be fed to
#' the \code{\link[success:funnel_plot]{funnel_plot()}}
#'  and \code{\link[success:bernoulli_cusum]{bernoulli_cusum()}} functions.}
#'  \item{coxphmod: }{A \code{\link[survival:coxph]{coxph()}} model which can be
#'  fed to the \code{\link[success:cgr_cusum]{cgr_cusum()}} and
#'  \code{\link[success:cgr_cusum]{cgr_cusum()}} functions.}
#'  \item{psi: }{Estimated Poisson arrival rate in \code{data}.}
#' }
#'
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @importFrom stats glm
#' @importFrom stats formula
#' @importFrom stats update
#' @importFrom stats family
#' @importFrom stats as.formula
#' @importFrom stats binomial
#'
#' @export
#'
#'
#' @author Daniel Gomon
#'
#'
#' @examples
#' require(survival)
#'
#' #Minimal example - no risk-adjustment
#' pars_min <- parameter_assist(baseline_data = surgerydat,
#' data = subset(surgerydat, unit == 1))
#'
#' #Specifying all parameters
#' pars <- parameter_assist(baseline_data = surgerydat,
#' data = subset(surgerydat, unit == 1),
#' formula = formula("survtime ~ age + sex + BMI"), followup = 100)






parameter_assist <- function(baseline_data, data, formula,
                             followup, theta = log(2), time, alpha = 0.05){
  call <- match.call()
  p0 <- NULL

  message("Checking provided data.")
  data <- check_data(data)
  message('Checking provided baseline_data.')
  baseline_data <- check_data(baseline_data)
  covariate_names <- NULL



  #Covariate checks
  if(!missing(formula)){
    covariate_names <- labels(terms(formula))
  }

  if(!is.null(covariate_names)){
    if(!all(covariate_names %in% colnames(data))){
      stop("Specified covariates not (all) present in data.")
    }
    if(!all(covariate_names %in% colnames(baseline_data))){
      stop("Specified covariates not (all) present in baseline_data.")
    }
  }
  if("unit" %in% colnames(data)){
    data <- data[,c("survtime", "entrytime", "censorid", "unit", covariate_names)]
  } else{
    data <- data[,c("survtime", "entrytime", "censorid", covariate_names)]
  }
  if("unit" %in% colnames(baseline_data)){
    baseline_data <- baseline_data[,c("survtime", "entrytime", "censorid", "unit",
                                      covariate_names)]
  } else{
    baseline_data <- baseline_data[,c("survtime", "entrytime", "censorid",
                                      covariate_names)]
  }



  #Attempt glmmod construction (only if followup is specified)
  if(!missing(followup)){
    #Don't construct glmmod if missing followup
    if(!missing(formula)){
      formulaglm <- update(formula, (survtime <= followup) & (censorid == 1) ~ .)
    } else if (!is.null(covariate_names)){
      formulaglm <- as.formula(paste0("(survtime <= followup) & (censorid == 1) ~ ",
                                   paste(covariate_names, collapse = "+")))
    } else {
      formulaglm <- as.formula("(survtime <= followup) & (censorid == 1) ~ 1")
    }
    environment(formulaglm) = environment()

    glmmod <- glm(formula = formulaglm, family = binomial(link = "logit"),
                  data = baseline_data)
    p0 <- length(which((baseline_data$survtime <= followup) &
                         (baseline_data$censorid == 1)))/nrow(baseline_data)
  } else{
    message("glmmod will not be determined: missing argument followup.")
    glmmod <- NULL
  }

  if(!missing(formula)){
    formulasurv <- update(formula, survival::Surv(survtime, censorid) ~ .)
  } else if(!is.null(covariate_names)){
    formulasurv <- as.formula(paste0("Surv(survtime, censorid) ~ ",
                                     paste(covariate_names, collapse = "+")))
  } else{
    formulasurv <- as.formula("Surv(survtime, censorid) ~ 1")
  }
  coxphmod <- coxph(formula = formulasurv, data = baseline_data, model = TRUE)



  #Determine arrival rate psi of specified data frame
  psi <- nrow(data)/max(data$entrytime)

  if(missing(time)){
    time <- max(baseline_data$entrytime)
  }


  #Create list to return
  parlist <- list(call = call,
                  data = data,
                  baseline_data = baseline_data,
                  glmmod = glmmod,
                  coxphmod = coxphmod,
                  theta = theta,
                  psi = psi,
                  time = time,
                  alpha = alpha)
  if(!missing(followup)){
    parlist$followup <- followup
  }
  if(!is.null(p0)){
    parlist$p0 <- p0
  }
  class(parlist) <- "assisted_success"

  return(parlist)
}
