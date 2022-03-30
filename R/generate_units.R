#'
#' @title Generate units with specified failure rate
#'
#' @description Generate \code{n_sim} units with subjects arriving according to
#' a Poisson
#' process with rate \code{psi} until \code{time}. Failure rate is determined
#' either from \code{cbaseh} or \code{inv_cbaseh}, or from specified
#' \code{coxphmod}. Covariates will be resampled from \code{baseline_data} if
#' specified.
#'
#' @details In a Poisson arrival process, inter-arrival times are exponentially
#' distributed with parameter \code{psi}. If \code{cbaseh} is specified,
#' the inverse baseline hazard will be determined using \code{\link[stats:uniroot]{uniroot()}}.
#' The times of failure are then determined using \code{\link[success:gen_surv_times]{gen_surv_times()}}.
#'
#'
#' @param time A numeric value indicating until what time from the start of the
#' study subjects can arrive.
#' @param psi Poisson arrival rate for subjects.
#' @param n_sim An integer indicating how many units should be generated.
#' Default is 20.
#' @param cbaseh A function returning the cumulative baseline hazard
#' at each relevant time point. Will be numerically inverted to generate
#' failure times. If
#' the inverse cumulative baseline hazard function is available, please specify
#' \code{inv_cbaseh} instead.
#' @param inv_cbaseh A function returning the inverse cumulative baseline
#' hazard at each relevant time point.
#' @param coxphmod (optional): A cox proportional hazards model generated using
#' \code{\link[survival:coxph]{coxph()}} or a list containing:
#' \describe{
#'   \item{\code{formula}:}{a \code{\link[stats:formula]{formula()}} in the form \code{~ covariates};}
#'   \item{\code{coefficients}:}{a named vector specifying risk adjustment coefficients
#'   for covariates. Names must be the same as in \code{formula} and colnames of \code{data}.}
#' } If both \code{cbaseh} and \code{inv_cbaseh} are missing, the hazard rate
#' will be determined from \code{coxphmod}.
#' @param baseline_data (optional): A \code{data.frame} used for covariate resampling
#' with rows representing subjects and columns containing covariates to use for risk-adjustment.
#' Should only be specified in combination with \code{coxphmod}.
#' @param interval (optional) A numeric vector of length 2 indicating in which
#' range of values the failure times of subjects should be determined. By default,
#' failure times will be restricted between 0 and 9e12.
#' @param mu (optional) The increased log hazard ratio at the generated units with
#' respect to the specified baseline hazard rate. Default is log(1) = 0.
#'
#'
#'
#'
#' @return A \code{data.frame} with rows representing subjects and the
#' following named columns:
#' \describe{
#'   \item{\code{entrytime}:}{time of subject entry into the study;}
#'   \item{\code{survtime}:}{survival time of subject;}
#'   \item{\code{censorid}:}{censoring indicator, 0 = censored, 1 = observed;}
#'   \item{\code{unit}:}{unit number;}
#'   \item{\code{expmu}:}{exponent of the log hazard ratio used to generate
#'   survival times;}
#'   \item{\code{psival}:}{arrival rate at unit;}
#'   \item{\code{covariates}:}{covariates resampled from \code{baseline_data}.}
#' }
#'
#'
#' @author Daniel Gomon
#' @export
#' @importFrom stats uniroot
#'
#' @examples
#' require(survival)
#' #Fit a Cox model
#' exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
#' tcoxmod <- coxph(exprfit, data= surgerydat)
#'
#' #Generate 30 hospitals with on average 2 patients per day arriving
#' #according to the Cox model determined above, with resampling from the
#' #original data set. The hazard rate at the hospitals is twice the baseline
#' #hazard.
#' generate_units(time = 50, psi = 2, n_sim = 30, coxphmod = tcoxmod,
#' baseline_data = surgerydat, mu = log(2))
#'
#'
#'


generate_units <- function(time, psi, n_sim = 20, cbaseh, inv_cbaseh,
                           coxphmod = NULL, baseline_data, interval = c(0, 9e12),
                           mu = 0){

  manualcbaseh <- FALSE

  #Determine amount of required significant digits
  if(!missing(baseline_data)){
    signif_dig <- tryCatch({max(sapply(baseline_data$entrytime, function(x) match(TRUE, round(x, 1:20) == x)))-1},
                           error = function(cond){2})
  } else{
    signif_dig <- 2
  }

  #--------------Step 1: Data initialization-------------

  #Generate n_sim instances
  #We create a data frame to contain all information
  entrytime_temp <- numeric(0)
  unit_temp <- numeric(0)
  for(i in 1:n_sim){
    arrivtimes <- round(gen_arriv_times(psi = psi, t = time), digits = signif_dig)
    entrytime_temp <- c(entrytime_temp, arrivtimes)
    unit_temp <- c(unit_temp, rep(i, length(arrivtimes)))
  }

  #Resample covariates from existing data if supplied
  if(!missing(baseline_data)){
    df_temp <- baseline_data[sample(nrow(baseline_data), length(entrytime_temp), replace = TRUE),]
    df_temp$entrytime <- entrytime_temp
    df_temp$unit <- unit_temp
  } else if(missing(baseline_data)){
    df_temp <- data.frame("entrytime" = entrytime_temp, "unit" = unit_temp)
  }

  df_temp$survtime <- rep(NA, nrow(df_temp))
  #df_temp is now the dataframe containing all entrytimes and covariates for subjects

  #------------------------Step 2: inv_cbaseh determination ----------------

  if(!missing(inv_cbaseh)){
    #Do nothing
  } else if(!missing(cbaseh)){
    #Determine the inverse cumulative baseline hazard H^{-1}_0(t)
    inv_cbaseh_temp <- function(y, lower = interval[1], upper = interval[2]){
      tryCatch({return(unname(unlist(uniroot((function (x) cbaseh(x) - y),
                                             lower = lower, upper = upper)$root)))},
               error = function(cond){ return(upper)})
    }

    #inv_cbaseh_temp2 <- function(y, lower = 0, upper = 9E12){
    #  uniroot((function (x) cbaseh(x) - y), lower = lower, upper = upper)$root
    #}
    inv_cbaseh <- Vectorize(inv_cbaseh_temp)
  } else if(!missing(coxphmod)){
    manualcbaseh <- TRUE
    haz_temp <- extract_hazard(coxphmod)
    inv_cbaseh <- haz_temp$inv_cbaseh
    max_time <- haz_temp$max_time
  }


  #-------------------Step 3: survival times determination-----------------

  #Determine survival times
  if(missing(baseline_data)){
    survtimes <- gen_surv_times(invchaz = inv_cbaseh, mu = mu, data = nrow(df_temp), coxphmod = coxphmod)
  } else{
    survtimes <- gen_surv_times(invchaz = inv_cbaseh, mu = mu, data = df_temp, coxphmod = coxphmod)
  }

  #Add them to the data frame.
  df_temp$survtime <- survtimes
  #If we determined cbaseh manually, censor all observations which were generated at maximum failure time
  if(manualcbaseh == TRUE){
    df_temp$censorid <- ifelse(survtimes == max(max_time), 0, 1)
  } else{
    #Add trivial censoring indicator
    df_temp$censorid <- rep(1, nrow(df_temp))
  }
  return(df_temp)
}
