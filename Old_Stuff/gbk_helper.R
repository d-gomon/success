#' Generalized BK-CUSUM (GBK-CUSUM) helper -
#' single time point
#'
#' @description This function calculates the value of the CGR-CUSUM at one
#' specified timepoint
#'
#'
#'
#' @param data data frame containing the following named columns:
#' \code{entrytime} (numeric - time of entry into study), \code{survtime}
#' (numeric - time from entry until event), \code{censorid} (integer - censoring
#' indicator: 0 - right censored, 1 - observed), \code{cause} (factor - cause of
#' event - competing risks).
#' @param ctime construction time (single) at which the value of the chart
#' should be determined.
#' @param coxphmod a cox proportional hazards regression model as produced by
#' the function \code{\link[survival:coxph]{coxph}}. Obtained using:
#' \code{coxph(Surv(survtime, censorid) ~ covariates, data = data)}.
#' Alternatively, a list with $formula (~ covariates)
#' and $coefficients (named vector specifying risk adjustment coefficients
#' for covariates - names must be the same as in $formula and colnames).
#' @param cbaseh a function which returns the non risk-adjusted cumulative
#' baseline hazard \eqn{h_0(t)}. If \code{cbaseh} is missing but
#' \code{coxphmod} has been
#' specified as a survival object, this baseline hazard rate will be determined
#' using the provided \code{coxphmod}.
#' @param displaypb (optional) boolean indicating whether a progress bar should be
#' displayed
#'
#' @return A list containing the following:
#' \itemize{
#'   \item $val value of CGR-CUSUM at specified time point
#'   \item $theta value at corresponding time of the MLE \eqn{\hat{\theta}_t}
#'   \item $starttime time from which individuals contribute to the chart \eqn{S_\nu}
#' }
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @noRd
#' @keywords internal
#'
#' @author Daniel Gomon
#'
#' @seealso \code{\link{bkcusum}},
#' \code{\link{bercusum}} (step 2)
#'
#' @examples
#' #TO-DO

gbk_helper <- function(data, ctime, coxphmod, cbaseh){
  entrytime <- NULL
  #check whether a construction time has been specified, otherwise take max
  if(missing(ctime)){
    ctime <- max(data$otime)
  }
  if(!"censorid" %in% colnames(data)){
    cat("No censoring mechanism specified. Assuming data is uncensored.")
    data$censorid <- rep(1, nrow(data))
  }
  #Calculate subject specific failure risk (exp(\beta Z))
  riskdat <- calc_risk(data, coxphmod)

  #Consider only patients with starting time larger than helperstimes[i]
  tdat <- subset(data, entrytime <= ctime)
  tdatidx <- which(data$entrytime <= ctime)
  #Determine amount of (relevant) failures in this subset
  NDT <- length(which(tdat$censorid == 1 & tdat$otime <= ctime))
  failnow <- length(which(tdat$censorid == 1 & tdat$otime == ctime))
  NDT_low <- NDT - failnow
  #Contribution of active subjects to A(t)
  activecbaseh <- cbaseh(ifelse(tdat$otime < ctime, tdat$otime, ctime)-tdat$entrytime)
  #Error check for empty contribution
  activecbaseh[is.na(activecbaseh)] <- 0
  #Determine A(t) (have to risk-adjust)
  AT <- sum(riskdat[tdatidx] * activecbaseh)
  #Determine \hat{\theta}_t
  thetat <- log(NDT/AT)
  thetat_low <- log(NDT_low/AT)
  if (is.infinite(thetat)){
    thetat <- 0
  }
  if(is.infinite(thetat_low)){
    thetat_low <- 0
  }
  #Determine value of CGI-CUSUM using only patients with S_i > helperstimes[i]
  val_down <- thetat_low * NDT_low - (exp(thetat_low) - 1) * AT
  val_up <- thetat* NDT - (exp(thetat)- 1) * AT
  #return list of relevant values
  return(list(val_down = val_down,
              val_up = val_up,
              theta = thetat))
}




