#' @title Continuous time BK-CUSUM
#'
#' @description This function implements the BK-CUSUM procedure based on the
#' Biswas & Kalbfleisch (2008) CUSUM. To construct the Biswas & Kalbfleisch
#' (2008) CUSUM, set \code{C = 1} (years) or \code{C = 365} (days).
#'  For detection purposes, it is sufficient
#' to determine the value of the chart at the times of failure. This can be
#'  achieved by leaving \code{ctimes} unspecified.
#' The function requires the specification of \code{theta} and
#' has two vital parameters, at least one of which must be specified:
#' \describe{
#' \item{\code{coxphmod}: }{Cox proportional hazards model to be used for
#' risk-adjustment. If \code{cbaseh} is not specified, it will be determined
#' from \code{coxphmod} numerically.}
#' \item{\code{cbaseh}: }{The cumulative baseline hazard rate to use for chart
#' construction. If specified with \code{coxphmod} missing, no risk-adjustment
#' will be performed}
#' }
#'
#' @references Biswas P. and Kalbfleisch J.D. (2008), A risk-adjusted CUSUM in
#' continuous time based on the Cox Model, Statistics in medicine 27, 3452-3452.
#' \doi{10.1002/sim.3216}
#'
#'
#'
#' @details The BK-CUSUM can be used to test the alternative hypothesis of an
#'  instant change of fixed size \eqn{e^\theta}{exp(\theta)}
#' in the subject specific hazard rate from \eqn{h_i(t)}{h_i(t)} to
#' \eqn{h_i(t) e^\theta}{h_i(t) exp(\theta)}. The parameter \code{C} can be used
#' to ignore the contributions of subjects, C time units after their entry
#' into the study.
#' The BK-CUSUM is constructed as
#' \deqn{G(t) = \max_{0 \leq k \leq t} \left( \theta N(k,t) - \left( e^\theta -1  \right) \Lambda(k,t)  \right),}{G(t) = max_{0 <= k <= t} (\theta N(k,t) - (e^\theta -1) \Lambda(k,t)),}
#' where \eqn{\theta}{\theta} is the log expected hazard ratio,
#' \deqn{N(k,t) = N(t) - N(k)}{N(k,t) = N(t)-N(k)}
#' with \eqn{N(t)}{N(t)} the counting process of all failures at time t, and \deqn{\Lambda(k,t) = \Lambda(t) - \Lambda(k)}{\Lambda(k,t) = \Lambda(t) - \Lambda(k)}
#' with \eqn{\Lambda(t)}{\Lambda(t)} the summed cumulative intensity of all
#'  subjects at time \eqn{t}{t}.
#'
#' @inheritParams cgr_cusum
#' @param theta The expected log-hazard ratio \eqn{\theta}{\theta} under the alternative hypothesis.
#'  If \eqn{\theta >= 0}{\theta >= 0}, the chart will try to detect an increase
#'   in hazard ratio (upper one-sided). If \eqn{\theta < 0}{\theta < 0},
#' the chart will look for a decrease in hazard ratio (lower one-sided).
#' @param twosided (optional): A boolean indicating whether a two-sided CUSUM
#'  should be constructed.
#' If \code{TRUE}, 2 CUSUM charts are determined. One to check for an increase
#' of \eqn{e^\theta}{exp(\theta)} and one for
#' a decrease of \eqn{e^{-\theta}}{exp(-\theta)} in the hazard rate w.r.t.
#'  the baseline hazard. Default is \code{FALSE}.
#'
#'
#' @return An object of class \code{bkcusum} containing:
#' \itemize{
#' \item \code{BK}: a \code{data.frame} containing the following named columns:
#' \describe{
#'   \item{\code{time}:}{times at which chart is constructed;}
#'   \item{\code{value}:}{value of the chart at corresponding times.}
#' }
#' \item \code{stopind}: indicator for whether the chart was stopped by
#' the control limit;
#' \item \code{call}: the call used to obtain output;
#' \item \code{h}: Specified value for the control limit.
#' }
#There are \code{\link[cgrcusum:plot.bkcusum]{plot}} and
# \code{\link[cgrcusum:runlength.bkcusum]{runlength}} methods for "bkcusum" objects.
#'
#' @import survival
#' @importFrom stats loess
#' @importFrom stats predict
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#'
#' @author Daniel Gomon
#' @family quality control charts
#'
#' @seealso \code{\link[success]{plot.bkcusum}}, \code{\link[success]{runlength.bkcusum}}
#'
#'
#' @examples
#' require(survival)
#' #Select only the data of the first hospital in the surgerydat data set
#' tdat <- subset(surgerydat, unit == 1)
#'
#' #We know that the cumulative baseline hazard in the data set is
#' #Exponential(0.01). If you don't know the cumulative baseline, we suggest
#' #leaving the cbaseh argument empty and determining a coxphmod (see help)
#' tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
#'
#' #Determine a risk-adjustment model using a Cox proportional hazards model.
#' #Outcome (survival) is regressed on the available covariates:
#' exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
#' tcoxmod <- coxph(exprfit, data= surgerydat)
#'
#' #Determine the values of the chart
#' bk <- bk_cusum(data = tdat, theta = log(2), coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE)
#' #plot the BK-CUSUM (exact hazard)
#' plot(bk)
#'
#' #Alternatively, cbaseh can be left empty when specifying coxphmod through coxph()
#' bk_cox <- bk_cusum(data = tdat, theta = log(2), coxphmod = tcoxmod, pb = TRUE)
#' #plot the BK-CUSUM (estimated hazard from coxph)
#' plot(bk_cox)


bk_cusum <- function(data, theta, coxphmod, cbaseh, ctimes, h, stoptime,
                    C, twosided = FALSE, pb = FALSE, assist){

  if(!missing(assist)){
    list2env(assist, envir = environment())
  }
  call = match.call()
  #-------------------------------DATA CHECKS-----------------------------#
  #Basic data checks (global for BK, CGR and Bernoulli)
  if(missing(data)){
    stop("Please provide data to construct chart.")
  } else{
    data <- check_data(data)
  }


  # determine chronological failure times
  data$otime <- data$entrytime + data$survtime
  if(!missing(C)){
    tempidx <- which(data$otime > data$entrytime + C)
    data$otime[tempidx] <- data$entrytime[tempidx] + C
    data$censorid[tempidx] <- rep(0, length(tempidx))
  }
  # determine the default construction times (all failtimes + entrytimes), if none specified
  #NOTE: WE NEED BOTH FAILTIMES AND ENTRYTIME, OTHERWISE THE VALUES MIGHT BE WRONG!!!
  if(missing(ctimes)){
    ctimes <- union(unique(data$otime), unique(data$entrytime))
  } else{
    ctimes <- union(ctimes, union(unique(data$entrytime[which(data$entrytime <= max(ctimes))]), unique(data$otime[which(data$otime <= max(ctimes))])))
  }
  ctimes <- sort(ctimes)
  if(missing(stoptime)){
    stoptime <- max(data$otime[is.finite(data$otime)])
  }
  checkcbase <- FALSE
  if(missing(coxphmod)){
    coxphmod <- NULL
  } else if(inherits(coxphmod, "coxph")){
    if(missing(cbaseh)){
      checkcbase <- TRUE
      #message("Missing cumulative baseline hazard. Determining using provided Cox PH model.")
      cbaseh <- extract_hazard(coxphmod)$cbaseh
    }
  } else if(is.list(coxphmod)){
    if(all(c("formula", "coefficients") %in% names(coxphmod))){
      checkcoxlist <- TRUE
    } else{
      stop("coxphmod does not contain $formula and/or $coefficients.")
    }
  } else if(!is.null(coxphmod)){
    stop("coxphmod is not a list or survival object.")
    }
  if(missing(cbaseh)){
    if(!checkcbase){
      stop("Please specify cbaseh (function) or coxphmod as Survival object.")
    }
  }
  if(missing(theta)){
    stop("Please specify a value for theta (ln(expected hazard ratio)).")
  }

  #Determine and sort the times at which to construct the chart
  ctimes <- sort(ctimes[which(ctimes <= stoptime)])


  #---------------------------FUNCTION BODY---------------------------#

  #How is the BK-CUSUM constructed?
  #First we calculate the value of the BK-CUSUM at first specified ctime
  #Afterwards, we calculate the increments dUt (see Biswas & Kalbfleisch (2008))
  #to calculate the value at each following ctime.
  #The problem with this is that between each 2 ctimes, we need to redetermine
  #which patients provide an active contribution to the chart.
  #Due to how the chart is constructed,
  #we require to have all patient entry and failure times to be included in ctimes.

  #A different way would be to consider the active cbaseh contribution between the ctimes separately.
  #Maybe for a future implementation?
  #Not a priority - BK-CUSUM is not computationally expensive.



  #If twosided chart is required, determine the chart in two directions
  if(twosided == TRUE){
    Gt <- matrix(0, nrow =0, ncol = 3)
    theta <- sort(c(theta, -theta))
    if(!missing(h) && length(h) == 1){
      h <- sort(c(-h, h))
    } else if(!missing(h) && length(h) == 2){
      if(!all(sign(sort(h)) == c(-1, 1))){
        stop("When specifying 2 control limits the two values should have reverse signs.")
      } else{
        h <- sort(h)
      }
    } else if(!missing(h) && length(h) > 2){
      stop("Please provide 1 or 2 values for the control limit.")
    }
  } else if (twosided == FALSE) {
    Gt <- matrix(0, nrow =0, ncol = 2)
    if(!missing(h)){
      if(length(h) > 1){
        stop("Please provide only 1 value for the control limit")
      }
      if(theta >= 0){
        h = abs(h)
      } else{
        h = -abs(h)
      }
    }
  }
  #Some booleans to keep track of events
  startval <- 0
  stopind <- FALSE
  #Progress bar
  if(pb){
    pb2 <- txtProgressBar(min= 0, max = length(ctimes), style = 3)
  }
  #Calculate subject specific risk
  riskdat <- calc_risk(data = data, coxphmod = coxphmod)

  #Keep track with index:
  j_temp <- 0

  #Calculate value at each of the required times
  for(j in seq_along(ctimes)){
    if(pb){
      setTxtProgressBar(pb2, value = j)
    }
    #As we do not have a previous value at the first time, we calculate it separately
    #We always take the first entry time of a patient
    if(j == 1){
      if(ctimes[j] >= min(data$entrytime)){
        #Which subjects are active (contribute to the cumulative intensity)
        active <- which(data$entrytime < ctimes[j] & data$otime > min(data$entrytime))
        if(length(active) > 0){
          tdat <- data[active,]
          #Determine their contribution to the total cumulative intensity Lambda(t - Si) - Lambda(prevt - Si)
          activecbaseh <- cbaseh(ifelse(tdat$otime < ctimes[j], tdat$otime, ctimes[j]))-cbaseh(min(data$entrytime))
          activecbaseh[is.na(activecbaseh)] <- 0
        } else{
          activecbaseh <- 0
        }
        #Risk-adjust the contribution in cumulative intensity
        dAT <- sum(riskdat[active] * activecbaseh)
        #How many subjects experience failure in the first time frame
        dNDT <- length(which(data$otime <= ctimes[j] & data$censorid == 1))
        #If the chart is not two-sided, construct only in one direction
      } else{
        dAT = 0
        dNDT = 0
        }
      if(twosided == FALSE){
        #Upper direction (theta >= 0), lower (theta < 0)
        #For upper, we first substract the Cumulative intensity dAT,
        #then we add the failures
        #For lower, we first add the Cumulative intensity dAT,
        #then we substract the failures.
        if(theta >= 0){
          newval_down <- max(0, 0 - (exp(theta) -1) * dAT)
          newval_up <- newval_down + theta*dNDT
          if(newval_down == newval_up){
            Gt <- rbind(Gt, c(ctimes[j], newval_up))
            j_temp <- j_temp + 1
          } else{
            Gt <- rbind(Gt, c(ctimes[j], newval_down))
            Gt <- rbind(Gt, c(ctimes[j], newval_up))
            j_temp <- j_temp + 2
          }

        } else if(theta < 0){
          newval_down <- 0 + (exp(theta)-1)* dAT
          newval_up <- newval_down - theta*dNDT
          newval_up <- min(0, newval_up)

          if(newval_down == newval_up){
            Gt <- rbind(Gt, c(ctimes[j], newval_down))
            j_temp <- j_temp + 1
          } else{
            Gt <- rbind(Gt, c(ctimes[j], newval_down))
            Gt <- rbind(Gt, c(ctimes[j], newval_up))
            j_temp <- j_temp + 2
          }
        }
      } else if(twosided == TRUE){ #Twosided BK-CUSUM
        newvalupper_down <- max(0, 0 - (exp(theta[2])-1)* dAT)
        newvalupper_up <- newvalupper_down + theta[2]*dNDT
        newvallower_down <- 0 + (exp(theta[1])-1)* dAT
        newvallower_up <- newvallower_down - theta[1]*dNDT
        newvallower_up <- min(0, newvallower_up)
        if((newvalupper_down == newvalupper_up) & (newvallower_down == newvallower_up)){
          Gt <- rbind(Gt, c(ctimes[j], newvalupper_up, newvallower_down))
          j_temp <- j_temp + 1
        } else{
          Gt <- rbind(Gt, c(ctimes[j], newvalupper_down, newvallower_down))
          Gt <- rbind(Gt, c(ctimes[j], newvalupper_up, newvallower_up))
          j_temp <- j_temp + 2
        }
      }
    } else{
      #Determine dUt from ctimes[j-1] to ctimes[j]
      active <- which(data$entrytime < ctimes[j] & data$otime > ctimes[j-1])
      tdat <- data[active,]

      if(length(active) > 0){
        #Determine amount of (relevant) failures in this subset
        dNDT <- length(which(tdat$otime <= ctimes[j] & tdat$censorid == 1))
        #Contribution of active subjects to A(t)
        activecbaseh <- cbaseh(ifelse(tdat$otime < ctimes[j], tdat$otime, ctimes[j]) - tdat$entrytime)-cbaseh(ctimes[j-1] - tdat$entrytime)
        #Error check for empty contribution
        activecbaseh[is.na(activecbaseh)] <- 0
        #Determine dA(t) (have to risk-adjust)
        dAT <- sum(riskdat[active] * activecbaseh)
      } else{
        dAT <- 0
        dNDT <- 0
      }

      #dUt <- theta * dNDT - (exp(theta) - 1) * dAT
      #Determine cusum value and rbind to previous values
      if(twosided == FALSE){
        if(theta >= 0){
          newval_down <- max(0, Gt[j_temp,2] - (exp(theta)-1)* dAT)
          newval_up <- newval_down + theta*dNDT
          if(newval_down == newval_up){
            Gt <- rbind(Gt, c(ctimes[j], newval_up))
            j_temp <- j_temp + 1
          } else{
            Gt <- rbind(Gt, c(ctimes[j], newval_down))
            Gt <- rbind(Gt, c(ctimes[j], newval_up))
            j_temp <- j_temp + 2
          }
        } else if(theta < 0){
          newval_down <- Gt[j_temp,2] + (exp(theta)-1)* dAT
          newval_up <- newval_down - theta*dNDT
          newval_up <- min(0, newval_up)
          if(newval_down == newval_up){
            Gt <- rbind(Gt, c(ctimes[j], newval_down))
            j_temp <- j_temp + 1
          } else{
            Gt <- rbind(Gt, c(ctimes[j], newval_down))
            Gt <- rbind(Gt, c(ctimes[j], newval_up))
            j_temp <- j_temp + 2
          }
        }
      } else if(twosided == TRUE){
        newvalupper_down <- max(0, Gt[j_temp,2] - (exp(theta[2])-1)* dAT)
        newvalupper_up <- newvalupper_down + theta[2]*dNDT
        newvallower_down <- Gt[j_temp,3] + (exp(theta[1])-1)* dAT
        newvallower_up <- newvallower_down - theta[1]*dNDT
        newvallower_up <- min(0, newvallower_up)
        if((newvalupper_down == newvalupper_up) & (newvallower_down == newvallower_up)){
          Gt <- rbind(Gt, c(ctimes[j], newvalupper_up, newvallower_down))
          j_temp <- j_temp + 1
        } else{
          Gt <- rbind(Gt, c(ctimes[j], newvalupper_down, newvallower_down))
          Gt <- rbind(Gt, c(ctimes[j], newvalupper_up, newvallower_up))
          j_temp <- j_temp + 2
        }
      }
    }
    if (!missing(h)){
      if(twosided == TRUE){
        if(length(h) == 2){
          if( (Gt[j_temp,2] >= h[2]) | (Gt[j_temp,3] <= h[1]) ) {stopind = TRUE; break}
        } else if(length(h) == 1){
          if( (abs(Gt[j_temp,2]) >= abs(h)) | (abs(Gt[j_temp,3]) >= abs(h)) ) {stopind = TRUE; break}
        }
      } else{
        if( abs(Gt[j_temp,2]) >= abs(h) ) {stopind = TRUE; break}
      }
    }
  }
  if(twosided == FALSE){
    colnames(Gt) <- c("time", "value")
  } else{
    colnames(Gt) = c("time", "val_up", "val_down")
  }
  if(pb){
    close(pb2)
  }
  Gt <- as.data.frame(Gt)
  bkcus <- list(BK = Gt,
              stopind = stopind,
              call = call)
  if(!missing(h)){bkcus$h <- h}
  class(bkcus) <- "bkcusum"
  bkcus
}
