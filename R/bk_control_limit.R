#' @title Determine control limits for BK-CUSUM by simulation
#'
#' @description This function can be used to determine control limits for the
#' BK-CUSUM (\code{\link[success]{bk_cusum}}) procedure by restricting the type I error \code{alpha} of the
#' procedure over \code{time}.
#'
#' @details This function performs 3 steps to determine a suitable control limit.
#' \itemize{
#' \item Step 1: Generates \code{n_sim} in-control units (failure rate as baseline).
#' If \code{data} is provided, subject covariates are resampled from the data set.
#' \item Step 2: Determines chart values for all simulated units.
#' \item Step 3: Determines control limits such that at most a proportion \code{alpha}
#' of all units cross the control limit.
#' } The generated data as well as the charts are also returned in the output.
#'
#'
#'
#' @inheritParams bk_cusum
#' @param time A numeric value over which the type I error \code{alpha} must be restricted.
#' @param alpha A proportion between 0 and 1 indicating the required maximal type I error.
#' Default is 0.05.
#' @param psi A numeric value indicating the estimated Poisson arrival rate of subjects
#' at their respective units. Can be determined using
#' \code{\link[success:parameter_assist]{parameter_assist()}}.
#' @param n_sim An integer value indicating the amount of units to generate for the
#' determination of the control limit. Larger values yield more precise control limits,
#' but greatly increase computation times. Default is 20.
#' @param inv_cbaseh (optional): A function that returns the unadjusted inverse cumulative
#' baseline
#' hazard \eqn{H^{-1}_0(t)}{H_0^-1(t)}. If \code{inv_cbaseh} is missing, it will be
#' determined from \code{cbaseh} numerically.
#' @param baseline_data (optional): A \code{data.frame} used for covariate resampling
#' with rows representing subjects and at least the
#' following named columns: \describe{
#'   \item{\code{entrytime}:}{time of entry into study (numeric);}
#'   \item{\code{survtime}:}{time from entry until event (numeric);}
#'   \item{\code{censorid}:}{censoring indicator (0 = right censored, 1 = observed),
#'    (integer).}
#' } and optionally additional covariates used for risk-adjustment. Can only be specified
#'  in combination with \code{coxphmod}.
#' @param interval (optional): Interval in which survival times should be solved for numerically.
#' @param h_precision (optional): A numerical value indicating how precisely the control limit
#' should be determined. By default, control limits will be determined up to 2 significant digits.
#' @param seed (optional): A numeric seed for survival time generation. Default
#' is 01041996 (my birthday).
#' @param pb (optional): A boolean indicating whether a progress bar should
#' be shown. Default is \code{FALSE}.
#' @param chartpb (optional): A boolean indicating whether progress bars should
#' be displayed for the constructions of the charts. Default is \code{FALSE}.
#'
#'
#' @return A list containing three components:
#' \itemize{
#' \item \code{call}: the call used to obtain output;
#' \item \code{charts}: A list of length \code{n_sim} containing the constructed charts;
#' \item \code{data}: A \code{data.frame} containing the in-control generated data.
#' \item \code{h}: Determined value of the control limit.
#' }
# There are \code{\link[cgrcusum:plot.cgrcusum]{plot}} and
#  \code{\link[cgrcusum:runlength.cgrcusum]{runlength}} methods for "cgrcusum" objects.
#'
#' @export
#'
#' @author Daniel Gomon
#' @family control limit simulation
#' @seealso \code{\link[success]{plot.bkcusum}}, \code{\link[success]{runlength.bkcusum}}
#'
#'
#' @examples
#' require(survival)
#' \dontrun{
#' require(survival)
#' exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
#' tcoxmod <- coxph(exprfit, data= surgerydat)
#'
#' a <- bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
#' cbaseh = function(t) chaz_exp(t, lambda = 0.02),
#' inv_cbaseh = function(t) inv_chaz_exp(t, lambda = 0.02), psi = 0.5,
#' n_sim = 10)
#'
#' b <- bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
#' coxphmod = tcoxmod, psi = 0.5, n_sim = 10)
#' }

















bk_control_limit <- function(time, alpha = 0.05, psi, n_sim = 20, theta,
                             coxphmod, cbaseh, inv_cbaseh, baseline_data,
                             interval = c(0, 9e12), h_precision = 0.01,
                             seed = 1041996, pb = FALSE, chartpb = FALSE){
  #This function consists of 3 steps:
  #1. Constructs n_sim instances (hospitals) with subject arrival rate psi and
  #   cumulative baseline hazard cbaseh. Possibly by resampling subject charac-
  #   teristics from data and risk-adjusting using coxphmod.
  #2. Construct the CGR-CUSUM chart for each hospital until timepoint time
  #3. Determine control limit h such that at most proportion alpha of the
  #   instances will produce a signal.
  call = match.call()
  set.seed(seed)
  manualcbaseh <- FALSE

  #First we generate the n_sim unit data
  message("Step 1/3: Generating in-control data.")
  df_temp <- generate_units(time = time, psi = psi, n_sim = n_sim, cbaseh = cbaseh,
                            inv_cbaseh = inv_cbaseh, coxphmod = coxphmod,
                            baseline_data = baseline_data, interval = interval)


  message("Step 2/3: Determining BK-CUSUM chart(s).")
  if(!missing(coxphmod) & missing(baseline_data)){
    #We don't want to use risk-adjustment if no baseline_data specified
    cbaseh <- extract_hazard(coxphmod)$cbaseh
    coxphmod <- NULL
  } else if(!missing(coxphmod)){
    #Otherwise, if we do want to use risk-adjust: calculate cbaseh once.
    cbaseh <- extract_hazard(coxphmod)$cbaseh
  } else if(missing(coxphmod)){
    coxphmod <- NULL
  }

  #Construct for each unit a CGR-CUSUM until time
  BK_CUSUMS <- list(length = n_sim)
  if(pb){
    pbbar <- pbapply::timerProgressBar(min = 1, max = n_sim)
    on.exit(close(pbbar))
  }

  for(j in 1:n_sim){
    if(pb){
      pbapply::setTimerProgressBar(pbbar, value = j)
    }
    BK_CUSUMS[[j]] <- bk_cusum(data = subset(df_temp, unit == j),
                                 coxphmod = coxphmod, theta = theta, cbaseh = cbaseh,
                                 stoptime = time,
                                 pb = chartpb)
  }


  message("Step 3/3: Determining control limits")

  #Create a sequence of control limit values h to check for
  #start from 0.1 to maximum value of all CGR-CUSUMS
  CUS_max_val <- 0
  for(k in 1:n_sim){
    temp_max_val <- max(BK_CUSUMS[[k]]$BK["value"])
    if(temp_max_val >= CUS_max_val){
      CUS_max_val <- temp_max_val
    }
  }
  hseq <- rev(seq(from = h_precision, to = CUS_max_val + h_precision, by = h_precision))

  #Determine control limits using runlength
  control_h <- CUS_max_val
  for(i in seq_along(hseq)){
    #Determine type I error using current h
    typ1err_temp <- sum(sapply(BK_CUSUMS, function(x) is.finite(runlength(x, h = hseq[i]))))/n_sim
    if(typ1err_temp <= alpha){
      control_h <- hseq[i]
    } else{
      break
    }
  }


  return(list(call = call,
              charts = BK_CUSUMS,
              data = df_temp,
              h = control_h))
}

