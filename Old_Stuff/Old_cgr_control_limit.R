#' @title Determine control limits for CGR-CUSUM by simulation
#'
#' @description This function can be used to determine control limits for the
#' CGR-CUSUM procedure by restricting the type I error \code{alpha} of the
#' procedure over \code{time}.
#'
#' @details This function performs 3 steps to determine a suitable control limit.
#' \itemize{
#' \item Step 1: Generates \code{n_sim} in-control units (failure rate as baseline).
#' If \code{data} is provided, subject covariates are resampled from the data set.
#' \item Step 2: Determines chart values for all simulated units.
#' \item Step 3: Determines control limits such that at most a proportion \code{alpha}
#' of all units are cross the control limit.
#' } The generated data as well as the charts are also returned in the output.
#'
#'
#' @param time A numeric value over which the type I error \code{alpha} must be restricted.
#' @param alpha A proportion between 0 and 1 indicating the required maximal type I error.
#' @param psi A numeric value indicating the estimated Poisson arrival rate of subjects
#' at their respective units. Can be determined using
#' \code{\link[cgrcusum:parameter_assist]{parameter_assist()}}.
#' @param n_sim An integer value indicating the amount of units to generate for the
#' determination of the control limit. Larger values yield more precise control limits,
#' but greatly increase computation times. Default is 20.
#' @param cbaseh (optional): A function that returns the unadjusted cumulative
#' baseline hazard \eqn{H_0(t)}{H_0(t)}. If \code{cbaseh} is missing but
#' \code{coxphmod} has been
#' specified as a survival object, this baseline hazard rate will be determined
#' using the provided \code{coxphmod}.
#' @param inv_cbaseh (optional): A function that returns the unadjusted inverse cumulative
#' baseline
#' hazard \eqn{H^{-1}_0(t)}{H_0^-1(t)}. If \code{inv_cbaseh} is missing, it will be
#' determined from \code{cbaseh} numerically.
#' @param coxphmod (optional): A cox proportional hazards regression model as
#' produced by
#' the function \code{\link[survival:coxph]{coxph()}}. Suggested: \cr
#' \code{coxph(Surv(survtime, censorid) ~ covariates, data = baseline_data)}. \cr
#' Alternatively, a list with the following elements:
#' \describe{
#'   \item{\code{formula}:}{a \code{\link[stats:formula]{formula()}} in the form \code{~ covariates};}
#'   \item{\code{coefficients}:}{a named vector specifying risk adjustment coefficients
#'   for covariates. Names must be the same as in \code{formula} and colnames of \code{data}.}
#' }
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
#' @param ncores (optional): Number of cores to use to parallelize the computation of the
#' CGR-CUSUM charts. If ncores = 1 (default), no parallelization is done. You
#' can use \code{\link[parallel:detectCores]{detectCores()}} to check how many
#' cores are available on your computer.
#' @param seed (optional): A numeric seed for survival time generation. Default = my birthday.
#' @param pb (optional): A boolean indicating whether a progress bar should
#' be shown. Default is \code{FALSE}.
#' @param chartpb (optional): A boolean indicating whether progress bars should
#' be displayed for the constructions of the charts.
#'
#'
#' @return A list containing three components:
#' \itemize{
#' \item \code{call}: the call used to obtain output;
#' \item \code{CHARTS}: A list of length \code{n_sim} containing the constructed charts;
#' \item \code{data}: A \code{data.frame} containing the in-control generated data.
#' \item \code{h}: Determined value of the control limit.
#' }
# There are \code{\link[cgrcusum:plot.cgrcusum]{plot}} and
#  \code{\link[cgrcusum:runlength.cgrcusum]{runlength}} methods for "cgrcusum" objects.
#'
#' @export
#'
#' @author Daniel Gomon
#' @family quality control charts
#' @seealso \code{\link[cgrcusum]{plot.cgrcusum}}, \code{\link[cgrcusum]{runlength.cgrcusum}}
#'
#'
#' @examples
#' require(survival)
#' \dontrun{
#' a <- cgr_control_limit(time = 500, alpha = 0.1, cbaseh = function(t) chaz_exp(t, lambda = 0.02),
#' inv_cbaseh = function(t) inv_chaz_exp(t, lambda = 0.02), psi = 0.5, n_sim = 10)
#'
#' b <- cgr_control_limit(time = 500, alpha = 0.1, coxphmod = tcoxmod, psi = 0.5, n_sim = 10,
#' data = subset(surgerydat, hosp_num == 1))
#' }






#inv_cbaseh should be specified TOGETHER with cbaseh
#If only inv_cbaseh is specified, then the CUSUM cannot be calculated
#(or only if coxphmod is present - but I don't want to program that)











cgr_control_limit <- function(time, alpha, psi, n_sim = 20, cbaseh, inv_cbaseh,
                              coxphmod = NULL, baseline_data, interval = c(0, 9e12),
                              h_precision = 0.01,
                              ncores = 1, seed = 1041996, pb = FALSE, chartpb = FALSE){
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

  message("Step 1/3: Generating in-control data.")

  #Generate n_sim instances
  #We create a data frame to contain all information
  entrytime_temp <- numeric(0)
  unit_temp <- numeric(0)
  #Determine amount of required significant digits
  if(!missing(baseline_data)){
    signif_dig <- tryCatch({max(sapply(baseline_data$entrytime, function(x) match(TRUE, round(x, 1:20) == x)))-1},
                           error = function(cond){2})
  } else{
    signif_dig <- 2
  }
  for(i in 1:n_sim){
    arrivtimes <- round(gen_arriv_times(psi = psi, t = time), digits = signif_dig)
    entrytime_temp <- c(entrytime_temp, arrivtimes)
    unit_temp <- c(unit_temp, rep(i, length(arrivtimes)))
  }

  #Resample covariates from existing data if supplied
  if(!missing(baseline_data) & !missing(coxphmod)){
    df_temp <- baseline_data[sample(nrow(baseline_data), length(entrytime_temp), replace = TRUE),]
    df_temp$entrytime <- entrytime_temp
    df_temp$unit <- unit_temp
  } else if(missing(baseline_data) & !missing(coxphmod)){ #If coxphmod is specified, but no data, determine cbaseh but don't resample
    manualcbaseh <- TRUE
    #Determine cbaseh
    cbase_temp <- basehaz(coxphmod, centered = FALSE)
    #Set the cumulative hazard at time 0 to 0
    if(length(which(cbase_temp$time == 0)) != 0){
      cbase_temp$hazard[which(cbase_temp$time == 0)] <- 0
    } else{
      rbind(c(0, 0), cbase_temp)
    }
    max_haz <- max(cbase_temp$hazard)
    max_time <- max(cbase_temp$time)
    # cbase_loess <- loess(cbase_temp$hazard~cbase_temp$time)
    # cbaseh <- function(t) predict(cbase_loess, t)
    cbaseh <- approxfun(x = cbase_temp$time, y = cbase_temp$hazard)
    interval[2] <- max(cbase_temp$time)

    inv_cbaseh_temp <- function(y){
      if(y < max_haz){
        return(RootLinearInterpolant(x = cbase_temp$time, y = cbase_temp$hazard, y0 = y))
      } else{
        return(max_time)
      }
    }
    inv_cbaseh <- Vectorize(inv_cbaseh_temp)
    #We don't want to use risk-adjustment as data is not specified.
    coxphmod <- NULL
    df_temp <- data.frame("entrytime" = entrytime_temp, "unit" = unit_temp)
  } else{ #Else just create a data frame with arrival times and unit number
    df_temp <- data.frame("entrytime" = entrytime_temp, "unit" = unit_temp)
  }
  df_temp$survtime <- rep(NA, nrow(df_temp))
  #df_temp is now the dataframe containing all entrytimes and covariates for subjects

  #Calculate relative risk to use for adjustment
  #relative_risk <- calc_risk(df_temp, coxphmod = coxphmod)


  #Determine survival times (3 possibilities)
  if(!missing(inv_cbaseh)){
    #We can determine outcome directly
    if(missing(baseline_data)){
      survtimes <- gen_surv_times(invchaz = inv_cbaseh, data = nrow(df_temp), coxphmod = coxphmod)
    } else{
      survtimes <- gen_surv_times(invchaz = inv_cbaseh, data = df_temp, coxphmod = coxphmod)
    }
  } else{ #if inv_cbaseh is not specified, we determine outcome times by solving S(t) - U
    if(!is.null(coxphmod) & missing(cbaseh)){
      manualcbaseh <- TRUE
      #Determine cbaseh
      cbase_temp <- basehaz(coxphmod, centered = FALSE)
      #Set the cumulative hazard at time 0 to 0
      #Set the cumulative hazard at time 0 to 0
      if(length(which(cbase_temp$time == 0)) != 0){
        cbase_temp$hazard[which(cbase_temp$time == 0)] <- 0
      } else{
        rbind(c(0, 0), cbase_temp)
      }
      max_haz <- max(cbase_temp$hazard)
      max_time <- max(cbase_temp$time)
      # cbase_loess <- loess(cbase_temp$hazard~cbase_temp$time)
      # cbaseh <- function(t) predict(cbase_loess, t)
      cbaseh <- approxfun(x = cbase_temp$time, y = cbase_temp$hazard)
      interval[2] <- max(cbase_temp$time)

      inv_cbaseh_temp <- function(y){
        if(y < max_haz){
          return(RootLinearInterpolant(x = cbase_temp$time, y = cbase_temp$hazard, y0 = y))
        } else{
          return(max_time)
        }
      }
      inv_cbaseh <- Vectorize(inv_cbaseh_temp)
    } else if(!missing(cbaseh)){
      #Determine the inverse cumulative baseline hazard H^{-1}_0(t)
      inv_cbaseh_temp <- function(y, lower = 0, upper = interval[2]){
        tryCatch({return(unname(unlist(uniroot((function (x) cbaseh(x) - y), lower = lower, upper = upper)$root)))}, error = function(cond){ return(upper)})
      }

      #inv_cbaseh_temp2 <- function(y, lower = 0, upper = 9E12){
      #  uniroot((function (x) cbaseh(x) - y), lower = lower, upper = upper)$root
      #}
      inv_cbaseh <- Vectorize(inv_cbaseh_temp)
    } else if(missing(cbaseh)){
      stop("Please specify either coxphmod, cbaseh or inv_cbaseh.")
    }
    #We now determine the survival time by solving S(t) - U = 0 or
    #log(S(t)) - log(U) = 0
    # rootfun <- function(t, x, cbaseh){
    #   return((-cbaseh(t) * relative_risk[x]) - log(stats::runif(1)))
    # }
    # find_root <- function(x, cbaseh){
    #   stats::uniroot(f = function(t) rootfun(t, x, cbaseh = cbaseh),
    #                  interval = interval, check.conv = TRUE)$root
    # }


    if(missing(baseline_data)){
      survtimes <- gen_surv_times(invchaz = inv_cbaseh, data = nrow(df_temp), coxphmod = coxphmod)
    } else{
      survtimes <- gen_surv_times(invchaz = inv_cbaseh, data = df_temp, coxphmod = coxphmod)
    }
  }
  df_temp$survtime <- survtimes
  df_temp$censorid <- rep(1, nrow(df_temp))
  #If we determined cbaseh manually, censor all observations which were generated at maximum failure time
  if(manualcbaseh == TRUE){
    df_temp$censorid <- ifelse(survtimes == max(cbase_temp$time), 0, 1)
  }

  message("Step 2/3: Determining CGR-CUSUM chart.")
  if(missing(coxphmod)){
    coxphmod <- NULL
  }

  #Construct for each unit a CGR-CUSUM until time
  CGR_CUSUMS <- list(length = n_sim)
  if(pb){
    pbbar <- pbapply::timerProgressBar(min = 1, max = n_sim)
    on.exit(close(pbbar))
  }

  for(j in 1:n_sim){
    if(pb){
      pbapply::setTimerProgressBar(pbbar, value = j)
    }
    CGR_CUSUMS[[j]] <- cgr_cusum(data = subset(df_temp, unit == j), coxphmod = coxphmod, cbaseh = cbaseh,
                                   stoptime = time, ncores = ncores,
                                 pb = chartpb)
  }


  message("Step 3/3: Determining control limits")

  #Create a sequence of control limit values h to check for
  #start from 0.1 to maximum value of all CGR-CUSUMS
  CUS_max_val <- 0
  for(k in 1:n_sim){
    temp_max_val <- max(CGR_CUSUMS[[k]]$CGR["value"])
    if(temp_max_val >= CUS_max_val){
      CUS_max_val <- temp_max_val
    }
  }
  hseq <- rev(seq(from = h_precision, to = CUS_max_val + h_precision, by = h_precision))

  #Determine control limits using runlength
  control_h <- CUS_max_val
  for(i in seq_along(hseq)){
    #Determine type I error using current h
    typ1err_temp <- sum(sapply(CGR_CUSUMS, function(x) is.finite(runlength(x, h = hseq[i]))))/n_sim
    if(typ1err_temp <= alpha){
      control_h <- hseq[i]
    } else{
      break
    }
  }


  return(list(call = call,
              CUSUMS = CGR_CUSUMS,
              data = df_temp,
              h = control_h))
}

# Original function to find cbaseh inverse.
# testfct <- function(t){
#   return(exp(-chaz_exp(t, lambda = 0.02)) - runif(1))
# }
#
# uniroot(f = function(t) testfct(t), interval = c(0, 8000))$root
#
# invcbaseh <- function(y, lower = 0, upper = 6 * 365){
#   tryCatch({return(unname(unlist(uniroot((function (x) cbaseh(x) - y), lower = lower, upper = upper)[1])))}, error = function(cond){ return(upper)})
# }


