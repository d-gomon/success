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
#' @inheritParams bernoulli_cusum
#' @param time A numeric value over which the type I error \code{alpha} must be restricted.
#' @param alpha A proportion between 0 and 1 indicating the required maximal type I error.
#' @param psi A numeric value indicating the estimated Poisson arrival rate of subjects
#' at their respective units. Can be determined using
#' \code{\link[success:parameter_assist]{parameter_assist()}}.
#' @param n_sim An integer value indicating the amount of units to generate for the
#' determination of the control limit. Larger values yield more precise control limits,
#' but greatly increase computation times. Default is 20.
#' @param baseline_data (optional): A \code{data.frame} used for covariate resampling
#' with rows representing subjects and at least the
#' following named columns: \describe{
#'   \item{\code{entrytime}:}{time of entry into study (numeric);}
#'   \item{\code{survtime}:}{time from entry until event (numeric);}
#'   \item{\code{censorid}:}{censoring indicator (0 = right censored, 1 = observed),
#'    (integer).}
#' } and optionally additional covariates used for risk-adjustment. Can only be specified
#'  in combination with \code{coxphmod}.
#' @param h_precision (optional): A numerical value indicating how precisely the control limit
#' should be determined. By default, control limits will be determined up to 2 significant digits.
#' @param seed (optional): A numeric seed for survival time generation.
#' Default is 01041996 (my birthday).
#' @param pb (optional): A boolean indicating whether a progress bar should
#' be shown. Default is \code{FALSE}.
#'
#'
#' @return A list containing three components:
#' \itemize{
#' \item \code{call}: the call used to obtain output;
#' \item \code{charts}: A list of length \code{n_sim} containing the constructed charts;
#' \item \code{data}: A \code{data.frame} containing the in-control generated data.
#' \item \code{h}: Determined value of the control limit.
#' }
# There are \code{\link[success:plot.success]{plot}} and
#  \code{\link[success:runlength.success]{runlength}} methods for "success" objects.
#'
#' @export
#'
#' @author Daniel Gomon
#' @family control limit simulation
#' @seealso \code{\link[success]{plot.bkcusum}}, \code{\link[success]{runlength.bkcusum}}
#'
#'
#' @examples
#' \dontrun{
#' #We consider patient outcomes 100 days after their entry into the study.
#' followup <- 100
#' #Determine a risk-adjustment model using a generalized linear model.
#' #Outcome (failure within 100 days) is regressed on the available covariates:
#' exprfitber <- as.formula("(survtime <= followup) & (censorid == 1)~ age + sex + BMI")
#' glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
#' a <- bernoulli_control_limit(time = 500, alpha = 0.1, followup = followup,
#'  psi = 0.5, n_sim = 10, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat)
#'
#' }





bernoulli_control_limit <- function(time, alpha = 0.05, followup, psi, n_sim = 20,
                                    theta, p0, p1, glmmod, baseline_data,
                                    h_precision = 0.01,
                                    seed = 1041996, pb = FALSE){
  #This function consists of 3 steps:
  #1. Constructs n_sim instances (hospitals) with subject arrival rate psi and
  #   cumulative baseline hazard cbaseh. Possibly by resampling subject charac-
  #   teristics from data and risk-adjusting using coxphmod.
  #2. Construct the CGR-CUSUM chart for each hospital until timepoint time
  #3. Determine control limit h such that at most proportion alpha of the
  #   instances will produce a signal.
  call = match.call()
  set.seed(seed)

  #First we generate the n_sim unit data
  message("Step 1/3: Generating in-control data.")
  df_temp <- generate_units_bernoulli(time = time, psi = psi, n_sim = n_sim,
                                      p0 = p0, p1 = p1, theta = theta, glmmod = glmmod,
                                      followup = followup,
                            baseline_data = baseline_data)


  message("Step 2/3: Determining Bernoulli CUSUM chart(s).")

  #Construct for each unit a Bernoulli CUSUM until time
  charts <- list(length = n_sim)
  if(pb){
    pbbar <- pbapply::timerProgressBar(min = 1, max = n_sim)
    on.exit(close(pbbar))
  }

  for(j in 1:n_sim){
    if(pb){
      pbapply::setTimerProgressBar(pbbar, value = j)
    }
    charts[[j]] <- bernoulli_cusum(data = subset(df_temp, unit == j),
                                   followup = followup, theta = theta,
                                   glmmod = glmmod,
                                   p0 = p0, p1 = p1, stoptime = time)
  }


  message("Step 3/3: Determining control limits")

  #Create a sequence of control limit values h to check for
  #start from 0.1 to maximum value of all CGR-CUSUMS
  CUS_max_val <- 0
  for(k in 1:n_sim){
    temp_max_val <- max(charts[[k]]$CUSUM["value"])
    if(temp_max_val >= CUS_max_val){
      CUS_max_val <- temp_max_val
    }
  }
  hseq <- rev(seq(from = h_precision, to = CUS_max_val + h_precision, by = h_precision))

  #Determine control limits using runlength
  control_h <- CUS_max_val
  for(i in seq_along(hseq)){
    #Determine type I error using current h
    typ1err_temp <- sum(sapply(charts, function(x) is.finite(runlength(x, h = hseq[i]))))/n_sim
    if(typ1err_temp <= alpha){
      control_h <- hseq[i]
    } else{
      break
    }
  }


  return(list(call = call,
              charts = charts,
              data = df_temp,
              h = control_h))
}

