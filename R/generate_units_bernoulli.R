#' @noRd
#' @keywords internal
#'
#' @author Daniel Gomon
#'
#' @importFrom stats rbinom
#'
#' @inheritParams bernoulli_control_limit
#'

generate_units_bernoulli <- function(time, psi, n_sim = 20, p0, p1, theta,
                           glmmod = NULL, followup, baseline_data){

  manualglm <- FALSE

  #Determine amount of required significant digits
  if(!missing(baseline_data)){
    signif_dig <- tryCatch({max(sapply(baseline_data$entrytime,
                                       function(x) match(TRUE, round(x, 1:20) == x)))-1},
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


  #------------------------Step 2: generate survtimes ----------------
  if(!missing(glmmod)){
    surv_probabilities <- predict(glmmod, newdata = df_temp, type = "response")
    survtimes <- ifelse(rbinom(n = nrow(df_temp), size = 1, prob = surv_probabilities) == 1,
                        followup/2, followup + 1)
  } else if(!missing(p0)){
    surv_probabilities <- rbinom(n = nrow(df_temp), size = 1, prob = p0)
    survtimes <- ifelse(rbinom(n = nrow(df_temp), size = 1, prob = surv_probabilities) == 1,
                        followup/2, followup + 1)
  }

  df_temp$survtime <- survtimes
  df_temp$censorid <- rep(1, nrow(df_temp))



  return(df_temp)
}
