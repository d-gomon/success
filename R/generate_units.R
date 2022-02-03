#' @noRd
#' @keywords internal
#'
#' @author Daniel Gomon
#'
#' @importFrom stats uniroot
#'
#'
#'
#'


generate_units <- function(time, psi, n_sim = 20, cbaseh, inv_cbaseh,
                           coxphmod = NULL, baseline_data, interval = c(0, 9e12)){

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
    survtimes <- gen_surv_times(invchaz = inv_cbaseh, data = nrow(df_temp), coxphmod = coxphmod)
  } else{
    survtimes <- gen_surv_times(invchaz = inv_cbaseh, data = df_temp, coxphmod = coxphmod)
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
