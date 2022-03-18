#' @title Extract (inverse) cumulative baseline hazard from Cox PH model
#'
#' @description Extracts a function which returns the (inverse) cumulative
#' baseline hazard from a \code{\link[survival:coxph]{coxph()}} call.
#'
#'
#' @param coxphmod A call to \code{\link[survival:coxph]{coxph()}}.
#'
#'
#' @author Daniel Gomon
#'
#' @seealso \code{\link[survival:coxph]{coxph}}
#'
#' @returns A list containing:
#' \itemize{
#' \item \code{cbaseh}: A function which returns the cumulative baseline hazard
#' at specified time;
#' \item \code{inv_cbaseh}: A function which returns the inverse cumulative
#' baseline hazard at specified time.
#' \item \code{max_time}: maximal time at which \code{cbaseh} is known;
#' \item \code{max_haz}: value of maximal hazard (at maximum time).
#' }
#'
#' @importFrom stats approxfun
#' @export
#'
#' @examples
#' require(survival)
#' exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
#' tcoxmod <- coxph(exprfit, data= surgerydat)
#' tcox_hazard_fcts <- extract_hazard(tcoxmod)
#'
#'



extract_hazard <- function(coxphmod){
  #We extract a list of time and chaz from the cox PH model
  cbase_temp <- basehaz(coxphmod, centered = FALSE)

  #We check to make sure that the cumulative hazard is defined at time 0
  if(length(which(cbase_temp$time == 0)) == 0){
    cbase_temp <- rbind(c(0, 0), cbase_temp)
  }

  #We do not know what happens outside the intervals, so determine max values.
  max_haz <- max(cbase_temp$hazard)
  max_time <- max(cbase_temp$time)


  cbaseh <- approxfun(x = cbase_temp$time, y = cbase_temp$hazard, yleft = 0, yright = max_haz)



  #Inverse cumulative hazard can be solved numerically from linear interpolant
  inv_cbaseh_temp <- function(y){
    #If we know value, determine inverse
    if(y < max_haz){
      out <- RootLinearInterpolant(x = cbase_temp$time, y = cbase_temp$hazard, y0 = y)
    } else{ #Otherwise, simply return largest value (can censor afterwards)
      out <- max_time
    }
    #Can happen if cbaseh does not reach 0!
    if(is.na(out)){
      out <- 0
    }
    return(out)
  }
  #Vectorize the function to handle vector inputs (important for surv time
  #generation)
  inv_cbaseh <- Vectorize(inv_cbaseh_temp)

  return(list(cbaseh = cbaseh,
              inv_cbaseh = inv_cbaseh,
              max_time = max_time,
              max_haz = max_haz))
}




