#' Determine run length of a CUSUM chart
#'
#' This function can be used to calculate the run length of a 'cgrcusum', 'bkcusum'
#' or 'bercusum'
#' chart when using control limit h
#'
#'
#' @param chart A \code{cgrcusum}, \code{bkcusum} or \code{bercusum} chart.
#' @param h Control limit h to be used when determining the run length
#' @param ... Other parameters
#'
#' @return The run length of the chart with the given control limit.
#'
#'
#' @author Daniel Gomon
#' @examples
#' exprfitber <- as.formula("(survtime <= 100) & (censorid == 1) ~ age + sex + BMI")
#' glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
#' bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 14), glmmod = glmmodber,
#'                    followup = 100, theta = log(2))
#' #Determine the run length of the above Bernoulli CUSUM when using a control limit
#' #of h = 1.
#' runlength(bercus, h = 1)


#' @export
runlength <- function(chart, h){
  UseMethod("runlength")
}


#' @describeIn runlength determines runlength of \code{cgrcusum} object
#' @export
runlength.cgrcusum <- function(chart, h, ...){
  if(missing(h)){
    stop("Please specify a control limit h.")
  } else if(!all(is.numeric(h), length(h) == 1)){
    stop("Control limit h must be a single numeric value.")
  }
  if(missing(chart)){
    stop("Please provide a 'cgrcusum', 'bkcusum' or 'bercusum' chart as input.")
  }
  ind <- which(abs(chart$CGR[,2]) >= abs(h))[1]
  if(is.na(ind)){
    return(Inf)
  } else{
    return(chart$CGR[ind,1] - chart$CGR[1, "time"])
  }
}





#' @describeIn runlength determines runlength of \code{bkcusum} object
#' @export
runlength.bkcusum <- function(chart, h, ...){
  if(missing(h)){
    stop("Please specify a control limit h.")
  } else if(!is.numeric(h)){
    stop("Control limit h must be numeric.")
  }
  if(missing(chart)){
    stop("Please provide a 'cgrcusum', 'bkcusum' or 'bercusum' chart as input.")
  }
  if(isFALSE(chart$call[["twosided"]]) | is.null(chart$call[["twosided"]])){
    if(!length(h) == 1){
      stop("Please provide only 1 control limit.")
    }
    ind <- which(abs(chart$BK[,2]) >= abs(h))[1]
    if(is.na(ind)){
      return(Inf)
    } else{
      return(chart$BK[ind,1] - chart$BK[1, "time"])
    }
  } else{ #twosided BK-CUSUM
    if(length(h) == 1){
      if(h < 0){
        ind <- which(chart$BK$val_down <= h)[1]
        if(is.na(ind)){
          return(Inf)
        } else{
          return(chart$BK[ind,1] - chart$BK[1, "time"])
        }
      } else if(h >= 0){
        ind <- which(chart$BK$val_up >= h)[1]
        if(is.na(ind)){
          return(Inf)
        } else{
          return(chart$BK[ind,1] - chart$BK[1, "time"])
        }
      }
    } else if(length(h) == 2){
      if(!all(sign(sort(h)) == c(-1, 1))){
        stop("When specifying 2 control limits the two values should have reverse signs.")
      }
      h <- sort(h)
      ind_down <- which(chart$BK$val_down <= h[1])[1]
      ind_up <- which(chart$BK$val_up >= h[2])[1]
      ind <- tryCatch(min(ind_down, ind_up, na.rm = TRUE), warning = function(cond){NA})
      if(is.na(ind)){
        return(Inf)
      } else{
        return(chart$BK[ind,1] - chart$BK[1, "time"])
      }
    } else{
      stop("Please provide either 1 or 2 control limits.")
    }
  }

}

#' @describeIn runlength determines runlength of \code{bercusum} object
#' @export
runlength.bercusum <- function(chart, h, ...){
  if(missing(h)){
    stop("Please specify a control limit h.")
  } else if(!is.numeric(h)){
    stop("Control limit h must be numeric.")
  }
  if(missing(chart)){
    stop("Please provide a 'cgrcusum', 'bkcusum' or 'bercusum' chart as input.")
  }
  if(isFALSE(chart$call[["twosided"]]) | is.null(chart$call[["twosided"]])){
    if(!length(h) == 1){
      stop("Please provide only 1 control limit.")
    }
    ind <- which(abs(chart$CUSUM[,2]) >= abs(h))[1]
    if(is.na(ind)){
      return(Inf)
    } else{
      return(chart$CUSUM[ind,1] - chart$CUSUM[1, "time"])
    }
  } else{ #twosided Bernoulli-CUSUM
    if(length(h) == 1){
      if(h < 0){
        ind <- which(chart$CUSUM$val_down <= h)[1]
        if(is.na(ind)){
          return(Inf)
        } else{
          return(chart$CUSUM[ind,1] - chart$CUSUM[1, "time"])
        }
      } else if(h >= 0){
        ind <- which(chart$CUSUM$val_up >= h)[1]
        if(is.na(ind)){
          return(Inf)
        } else{
          return(chart$CUSUM[ind,1] - chart$CUSUM[1, "time"])
        }
      }
    } else if(length(h) == 2){
      if(!all(sign(sort(h)) == c(-1, 1))){
        stop("When specifying 2 control limits the two values should have reverse signs.")
      }
      h <- sort(h)
      ind_down <- which(chart$CUSUM$val_down <= h[1])[1]
      ind_up <- which(chart$CUSUM$val_up >= h[2])[1]
      ind <- tryCatch(min(ind_down, ind_up, na.rm = TRUE), warning = function(cond){NA})
      if(is.na(ind)){
        return(Inf)
      } else{
        return(chart$CUSUM[ind,1] - chart$CUSUM[1, "time"])
      }
    } else{
      stop("Please provide either 1 or 2 control limits.")
    }
  }
}




