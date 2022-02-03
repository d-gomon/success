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
#' bercus <- bernoulli_cusum(data = subset(surgerydat, hosp_num == 14), glmmod = glmmodber,
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
  }
  if(missing(chart)){
    stop("Please provide a 'cgrcusum', 'bkcusum' or 'bercusum' chart as input.")
  }
  ind <- which(chart$CGR[,2] >= h)[1]
  if(is.na(ind)){
    return(Inf)
  } else{
    return(chart$CGR[ind,1])
  }
}





#' @describeIn runlength determines runlength of \code{bkcusum} object
#' @export
runlength.bkcusum <- function(chart, h, ...){
  if(missing(h)){
    stop("Please specify a control limit h.")
  }
  if(missing(chart)){
    stop("Please provide a 'cgrcusum', 'bkcusum' or 'bercusum' chart as input.")
  }
  ind <- which(chart$BK[,2] >= h)[1]
  if(is.na(ind)){
    return(Inf)
  } else{
    return(chart$BK[ind,1])
  }
}

#' @describeIn runlength determines runlength of \code{bercusum} object
#' @export
runlength.bercusum <- function(chart, h, ...){
  if(missing(h)){
    stop("Please specify a control limit h.")
  }
  if(missing(chart)){
    stop("Please provide a 'cgrcusum', 'bkcusum' or 'bercusum' chart as input.")
  }
  ind <- which(chart$CUSUM[,2] >= h)[1]
  if(is.na(ind)){
    return(Inf)
  } else{
    return(chart$CUSUM[ind,1])
  }
}




