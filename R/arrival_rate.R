#' @title Estimate arrival rate of a Poisson points process
#'
#' @description In a Poisson point process, subjects arrive with exponentially
#' distributed inter-arrival times with rate \eqn{\psi}{\psi}. This function
#' can be used to estimate the parameter \eqn{\psi}{\psi}.
#'
#' @param data A \code{data.frame} containing the following named column for each subject:
#' \describe{
#'   \item{\code{entrytime}:}{time of entry into study (numeric);}
#' } If the \code{data.frame} also contains a column named \code{unit}, the arrival
#' rate will be determined for each unit separately.
#'
#'
#'
#' @export
#'
#' @return A (named) vector containing the estimated arrival rate in the data,
#' or for each unit in the data.
#'
#' @author Daniel Gomon
#'
#' @examples
#' arrival_rate(surgerydat)
#'
#'






arrival_rate <- function(data){
  unit <- NULL
  #Transform to df if vector is given
  if(is.vector(data)){
    warning("Provided vector is not named, assuming that specified data represents entrytimes.")
    data <- data.frame(entrytime = data)
  }
  #Check whether entrytime is in column names
  if(!all("entrytime" %in% colnames(data), is.numeric(data$"entrytime"))){
    stop("Argument data should contain at least a numeric column named 'entrytime'")
  }
  get_arr_rate <- function(dat){
    (nrow(dat)-1)/(max(dat$"entrytime") - min(dat$"entrytime"))
  }
  if("unit" %in% colnames(data)){
    units <- unique(data$"unit")
    arr_rates <- numeric(length(units))
    for(i in seq_along(units)){
      tdat <- subset(data, unit == units[i])
      arr_rates[i] <- get_arr_rate(tdat)
    }
    names(arr_rates) <- units
  } else{
    arr_rates <- get_arr_rate(data)
  }
  return(sort(arr_rates))
}
