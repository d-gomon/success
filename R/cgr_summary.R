#' @title Summarizes S3 objects in this package.
#'
#' @description Prints a summary of the
#'  \code{funnelplot} object.
#'
#'
#' @param object S3 object to summarize
#' @param ... extra parameters
#'
#' @return A \code{data.frame} with:
#' \itemize{
#' \item \code{unit}: unit number/identifier;
#' \item \code{observed}: the observed amount of failures at respective unit;
#' \item \code{expected}: the expected amount of failures at respective unit,
#' given that the unit is performing at target;
#' \item \code{numtotal}: total number of subjects at respective unit;
#' \item \code{p}: estimated probability of failure at unit;
#' \item \code{'0.xx'}: better/normal/worse proportion of failure at specified
#' confidence levels.
#' }
#'
#' @seealso \code{\link[success:funnel_plot]{funnel_plot}}
#'
#' @describeIn summary summarize instances detected by the
#' \code{funnelplot} object
#' @export
summary.funnelplot <- function(object, ...){
  # k <- object$conflev
  # outp <- list("call" = object$call)
  # for(i in k){
  #   temp <- sort(object$data[which(object$data[, as.character(i)] == "worse"), "unit"])
  #   outp[[as.character(i)]] = droplevels(as.factor(temp))
  # }
  # return(outp)
  return(object$data)
}
