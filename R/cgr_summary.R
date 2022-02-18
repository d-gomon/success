#' @title Summarizes (or elaborates on) S3 objects in this package.
#'
#' @description Prints the (name of the) instances performing worse than expected
#' in a \code{funnelplot} object at the confidence levels specified in the
#' \code{funnelplot} object.
#'
#'
#' @param object S3 object to summarize
#' @param ... extra parameters
#'
#' @return A list with:
#' \itemize{
#' \item \code{call}: the call used to obtain the input object,
#' \item \code{'0.xx'}: the detected instances at specified confidence level(s).
#' }
#'
#' @seealso \code{\link[success:funnel_plot]{funnel_plot}}
#'
#' @describeIn summary summarize instances detected by the
#' \code{funnelplot} object
#' @export
summary.funnelplot <- function(object, ...){
  k <- object$conflev
  outp <- list("call" = object$call)
  for(i in k){
    temp <- sort(object$data[which(object$data[, as.character(i)] == "worse"), "unit"])
    outp[[as.character(i)]] = droplevels(as.factor(temp))
  }
  return(outp)
}
