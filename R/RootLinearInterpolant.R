#' @noRd
#' @keywords internal
#'
#' @details Code for this function was obtained from stackoverflow.
#' The original topic can be found at
#' \url{https://stackoverflow.com/questions/52650467/how-to-estimate-x-value-from-y-value-input-after-approxfun-in-r}.
#' To the best of our knowledge, the author of the code is Zheyuan Li
#' \url{https://stackoverflow.com/users/4891738/zheyuan-li}.
#'
#' @author Zheyuan Li
#'
#' @param x x values of original data
#' @param y y values of original data
#' @param y0 point at which to find the inverse of linear interpolant of x and y
#'
RootLinearInterpolant <- function (x, y, y0 = 0) {
  if (is.unsorted(x)) {
    ind <- order(x)
    x <- x[ind]; y <- y[ind]
  }
  z <- y - y0
  ## which piecewise linear segment crosses zero?
  k <- which(z[-1] * z[-length(z)] < 0)
  if(length(k) == 0){
    return(NA)
  }
  ## analytically root finding
  xk <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
  xk
}
