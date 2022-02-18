#' @describeIn plot Plot a GBK-CUSUM
#' @import ggplot2
#' @noRd
plot.gbkcusum <- function(x, h = x$h, ...){
  times <- nrow(x$GBK)
  pltcusum <- data.frame(time = numeric(), value = numeric())
  for (i in 1:times){
    pltcusum <- rbind(pltcusum, c(x$GBK[i, "time"], x$GBK[i, "val_down"]))
    pltcusum <- rbind(pltcusum, c(x$GBK[i, "time"], x$GBK[i, "val_up"]))
  }
  colnames(pltcusum) = c("time", "value")
  g <- ggplot(data = pltcusum, mapping= aes(x = time, y = value)) + geom_line() + geom_hline(yintercept = h, color = "red")
  return(g)
}
