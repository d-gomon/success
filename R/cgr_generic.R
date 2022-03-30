#' @title Plot a quality control chart
#'
#' @description Plot a \code{cgrcusum}, \code{bkcusum},
#'  \code{bercusum} or \code{funnelplot} chart, or a list containing a combination of
#'  \code{'bercusum'}, \code{'bkcusum'} and \code{'cgrcusum'} charts.
#'
#'
#' @param x Chart to plot
#' @param h Control limit to display for \code{cgrcusum}, \code{bkcusum} or \code{bercusum}
#' @param percentage Should output be shown in percentages? Default is \code{TRUE}.
#' @param ... Further plotting parameters
#'
#'
#' @author Daniel Gomon
#'
#' @seealso \code{\link[success]{cgr_cusum}}, \code{\link[success]{bk_cusum}}, \code{\link[success]{bernoulli_cusum}}, \code{\link[success]{funnel_plot}}
#'
#' @returns A plot of the associated chart is displayed in the current graphics device.
#'
#' @describeIn plot Plot a CGR-CUSUM
#' @import ggplot2
#' @export
plot.cgrcusum <- function(x, h, ...){
  time <- value <- NULL
  requireNamespace('ggplot2')
  g <- ggplot(as.data.frame(x$CGR),
                       mapping = aes(x = time, y = value)) + geom_line()
  if(!missing(h)){
    g <- g + geom_hline(yintercept = h, color = "red")
  } else if("h" %in% names(x)){
    g <- g + geom_hline(yintercept = x$h, color = "red")
  }
  g <- g + labs(x = "Time", y = "Value")
  return(g)
}

#' @describeIn plot Plot a BK-CUSUM
#' @import ggplot2
#' @export
plot.bkcusum <- function(x, h, ...){
  time <- value <- val_up <- val_down <- NULL
  requireNamespace('ggplot2')
  if(isFALSE(x$call[["twosided"]]) | is.null(x$call[["twosided"]])){
    g <- ggplot(as.data.frame(x$BK),
                mapping = aes(x = time, y = value)) + geom_line()
    if(!missing(h)){
      g <- g + geom_hline(yintercept = h, color = "red")
    } else if("h" %in% names(x)){
      g <- g + geom_hline(yintercept = x$h, color = "red")
    }
  }else if(isTRUE(x$call[["twosided"]])){
    g <- ggplot()
    g <- g + geom_line(as.data.frame(x$BK),
                mapping = aes(x = time, y = val_up), color = "black")
    g <- g + geom_line(as.data.frame(x$BK),
                       mapping = aes(x = time, y = val_down), color = "blue")
    if(!missing(h)){
      g <- g + geom_hline(yintercept = h, color = "red")
    } else if("h" %in% names(x)){
      g <- g + geom_hline(yintercept = x$h, color = "red")
    }
  }
  g <- g + labs(x = "Time", y = "Value")
  return(g)
}



#' @describeIn plot Display a funnel plot
#' @import ggplot2
#' @importFrom grDevices palette
#' @importFrom grDevices palette.colors
#' @export
plot.funnelplot <- function(x, percentage = TRUE, ...){
  numtotal <- lower <- conflev <- upper <- p <- unit <- detection <- NULL

  #Supply plot.FUNNEL with output from FUNNELsim or a data frame with $unitdata and $p0 and $conflev
  if(percentage == TRUE){
    x$plotdata[, c("lower", "upper")] <- x$plotdata[, c("lower", "upper")] * 100
    x$data$p <- x$data$p * 100
    x$p0 <- x$p0*100
  }
  if(percentage == TRUE){
    maxy <- min(100,max(x$data$p) + 0.1*max(x$data$p))
    miny <- max(0, min(x$data$p) - 0.1*max(x$data$p))
  } else{
    maxy <- min(1,max(x$data$p) + 0.1*max(x$data$p))
    miny <- max(0, min(x$data$p) - 0.1*max(x$data$p))
  }

  #Determine the number of required colours
  numcolours <- length(x$conflev) + 1
  cols <- palette.colors(n = numcolours)
  #Determine which colour to use for point:
  cols_columns <- ncol(x$data) - 5
  finalcols <- rep("in-control", length = nrow(x$data))
  for (k in rev(1:cols_columns)){
    t_col_data <- x$data[,ncol(x$data) - (k-1)]
    finalcols[which(t_col_data == "worse" | t_col_data == "better")] <- rev(x$conflev)[k]
  }
  x$data$detection <- finalcols
  t <- ggplot() + geom_line(data = x$plotdata, mapping= aes(x = numtotal, y = lower, group = as.factor(conflev)),
                            colour = "blue", linetype = "dashed") +
  geom_line(data = x$plotdata,aes(x = numtotal, y = upper, group = as.factor(conflev)),
            colour = "blue", linetype = "dashed") +
  geom_point(data = x$data, mapping= aes(x = numtotal, y = p, colour = detection, group = unit)) +
    #theme(legend.position = "none") +
    geom_hline(yintercept = x$p0, colour = "grey") + ylim(miny, maxy)
  t <- t + labs(x = "Number of outcomes", y = paste0("(Risk-adjusted) Proportion of failure (%)"))
  return(t)
}


#' @describeIn plot Plot a Bernoulli CUSUM
#' @import ggplot2
#' @export
plot.bercusum <- function(x, h = x$h, ...){
  time <- value <- NULL
  g <- ggplot(as.data.frame(x$CUSUM), mapping = aes(x = time, y = value)) + geom_line() + geom_hline(yintercept = h, color = "red")
  return(g)
}














