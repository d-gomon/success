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
#' @param unit_label Should unit labels be displayed next to detected units in the funnel plot?
#' Default is \code{TRUE}.
#' @param label_size Size of the labels when \code{unit_label} is \code{TRUE}. Default is 3.
#' @param col_fill Single fill colour for the prediction intervals in the funnel plot.
#' In any format that \code{\link[grDevices]{col2rgb}} accepts. Default is \code{"blue"}.
#' @param ... Further plotting parameters
#'
#'
#' @author Daniel Gomon
#'
#' @seealso \code{\link[success]{cgr_cusum}}, \code{\link[success]{bk_cusum}}, \code{\link[success]{bernoulli_cusum}}, \code{\link[success]{funnel_plot}}
#'
#' @returns A plot of the associated chart is displayed in the current graphics device.
#' @importFrom ggplot2 ggplot aes geom_line geom_hline labs
#' @describeIn plot Plot a CGR-CUSUM
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
#' @importFrom ggplot2 scale_colour_manual geom_point ylim scale_fill_manual geom_ribbon
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats relevel
#' @importFrom grDevices palette palette.colors colorRamp adjustcolor
#' @export
plot.funnelplot <- function(x, percentage = TRUE, unit_label = TRUE,
                            label_size = 3, col_fill = "blue", ...){
  numtotal <- lower <- predlim <- upper <- p <- unit <- detection <- NULL

  #Supply plot.FUNNEL with output from FUNNELsim or a data frame with $unitdata and $p0 and $predlim
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


  #X-axis for ribbons
  plotseq <- seq(max(1, min(x$data$numtotal)-10), max(x$data$numtotal) +10, by = 1)

  #Determine the number of required colours for the points!
  numcolours <- length(x$predlim) + 1
  #Gradient from green to red displaying detection levels
  cols <- grDevices::colorRampPalette(c("green", "red"))(numcolours)
  #Determine which colour to use for point:
  names(cols) = c("in-control", sort(as.numeric(x$predlim)))
  #Create color scale for ggplot
  colScale <- scale_colour_manual(name = "Detection", values = cols)

  #Now we create a colScale for the fill of the prediction intervals.
  colsfill <- sapply(seq(0.5, 0.1, length.out = numcolours-1), function(x) adjustcolor(col_fill, alpha.f = x))
  names(colsfill) = c(sort(as.numeric(x$predlim)))
  colScaleFill <- scale_fill_manual(name = "Prediction interval", values = colsfill)

  #Determine highest detection level.
  finalcols <- rep("in-control", length = nrow(x$data))
  for (k in rev(1:(numcolours-1))){
    t_col_data <- x$data[,ncol(x$data) - (k-1)]
    finalcols[which(t_col_data == "worse" | t_col_data == "better")] <- rev(x$predlim)[k]
  }
  finalcols <- as.factor(finalcols)
  finalcols <- stats::relevel(finalcols, "in-control")
  x$data$detection <- finalcols
  t <- ggplot(x$data, mapping = aes(x = numtotal, y = p)) +
    #Old code: lines used to display prediction intervals.
  #geom_line(data = x$plotdata, mapping= aes(x = numtotal, y = lower, group = as.factor(predlim)),
  #                          colour = "blue", linetype = "dashed") +
  #geom_line(data = x$plotdata, aes(x = numtotal, y = upper, group = as.factor(predlim)),
  #          colour = "blue", linetype = "dashed") +
  geom_ribbon(data = x$plotdata, aes(x = numtotal, ymin = lower, ymax = upper, group = as.factor(predlim), fill = as.factor(predlim)), outline.type = "both") +
  colScaleFill +
  geom_point(data = x$data, mapping= aes(x = numtotal, y = p, colour = as.factor(detection), group = unit)) +
    #theme(legend.position = "none") +
    geom_hline(yintercept = x$p0, colour = "grey") + ylim(miny, maxy) + colScale
  if(isTRUE(unit_label)){
    t <- t + geom_label_repel(aes(color = factor(detection), label = ifelse(detection!="in-control", unit, "")), size = label_size, show.legend = FALSE)
  }
  t <- t + labs(x = "Number of outcomes", y = paste0("(Risk-adjusted) Proportion of failure (%)"))
  return(t)
}


#' @describeIn plot Plot a Bernoulli CUSUM
#' @export
plot.bercusum <- function(x, h = x$h, ...){
  time <- value <- val_up <- val_down <- NULL
  requireNamespace('ggplot2')
  if(isFALSE(x$call[["twosided"]]) | is.null(x$call[["twosided"]])){
    g <- ggplot(as.data.frame(x$CUSUM),
                mapping = aes(x = time, y = value)) + geom_line()
    if(!missing(h)){
      g <- g + geom_hline(yintercept = h, color = "red")
    } else if("h" %in% names(x)){
      g <- g + geom_hline(yintercept = x$h, color = "red")
    }
  }else if(isTRUE(x$call[["twosided"]])){
    g <- ggplot()
    g <- g + geom_line(as.data.frame(x$CUSUM),
                       mapping = aes(x = time, y = val_up), color = "black")
    g <- g + geom_line(as.data.frame(x$CUSUM),
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














