#' @describeIn plot Plot a list of CUSUM charts
#' @import plotly
#' @export
plot.cusumlist <- function(x, unit_names, ...){
  n <- length(x)
  if(!missing(unit_names)){
    if(length(unit_names) != n){
      stop("Provided list of names must be equally long as the list of CUSUM charts.")
    }
  }
  plotframe <- data.frame(time = numeric(), value = numeric(), unit = factor(), type = factor())
  for(i in 1:n){
    chart <- x[[i]][[1]]
    chart <- chart[, c("time", "value")]
    if("h" %in% names(x[[i]])){
      chart$value <- chart$value/x[[i]]$h
    } else{
      chart$value <- chart$value/max(chart$value)
    }
    if(!missing(unit_names)){
      chart$unit <- rep(unit_names[i], nrow(chart))
    } else{
      chart$unit <- as.factor(rep(i, nrow(chart)))
    }
    chart$type <- rep(as.factor(class(x[[i]])), nrow(chart))
    plotframe <- rbind(plotframe, chart)
  }

  #Function for adding lines to plotly - not currently used
  hline <- function(y = 0, color = "black") {
    list(
      type = "line",
      x0 = min(plotframe$time),
      x1 = max(plotframe$time),
      xref = "paper",
      y0 = y,
      y1 = y,
      line = list(color = color)
    )

  }

  #g <- plotly::plot_ly(data = plotframe,  type = 'scatter', mode = 'lines')
  #We want to highlight 1 unit
  d <- highlight_key(plotframe, ~unit  )
  p <- ggplot( d, aes(x = time, y = value, group = unit)) + theme_bw() + geom_line() +geom_hline(yintercept = 1, color = "red")
  gg <- ggplotly( p, tooltip = c("unit", "time", "type")) # Add "value" here if you want to see value
  #add horizontal line
  plotly::layout(gg, showlegend = FALSE)
  return(highlight( gg, on = "plotly_hover", off = "plotly_deselect", color = "green" ))
}
