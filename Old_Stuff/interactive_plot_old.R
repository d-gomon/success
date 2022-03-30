interactive_plot_old <- function(x, unit_names, scale = FALSE, highlight = FALSE, ...){
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
    if(isTRUE(scale)){
      if("h" %in% names(x[[i]])){
        chart$value <- chart$value/x[[i]]$h
      } else{
        chart$value <- chart$value/max(chart$value)
      }
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
      x0 = 0,
      x1 = 1,
      xref = "paper",
      y0 = y,
      y1 = y,
      line = list(color = color)
    )

  }
  if(isTRUE(scale)){
    if(isTRUE(highlight)){
      a <- highlight_key(plotframe, key = ~unit)
      a <- plot_ly(a, x = ~time, y = ~value, type = 'scatter', mode = 'lines', split = ~unit, linetype = ~type, color = I("gray"))
      a <- plotly::layout(a, shapes = list(hline(1)))
      return(highlight(a, on = "plotly_hover", off = "plotly_deselect", color = "green"))
    } else{
      a <- plot_ly(plotframe, x = ~time, y = ~value, type = 'scatter', mode = 'lines', linetype = ~type, color = ~unit, colors = "Set2")
      a <- plotly::layout(a, shapes = list(hline(1)))
      return(a)
    }
  } else{
    if(isTRUE(highlight)){
      a <- highlight_key(plotframe, key = ~unit)
      a <- plot_ly(a, x = ~time, y = ~value, type = 'scatter', mode = 'lines', split = ~unit, linetype = ~type, color = I("gray"))
      return(highlight(a, on = "plotly_hover", off = "plotly_deselect", color = "green"))
    } else{
      a <- plot_ly(plotframe, x = ~time, y = ~value, type = 'scatter', mode = 'lines', linetype = ~type, color = ~unit, colors = "Set2")
      return(a)
    }
  }

}
