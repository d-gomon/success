#' @title Plot a list of CUSUM charts (interactive)
#'
#' @description Create an interactive plot visualizing a combination of control charts
#' which can be created using this package.
#'
#' @param x A list with each item containing one cumulative sum chart.
#' @param unit_names The unit names to be displayed in the interactive plot.
#' Must be of equal length as \code{x}.
#' @param scale Should charts be scaled with respect to their control limit/
#' maximum value? Default is \code{FALSE}.
#' @param highlight Should charts be highlighted on hover? Default is \code{FALSE}.
#' @param ... Further plotting parameters
#'
#'
#' @return An interactive plot will be produced in the current graphics device.
#' For more information on the possibilities for interaction, see \url{https://plotly.com/r/}.
#'
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom plotly highlight_key
#' @export
#'
#' @seealso \code{\link[success]{cgr_cusum}}, \code{\link[success]{bk_cusum}}, \code{\link[success]{bernoulli_cusum}}, \code{\link[success]{funnel_plot}}
#'
#' @examples
#' \dontrun{
#' require(survival)
#' #Extract data to construct CUSUM charts on
#' tdat <- subset(surgerydat, unit == 1 & entrytime < 365)
#' tdat2 <- subset(surgerydat, unit == 2 & entrytime < 365)
#'
#' #Determine model parameters
#' followup <- 100
#' tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
#' exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
#' tcoxmod <- coxph(exprfit, data= surgerydat)
#' exprfitber <- as.formula("(survtime <= followup) & (censorid == 1)~ age + sex + BMI")
#' glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
#'
#'
#' #Construct the charts
#' cgr <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE)
#' cgr$h <- 8.29
#' bk <- bk_cusum(data = tdat, theta = log(2), coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE)
#' bk$h <- 6.23
#' bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
#' followup = followup, theta = log(2))
#' bercus$h <- 3.36
#' bk2 <- bk_cusum(data = tdat2, theta = log(2), coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE)
#' bk2$h <- 6.23
#'
#' #Create the plot
#' interactive_plot(list(cgr, bk, bercus, bk2), unit_names =
#' c("hosp1", "hosp1", "hosp1", "hosp2"))
#' }



interactive_plot <- function(x, unit_names, scale = FALSE, highlight = FALSE, ...){
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
      a <- plot_ly(plotframe, x = ~time, y = ~value, type = 'scatter', mode = 'lines', split = ~unit, linetype = ~type, color = I("gray"))
      a <- plotly::layout(a, shapes = list(hline(1)))
      return(a)
    }
  } else{
    if(isTRUE(highlight)){
      a <- highlight_key(plotframe, key = ~unit)
      a <- plot_ly(a, x = ~time, y = ~value, type = 'scatter', mode = 'lines', split = ~unit, linetype = ~type, color = I("gray"))
      return(highlight(a, on = "plotly_hover", off = "plotly_deselect", color = "green"))
    } else{
      a <- plot_ly(plotframe, x = ~time, y = ~value, type = 'scatter', mode = 'lines', split = ~unit, linetype = ~type, color = I("gray"))
      return(a)
    }
  }

}
