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
#' @param group_by Character indicating how to group the CUSUM charts in the plot.
#' Possible values are \code{c("none", "unit", "type")}. Default is \code{"none"}.
#' @param manual_colors A character vector specifying which colors to use
#' for the units in the data. By default, the "Dark2" color set from
#' \code{\link[RColorBrewer:brewer.pal]{brewer.pal()}} will be used.
#'
#'
#' @return An interactive plot will be produced in the current graphics device.
#' For more information on the possibilities for interaction, see \url{https://plotly.com/r/}.
#'
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom plotly highlight_key
#' @importFrom plotly add_trace
#' @importFrom plotly add_lines
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @seealso \code{\link[success]{cgr_cusum}}, \code{\link[success]{bk_cusum}}, \code{\link[success]{bernoulli_cusum}}, \code{\link[success]{funnel_plot}}
#'
#' @examples
#' \donttest{
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


interactive_plot <- function(x, unit_names, scale = FALSE,
                             group_by = c("none", "unit", "type"), highlight = FALSE,
                             manual_colors = c(), ...){
  id <- NULL
  #https://plotly-r.com/index.html
  n <- length(x)
  if(!missing(unit_names)){
    if(length(unit_names) != n){
      stop("Provided list of names must be equally long as the list of CUSUM charts.")
    }
  } else{
    unit_names <- as.factor(1:n)
  }

  if(length(group_by > 1)){
    group_by <- group_by[1]
  }

  if(isTRUE(scale)){
    for(i in 1:n){
      if("h" %in% names(x[[i]])){
        x[[i]][[1]]$value <- x[[i]][[1]]$value/x[[i]]$h
      } else{
        x[[i]][[1]]$value <- x[[i]][[1]]$value/max(x[[i]][[1]]$value)
      }
    }
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


  ########################HIGHLIGHT############################
  if(isTRUE(highlight)){
    plotframe <- data.frame(time = numeric(), value = numeric(), unit = factor(), type = factor())
    for(i in 1:n){
      chart <- x[[i]][[1]]
      chart <- chart[, c("time", "value")]
      chart$unit <- rep(unit_names[i], nrow(chart))
      chart$type <- rep(as.factor(class(x[[i]])), nrow(chart))
      chart$id <- paste(chart$unit, chart$type)
      plotframe <- rbind(plotframe, chart)
    }
    tx <- highlight_key(plotframe, key = ~id)
    a <- group_by(plot_ly(tx, color = I("black")),id)
    a <- add_lines(group_by(a, id), x = ~ time, y = ~value)
    if(isTRUE(scale)){
      a <- plotly::layout(a, shapes = list(hline(1)))
    }
    return(highlight(a, on = "plotly_click", off = "plotly_deselect", selectize = TRUE, dynamic = TRUE, persistent = TRUE))
  }


  #######################NO HIGHLIGHT############################
  #Line types
  #cusumlty <- c("bercusum" = 3, "bkcusum" = 5, "cgrcusum" = 1)
  cusumlty <- c("bercusum" = "dot", "bkcusum" = "dash", "cgrcusum" = "solid")

  #Colours for different units
  ncols <- length(unique(unit_names))
  if(!missing(manual_colors)){
    if(length(manual_colors) != ncols){
      manual_colors <- rep(manual_colors, ceiling(ncols/length(manual_colors)))
    }
    unitcols <- manual_colors[1:ncols]
  } else{
    unitcols <- colorRampPalette(brewer.pal(n = 8, "Set2"))(ncols)
  }

  names(unitcols) <- unique(unit_names)
  for(i in 1:n){
    x[[i]][[1]]$col <- unit_names[i]
  }

  a <- plot_ly( )
  for(i in 1:n){
    #group_by checking
    if(group_by == "none"){
      tlegendgroup = ""
      tlegendgrouptitle = ""
      tname = paste(unit_names[i], class(x[[i]]), sep = "\n")
    } else if(group_by == "unit"){ #Grouping by unit
      tlegendgroup = unit_names[i]
      tlegendgrouptitle = list(text = unit_names[i])
      tname = class(x[[i]])
    } else if(group_by == "type"){ #Grouping by chart type
      tlegendgroup = class(x[[i]])
      tlegendgrouptitle = list(text = class(x[[i]]))
      tname = unit_names[i]
    }
    a <- add_lines(a, data = x[[i]][[1]], x = ~time, y = ~value,
                   colors = unitcols, color = ~col, line = list(dash = cusumlty[class(x[[i]])]),
                  mode = 'lines', name = tname, legendgroup = tlegendgroup,
                  legendgrouptitle = tlegendgrouptitle)
    #, text = class(x[[i]]) -> can add extra info into hoverbox   linetype = cusumlty[class(x[[i]])],
    #https://plotly.com/r/reference/scatter/#scatter-hovertemplate
  }
  if(isTRUE(scale)){
    a <- plotly::layout(a, shapes = list(hline(1)))
  }
  return(a)
}














