#' Risk-adjusted funnel plot
#'
#' @description This function allows to construct a risk-adjusted funnel plot
#' for comparing survival proportion between units, see Spiegelhalter (2005).
#'
#' @param data A \code{data.frame} with rows representing subjects and the
#' following named columns: \describe{
#'   \item{\code{entrytime}:}{time of entry into study (numeric);}
#'   \item{\code{survtime}:}{time from entry until event (numeric);}
#'   \item{\code{censorid}:}{censoring indicator (0 = right censored, 1 = observed),
#'    (integer);}
#'   \item{\code{unit}:}{integer or character indicating which unit
#'   (f.e. hospital) the observation belongs to.}
#' } and optionally additional covariates used for risk-adjustment.
#' @param ctime Construction time at which the funnel plot
#' should be determined. Maximum possible time used when not specified.
#' @param p0 The baseline failure probability at \code{entrytime + followup} for individuals.
#' If not specified, average failure proportion over whole data is used instead.
#' @param glmmod A generalized linear regression model as produced by
#' the function \code{\link[stats:glm]{glm()}}. Recommended: \cr
#' \code{glm(as.formula("(survtime <= followup) & (censorid == 1) ~ covariates"), data = data)}. \cr
#' Alternatively, a list with the following elements:
#' \describe{
#'   \item{\code{formula}:}{a \code{\link[stats:formula]{formula()}} in the form \code{~ covariates};}
#'   \item{\code{coefficients}:}{a named vector specifying risk adjustment coefficients
#'   for covariates. Names must be the same as in \code{formula} and colnames of \code{data}.}
#' }
#' @param followup The followup time for every individual. At what time
#' after subject entry do we consider the outcome?
#' @param predlim A vector of confidence levels for the prediction limits of interest. Default is c(0.95, 0.99).
#' @param assist (optional): Output of the function \code{\link[success:parameter_assist]{parameter_assist()}}
#'
#' @return An object of class "funnelplot" containing:
#' \itemize{
#' \item \code{data}: A \code{data.frame} containing:
#' \describe{
#'   \item{\code{unit}:}{unit number/name;}
#'   \item{\code{observed}:}{observed number of failures at unit;}
#'   \item{\code{expected}:}{expected (risk-adjusted) number of failures at unit;}
#'   \item{\code{numtotal}}{total number of individuals considered at this unit;}
#'   \item{\code{p}:}{(risk-adjusted) proportion of failure at unit;}
#'   \item{\code{predlimels}:}{worse/in-control/better performance than expected at
#'   specified confidence levels.}
#' }
#' \item \code{call}: the call used to obtain output
#' \item \code{plotdata}: data used for plotting confidence intervals
#' \item \code{predlim}: specified confidence level(s)
#' \item \code{p0}: (Estimated) baseline failure probability
#' }
#'
#' @references Spiegelhalter D. J. (2005). Funnel plots for comparing
#' institutional performance. Statistics in medicine, 24(8), 1185-1202.
#' \doi{10.1002/sim.1970}
#'
#' @importFrom stats predict.glm
#' @importFrom stats model.matrix
#' @importFrom stats qnorm
#' @export
#'
#' @author Daniel Gomon
#' @family quality control charts
#' @seealso \code{\link[success]{plot.funnelplot}}, \code{\link[success]{summary.funnelplot}}
#'
#'
#' @examples
#' #Determine a risk-adjustment model using a generalized linear model.
#' #Outcome (survival in first 100 days) is regressed on the available covariates:
#' exprfitfunnel <- as.formula("(survtime <= 100) & (censorid == 1)~ age + sex + BMI")
#' glmmodfun <- glm(exprfitfunnel, data = surgerydat, family = binomial(link = "logit"))
#' #Determine the necessary values to produce a funnel plot
#' funnel <- funnel_plot(data = surgerydat, ctime = 3*365, glmmod = glmmodfun, followup = 100)
#' #Produce a funnel plot!
#' plot(funnel)
#' \dontrun{
#' require(plotly)
#' #Create an interactive plot!
#' ggplotly(plot(funnel))
#' }






funnel_plot <- function(data, ctime, p0, glmmod, followup, predlim = c(0.95, 0.99),
                        assist){
  entrytime <- unit <- NULL
  call <- match.call()
  if(!missing(assist)){
    list2env(assist, envir = environment())
    data <- assist$baseline_data
  }


  #Check that predlim is numerical vector with values between 0 and 1
  if(!all(is.numeric(predlim), all(predlim > 0), all(predlim < 1))){
    stop("Argument predlim must be numeric vector with values between 0 and 1.")
  }
  #Sort predlim
  predlim <- sort(predlim)

  #Check that followup is a numeric value greater than 0
  if(!all(is.numeric(followup), length(followup) == 1, followup > 0)){
    stop("Argument followup must be a single numeric value larger than 0.")
  }

  #Basic data checks (global for BK, CGR and Bernoulli & funnel)
  if(missing(data)){
    stop("Please provide data to construct chart.")
  } else{
    data <- check_data(data)
  }

  #First perform a logistic regression to obtain coefficients for the Risk-adjustment, then specify them later
  #dat has to be a dataframe containing at least entrytime (time of entry), survtime (time until failure), and additionally the    covariates to RA   on, it also has to contain $unit indicating the hospital in question
  #glmmodel is the risk-adjustment model, either an object of class "glm" or $formula and $coefficients
  #followup is time until which we consider outcomes, usually 365 (1 year) as we consider 1 year post transplant
  #Specify institute name or number in dat$unit
  #predlim indicates the confidence levels at which to plot the boundaries
  #time is the chronological time at which the FUNNEL chart should be constructed, we remove non-qualifying cases
  if(!missing(ctime)){
    newdata <- subset(data, entrytime + followup <= ctime)
  } else{
    newdata <- data
  }
  if(missing(p0)){
    warning("No value provided for null (hypothesis) failure probability. Determining using average over whole data set.", immediate. = TRUE)
    p0 <- length(which((newdata$survtime <= followup) & (newdata$censorid == 1)))/nrow(newdata)
  }
  plotframe <- data.frame(unit = character(), observed = double(), expected = double(), numcases = double())
  for(j in unique(newdata$unit)){
    tempdata <- subset(newdata, unit == j)
    tempnum <- length(tempdata$survtime)
    if(!missing(glmmod)){
      if(inherits(glmmod, "glm")){
        tempprobs <-  predict(glmmod, newdata = tempdata, type = "response")
        tempexpec <- sum(tempprobs)
      } else{
        mmatrix <- model.matrix(glmmod$formula, tempdata)
        coeffs <- glmmod$coefficients[colnames(mmatrix)]
        tempexpec <- sum(c(1/(1 + exp(-mmatrix %*% coeffs))))
      }
    } else{
      tempexpec <- tempnum * p0
    }
    tempexpec <- ifelse(tempexpec > nrow(tempdata), nrow(tempdata), tempexpec)
    tempobs <- length(which((tempdata$survtime <= followup) & (tempdata$censorid == 1)))
    temprow <- data.frame(as.character(j), tempobs, tempexpec,  tempnum)
    plotframe <- rbind(plotframe, temprow)
  }
  colnames(plotframe) = c("unit", "observed", "expected", "numtotal")

  plotframe$p <- plotframe$observed/plotframe$expected * p0
  plotframe$p <- ifelse(plotframe$p > 1, 1, plotframe$p)
  boundplotframe <- data.frame(number = double(),predlim = double(),lower = double(), upper = double())
  plotseq <- seq(max(1, min(plotframe$numtotal)-10), max(plotframe$numtotal) +10, by = 1)
  findbounds <- function(t, predlim){
    return(c(p0 - qnorm(predlim) * sqrt((p0*(1-p0))/t),p0 + qnorm(predlim) * sqrt((p0*(1-p0))/t)))
  }
  for(k in 1:length(predlim)){
    temprow2 <- data.frame(plotseq, rep(predlim[k], length(plotseq)), t(sapply(plotseq, function(t) findbounds(t, predlim[k]))))
    boundplotframe <- rbind(boundplotframe, temprow2)
    tempchar <- character(length = nrow(plotframe))
    for(i in 1:nrow(plotframe)){
      tempbounds <- findbounds(plotframe$numtotal[i], predlim = predlim[k])
      if(plotframe$p[i] > tempbounds[2]){
        tempchar[i] <- "worse"
      } else if(plotframe$p[i] < tempbounds[1]){
        tempchar[i] <- "better"
      } else{
        tempchar[i] <- "in-control"
      }
    }
    plotframe <- cbind(plotframe, tempchar)
    colnames(plotframe)[ncol(plotframe)] <- as.character(predlim[k])
  }
  colnames(boundplotframe) = c("numtotal", "predlim", "lower", "upper")
  funnelp <- list(data = plotframe,
              call = call,
              plotdata = boundplotframe,
              predlim = predlim,
              p0 = p0)
  class(funnelp) <- "funnelplot"
  funnelp
}







