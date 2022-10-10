#' @title Risk-adjusted Bernoulli CUSUM
#'
#' @description This function can be used to construct a risk-adjusted Bernoulli
#' CUSUM chart for survival data.
#' It requires the specification of one of the following combinations of parameters
#' as arguments to the function:
#' \itemize{
#' \item \code{glmmod} & \code{theta}
#' \item \code{p0} & \code{theta}
#' \item \code{p0} & \code{p1}
#' }
#'
#'
#' @param data A \code{data.frame} containing the following named columns for each subject:
#' \describe{
#'   \item{\code{entrytime}:}{time of entry into study (numeric);}
#'   \item{\code{survtime}:}{time from entry until event (numeric);}
#'   \item{\code{censorid}:}{censoring indicator (0 = right censored, 1 = observed)
#'    (integer);}
#' } and optionally additional covariates used for risk-adjustment.
#' @param followup The value of the follow-up time to be used to determine event time.
#' Event time will be equal to \code{entrytime + followup} for each subject.
#' @param glmmod Generalized linear regression model used for risk-adjustment as produced by
#' the function \code{\link[stats:glm]{glm()}}. Suggested: \cr
#' \code{glm(as.formula("(survtime <= followup) & (censorid == 1) ~ covariates"), data = data)}. \cr
#' Alternatively, a list containing the following elements:
#' \describe{
#'   \item{\code{formula}:}{a \code{\link[stats:formula]{formula()}} in the form \code{~ covariates};}
#'   \item{\code{coefficients}:}{a named vector specifying risk adjustment coefficients
#'   for covariates. Names must be the same as in \code{formula} and colnames of \code{data}.}
#' }
#' @param theta The \eqn{\theta}{\theta} value used to specify the odds ratio
#'  \eqn{e^\theta}{e^\theta} under the alternative hypothesis.
#'  If \eqn{\theta >= 0}{\theta >= 0}, the chart will try to detect an increase
#'   in hazard ratio (upper one-sided). If \eqn{\theta < 0}{\theta < 0},
#' the chart will look for a decrease in hazard ratio (lower one-sided).
#' Note that \deqn{p_1 = \frac{p_0 e^\theta}{1-p_0 +p_0 e^\theta}.}{p1 = (p0 * e^\theta)/(1-p0+p0 * e^\theta).}
#' @param p0 The baseline failure probability at \code{entrytime + followup} for individuals.
#' @param p1 The alternative hypothesis failure probability at \code{entrytime + followup} for individuals.
#' @param h (optional): Control limit to be used for the procedure.
#' @param stoptime (optional): Time after which the value of the chart should no longer be determined.
#' @param assist (optional): Output of the function \code{\link[success:parameter_assist]{parameter_assist()}}
#' @param twosided (optional): Should a two-sided Bernoulli CUSUM be constructed?
#' Default is \code{FALSE}.
#'
#'
#' @details The Bernoulli CUSUM chart is given by
#' \deqn{S_n = \max(0, S_{n-1} + W_n),}{S_n = max(0, S_{n-1} + W_n),} where
#' \deqn{W_n = X_n \ln \left( \frac{p_1 (1-p_0)}{p_0(1-p_1)}  \right) + \ln \left( \frac{1-p_1}{1-p_0} \right)}{W_n = X_n ln((p_1 * (1-p_0))/(p_0 * (1-p_1))) + ln((1-p_1)/(1-p_0))}
#' and \eqn{X_n}{X_n} is the outcome of the \eqn{n}{n}-th (chronological) subject in the data.
#' Instead of the standard practice of displaying patient numbering on the
#' x-axis, the time of outcome is displayed.
#'
#'
#' @return An object of class \code{bercusum} containing:
#' \itemize{
#' \item \code{CUSUM}: A \code{data.frame} containing the following named columns:
#' \describe{
#'   \item{\code{time}:}{times at which chart is constructed;}
#'   \item{\code{value}:}{value of the chart at corresponding times;}
#'   \item{\code{numobs}:}{number of observations at corresponding times.}
#' }
#' \item \code{call}: the call used to obtain output;
#' \item \code{glmmod}: coefficients of the \code{\link[stats:glm]{glm()}} used
#'  for risk-adjustment, if specified;
#' \item \code{stopind}: indicator for whether the chart was stopped by the
#' control limit.
#' }
#'
# There are \code{\link[cgrcusum:plot.bercusum]{plot}} and
#   \code{\link[cgrcusum:runlength.bercusum]{runlength}} methods for "bercusum" objects.

#'
#' @importFrom stats predict.glm
#' @importFrom stats model.matrix
#' @export
#'
#' @author Daniel Gomon
#'
#' @seealso \code{\link[success]{plot.bercusum}}, \code{\link[success]{runlength.bercusum}}
#'
#'
#' @examples
#' #We consider patient outcomes 100 days after their entry into the study.
#' followup <- 100
#' #Determine a risk-adjustment model using a generalized linear model.
#' #Outcome (failure within 100 days) is regressed on the available covariates:
#' exprfitber <- as.formula("(survtime <= followup) & (censorid == 1)~ age + sex + BMI")
#' glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
#' #Construct the Bernoulli CUSUM on the 1st hospital in the data set.
#' bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
#'  followup = followup, theta = log(2))
#' #Plot the Bernoulli CUSUM
#' plot(bercus)




bernoulli_cusum <- function(data, followup, glmmod, theta, p0, p1, h, stoptime,
                   assist, twosided = FALSE){
  entrytime <- otime <- NULL

  if(!missing(assist)){
    list2env(assist, envir = environment())
  }
  call <- match.call()
  #exp(theta) is the Odds Ratio under the alternative hypothesis
  #Supply either of the following combinations:
  #1. glmmodel + theta  2. p0 + theta   3. p0 + p1
  #Relationship between p1 and theta: p1 = (p0*exp(theta))/((1-p0)*(1+p0*exp(theta)))
  #First perform a logistic regression model on in-control data to obtain RA probs
  #data must contain entrytime, survtime and the covariates to RA on
  #glmmodel must be either of class glm or contain $formula and $coefficients
  #stoptime is the time until which the CUSUM chart should be constructed

  #------------------------------DATA CHECKS----------------------------------------#
  #Basic data checks (global for BK, CGR and Bernoulli)
  if(missing(data)){
    stop("Please provide data to construct chart.")
  } else{
    data <- check_data(data)
  }



  if(!missing(stoptime)){
    data <- subset(data, entrytime + followup <= stoptime)
  }
  #---------------------------FUNCTION BODY---------------------------#
  #Boolean indicating whether chart has been stopped by control limit h
  stopind = FALSE
  hnull <- missing(h)

  #Some checks:
  if(nrow(data) == 0){
    warning("No failures observed in specified time frame.
Decrease 'followup' or consider a larger time frame for construction.
Returning trivial chart.")
    Gt <- data.frame(time = c(0), value = c(0), numobs = c(0))
    colnames(Gt) = c("time", "value", "numobs")
    Ber <- list(CUSUM = Gt,
                call = call,
                stopind = stopind)
    if(!missing(glmmod)){
      Ber$glmmod <- glmmod$coefficients
    }
    if(!missing(h)){Ber$h <- h}
    class(Ber) <- "bercusum"
    Ber
    return(Ber)
  } else{
    min_entrytime <- min(data$entrytime)
  }
  #If twosided chart is required, determine the chart in two directions
  if(isTRUE(twosided)){
    Gt <- data.frame(time = c(min_entrytime), val_up = c(0), val_down = c(0), numobs = c(0))
    Gtval_up <- 0
    Gtval_down <- 0
    if(!hnull && length(h) == 1){
      h <- sort(c(-h, h))
    } else if(!hnull && length(h) == 2){
      if(!all(sign(sort(h)) == c(-1, 1))){
        stop("When specifying 2 control limits the two values should have reverse signs.")
      } else{
        h <- sort(h)
      }
    } else if(!hnull && length(h) > 2){
      stop("Please provide 1 or 2 values for the control limit.")
    }
  } else if(isFALSE(twosided)){
    Gt <- data.frame(time = c(min_entrytime), value = c(0), numobs = c(0))
    Gtval <- 0
    if(!hnull){
      if(length(h) > 1){
        stop("Please provide only 1 value for the control limit")
      }
      if(theta >= 0){
        h = abs(h)
      } else{
        h = -abs(h)
      }
    }
  }
  #Order the data by subject entry time
  data <- data[order(data$entrytime),]
  #Determine whether patient had failure. Censored observations do
  #not count as failures
  data$outcome <- as.integer((data$survtime <= followup) & (data$censorid == 1))
  data$otime <- data$entrytime + followup
  j <- 1
  numobs <- 0
  if(!missing(p1)){
    if(missing(p0) & missing(theta)){
      stop("Please also provide a value for p0 or theta.")
    }
    if(missing(p0) & !missing(theta)){
      p0 <- p1/(exp(theta) - exp(theta)*p1 + p1)
    }
    theta <- log((p1*(1-p0))/(p0*(1-p1)))

  }
  if(isTRUE(twosided)){
    theta = abs(theta)
  }


  #pre-calculate risk-adjustment
  if(!missing(glmmod)){
    if(inherits(glmmod, "glm")){
      fixprobs <-  predict(glmmod, newdata = data, type = "response")
    } else{
      mmatrix <- model.matrix(glmmod$formula, data)
      coeffs <- glmmod$coefficients[colnames(mmatrix)]
      fixprobs <- c(1/(1 + exp(-mmatrix %*% coeffs)))
    }
  }

  #Loop over all unique observation times and determine value of the chart
  for(i in unique(data$otime)){
    #Only interested in subjects with failure time at single time point
    tempdata_ind <- which(data$otime == i)
    tempdata <- data[tempdata_ind,]
    numobs <- numobs + nrow(tempdata)
    #If risk-adjustment has been specified, calculate risk scores.
    #Otherwise, use specified failure probabilities
    if(!missing(glmmod)){
      tempprobs <- fixprobs[tempdata_ind]
      tempsecondval <- sum(log(1/(1-tempprobs + exp(theta)*tempprobs)))
      if(isTRUE(twosided)){
        tempsecondval_down <- sum(log(1/(1-tempprobs + exp(-theta)*tempprobs)))
      }
    }else if(!missing(p0)){
      if(!missing(theta)){
        tempsecondval <- log((1/(1-p0 + exp(theta)*p0))^(nrow(tempdata)))
        if(isTRUE(twosided)){
          tempsecondval_down <- log((1/(1-p0 + exp(-theta)*p0))^(nrow(tempdata)))
        }
      } else if(!missing(p1)){
        tempsecondval <- log(((1-p1)/(1-p0))^(nrow(tempdata)))
      }
    } else{ stop("Please supply a value of theta or p1 or a glmmod")}


    #Determine chart values
    if(isTRUE(twosided)){
      #Determine W_n to update the value of the CUSUM with
      Wn_upper <- sum(tempdata$outcome)*theta +  tempsecondval
      Wn_lower <- -sum(tempdata$outcome)*theta +  tempsecondval_down
      Gtval_up <- max(0, Gtval_up + Wn_upper)
      Gtval_down <- min(0, Gtval_down - Wn_lower)
      Gt <- rbind(Gt, c(i, Gtval_up, Gtval_down, numobs))
    } else if(isFALSE(twosided)){
      #Determine W_n to update the value of the CUSUM with
      Wn <- sum(tempdata$outcome)*theta +  tempsecondval
      if(theta >= 0){
        Gtval <- max(0, Gtval + Wn)
      } else if(theta < 0){
        Gtval <- min(0, Gtval - Wn)
      }
      Gt <- rbind(Gt, c(i, Gtval, numobs))
    }


    #Determine whether to stop if h is specified
    if (!hnull){
      if(isTRUE(twosided)){
        if(length(h) == 2){
          if( (Gtval_up >= h[2]) | (Gtval_down <= h[1]) ) {stopind = TRUE; break}
        } else if(length(h) == 1){
          if( (abs(Gtval_up) >= abs(h)) | (abs(Gtval_down) >= abs(h)) ) {stopind = TRUE; break}
        }
      } else if(isFALSE(twosided)){
        if( abs(Gtval) >= abs(h) ) {stopind = TRUE; break}
      }
    }

    j <- j+1
  }
  if(isTRUE(twosided)){
    colnames(Gt) <- c("time", "val_up", "val_down", "numobs")
  } else if(isFALSE(twosided)){
    colnames(Gt) = c("time", "value", "numobs")
  }
  Ber <- list(CUSUM = Gt,
              call = call,
              stopind = stopind)
  if(!missing(glmmod)){
    Ber$glmmod <- glmmod$coefficients
  }
  if(!missing(h)){Ber$h <- h}
  class(Ber) <- "bercusum"
  Ber
}

