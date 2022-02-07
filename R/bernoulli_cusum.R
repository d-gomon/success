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
#' @param theta The \eqn{theta}{theta} value used to specify the odds ratio
#'  \eqn{e^\theta}{e^\theta} under the alternative hypothesis.
#' Note that \deqn{p_1 = \frac{p_0 e^\theta}{(1-p_0) (1+p_0 e^\theta)}.}{p1 = (p0 * e^\theta)/((1-p0) * (1+p0 e^\theta)).}
#' @param p0 The baseline failure probability at \code{entrytime + followup} for individuals.
#' @param p1 The alternative hypothesis failure probability at \code{entrytime + followup} for individuals.
#' @param h (optional): Control limit to be used for the procedure.
#' @param stoptime (optional): Time after which the value of the chart should no longer be determined.
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
#' @export
#'
#' @author Daniel Gomon
#' @family quality control charts
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




bernoulli_cusum <- function(data, followup, glmmod, theta, p0, p1, h, stoptime){
  entrytime <- otime <- NULL
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
  #Order the data by subject entry time
  data <- data[order(data$entrytime),]
  #Determine whether patient had failure. Censored observations do
  #not count as failures
  data$outcome <- as.integer((data$survtime <= followup) & (data$censorid == 1))
  data$otime <- data$entrytime + followup
  Gt <- data.frame(time = double(), value = double(), numobs = double())
  Gtval <- 0
  j <- 1
  numobs <- 0
  if(!missing(p1)){
    OR = (p1*(1-p0))/(p0*(1-p1))
  }
  #Loop over all unique observation times and determine value of the chart
  for(i in unique(data$otime)){
    #Only interested in subjects with failure time at single time point
    tempdata <- subset(data, otime == i)
    numobs <- numobs + nrow(tempdata)
    #If risk-adjustment has been specified, calculate risk scores.
    #Otherwise, use specified failure probabilities
    if(!missing(glmmod)){
      if(inherits(glmmod, "glm")){
        tempprobs <-  predict(glmmod, newdata = tempdata, type = "response")
      } else{
        mmatrix <- model.matrix(glmmod$formula, tempdata)
        coeffs <- glmmod$coefficients[colnames(mmatrix)]
        tempprobs <- c(1/(1 + exp(-mmatrix %*% coeffs)))
      }
      tempsecondval <- sum(log(1/(1-tempprobs + exp(theta)*tempprobs)))
    }else if(!missing(p0)){
      if(!missing(theta)){
        tempsecondval <- log((1/(1-p0 + exp(theta)*p0))^(nrow(tempdata)))
      } else if(!missing(p1)){
        tempsecondval <- log(((1-p1)/(1-p0))^(nrow(tempdata)))
      }
    } else{ stop("Please supply a value of theta or p1 or a glmmod")}
    #Determine W_n to update the value of the CUSUM with
    if(!missing(theta)){
      Wn <- sum(tempdata$outcome)*theta +  tempsecondval
    } else{
      Wn <- sum(tempdata$outcome)*log(OR) +  tempsecondval
    }
    Gtval <- max(0, Gtval + Wn)
    #Push the new value to the data frame
    Gt <- rbind(Gt, c(i, Gtval, numobs))
    if(!hnull){if(Gtval >= h){Ber <- list(CUSUM = Gt,
                                          call = call,
                                          h = h,
                                          stopind = TRUE)}
      if(!missing(glmmod)){
        Ber$glmmod <- glmmod$coefficients
      }
      class(Ber) <- "bercusum"
      return(Ber)
      }
    j <- j+1
  }
  colnames(Gt) = c("time", "value", "numobs")
  Ber <- list(CUSUM = Gt,
              call = call,
              stopind = stopind)
  if(!missing(glmmod)){
    Ber$glmmod <- glmmod$coefficients
  }
  class(Ber) <- "bercusum"
  Ber
}

