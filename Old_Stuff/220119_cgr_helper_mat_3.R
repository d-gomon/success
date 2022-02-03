#' Continuous time Generalized Rapid response CUSUM (CGR-CUSUM) helper - matrix
#' formulation of the problem - version 3
#'
#' @description This function calculates the value of the CGR-CUSUM using a
#' matrix formulation of the problem - reduce calculations by specifying control limit.
#'
#'
#' @inheritParams cgr_cusum
#' @param data \code{data.frame} containing the following named columns:
#' \itemize{
#' \item \code{entrytime} numeric - time of entry into study,
#' \item \code{otime} numeric - time from entry until event,
#' \item \code{censorid} integer - (optional) censoring indicator (0 = right censored, 1 = observed),
#\item \code{cause} factor - cause of event - competing risks.
#' } and optionally additional covariates used for risk-adjustment.
#' @param displaypb boolean Display a progress bar?
#'
#' @return A matrix with 4 named columns:
#' \itemize{
#'   \item $time time at which the value of the CGR-CUSUM was determined
#'   \item $value value at corresponding time of the CGR-CUSUM
#'   \item $exp_theta_t value at corresponding time of the MLE \eqn{\hat{\theta}_t}{\theta_t}
#'   \item $S_nu time from which individuals contribute to the chart \eqn{S_\nu}{S_\nu}
#' }
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @noRd
#' @keywords internal
#'
#' @author Daniel Gomon
#'
#' @seealso \code{\link{cgrcusum}}
#'
#' @examples
#' \dontrun{
#' require(survival)
#' tdat <- subset(surgerydat, hosp_num == 1)
#' tdat$otime <- tdat$entrytime + tdat$survtime
#' tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
#' varsanalysis <- c("age", "sex", "BMI")
#' exprfit <- as.formula(paste("Surv(survtime, censorid) ~" ,paste(varsanalysis, collapse='+')))
#' tcoxmod <- coxph(exprfit, data= surgerydat)
#' #Alternatively, cbaseh can be left empty when specifying coxphmod through coxph()
#' cgr3 <- cgr_helper_mat_3(data = tdat, ctimes = unique(tdat$entrytime + tdat$survtime),
#'                          coxphmod = tcoxmod, cbaseh = tcbaseh, displaypb = TRUE)
#' }

cgr_helper_mat_3 <- function(data, ctimes, coxphmod, cbaseh, h, displaypb = FALSE){
  #WITH HCHECK!
  #This function consists of the following functions:
  #maxoverk: calculate CGI value at specific ctime, considering only patients with entrytime S_i >= k
  #maxoverj: Calculate CGR value at specific ctime by applying maxoverk on all relevant starting times (helperstimes)
  #fin = apply maxoverj to all required construction times to determine the values of CGR(t) for all ctimes
  #TRIPLE APPLY CONSTRUCTION.

  #Lambdamat is determined using an mapply on the columns of a matrix which is structured as:
  #      ctime_1  ctime_2 ctime_3  ctime_4 ................ ctime_k
  # i=1 (                                                            )
  # i=2 (     accumulated hazard of subject 2 at time in column      )
  # i=3 (                           subject 3                        )
  # ... (                             ...                            )
  # i=n (                                                            )
  #
  # Then we determine all patients which are relevant for constructing at a certain timepoint
  # by indexing over the data (checking their Stimes) and determine their active
  # contribution by summing over the column at the desired ctime.
  # The number of failures at that timepoint is determined the same way



  #Calculate Lambda_i for a single person at times - can be applied to matrix columns
  Lambdafun <- function(entrytime, otime, times, cbaseh){
    lamval <- sapply(times, FUN = function(x) {if(x >= entrytime){cbaseh(min(x, otime) - entrytime)}else{0}})
  }

  #-----------------------------------FUNCTION BODY------------------------------------#
  #Calculate subject specific failure risk (exp(\beta Z))
  riskdat <- calc_risk(data, coxphmod)
  #ctimes are already pre-sorted in cgrcusum
  #Create matrix containing Lambda_i(t) per patient. Rows represent subject i, columns represent ctimes
  #the columns are exactly ctimes
  lambdamat <- t(mapply(FUN = Lambdafun, entrytime = data$entrytime, otime = data$otime,
                        MoreArgs = list( times = ctimes, cbaseh = cbaseh)))
  #Risk-adjust the cumulative intensity calculated at times
  lambdamat <- lambdamat * as.vector(riskdat)
  #Determine times from which to construct the CGR
  helperstimes <- sort(unique(data$entrytime))

  #Function used for maximizing over starting points (patients with starting time >= k)
  maxoverk <- function(helperstime, ctime, ctimes, data, lambdamat){
    #Determine part of data that is active at required times
    matsub <- which(data$entrytime >= helperstime & data$entrytime <= ctime)
    #The cumulative intensity at that time is the column sum of the specified ctime
    AT <- sum(lambdamat[matsub, which(ctimes == ctime)])
    #THIS COULD BE SLOW, OTHERWISE ASSIGN TDAT <- subset(data, matsub)
    #Determine amount of failures at ctime.
    NDT <- length(which(data[matsub, ]$censorid == 1 & data[matsub,]$otime <= ctime))
    #Determine MLE of theta
    thetat <- log(NDT/AT)
    if (is.finite(thetat)){
      thetat <- max(0,thetat)
    } else {thetat <- 0}
    #Determine value of CGI-CUSUM using only patients with S_i > helperstimes[i]
    CGIvalue <- thetat* NDT - (exp(thetat)- 1) * AT
    #Return both the value of CGI and the MLE (to be used later)
    return(c(CGIvalue, thetat))
  }
  #Progress bar
  if(displaypb){
    pb <- txtProgressBar(min= 0, max = length(ctimes), style = 3)
    #Function to calculate CGR value at one ctime by ways of maximizing over all CGI(t) values.
    #Achieved by applying the maxoverk function and determining maxima.
    maxoverj <- function(y, h){
      setTxtProgressBar(pb, value = match(y, ctimes))
      #We first apply the maxoverk function to determine CGI values
      #starting from all relevant helper S_i times (patient entry times)
      #We only need to calculate the CGI value when the the starting time of patients
      #is before our construction time
      a <- sapply(helperstimes[which(helperstimes <= y)], function(x) maxoverk(helperstime = x,  ctime = y, ctimes = ctimes, data = data, lambdamat = lambdamat))
      #First row is value of chart, second row associated value of theta
      #Determine which entry is the largest (largest CGI value)
      tidmax <- which.max(a[1,])
      #Determine the corresponding value of CGI(t)
      atemp <- a[,tidmax]
      if(atemp[1] >= h){
        hcheck <<- TRUE
        stopind <<- TRUE
        stopctime <<- match(y, ctimes)
      }
      #Returns c(chartval, thetaval) at maximum of CGI(t)
      return(c(atemp, tidmax))
    }
  } else{
    #Same as above, but without progress bar
    maxoverj <- function(y, h){
      a <- sapply(helperstimes[which(helperstimes <= y)], function(x) maxoverk(helperstime = x,  ctime = y,
                                                     ctimes = ctimes, data = data, lambdamat = lambdamat))
      tidmax <- which.max(a[1,])
      atemp <- a[,tidmax]
      #Check whether chart has surpassed control limit h
      if(atemp[1] >= h){
        hcheck <<- TRUE
        stopind <<- TRUE
        #Time at which the chart was stopped
        stopctime <<- match(y, ctimes)
      }
      return(c(atemp, tidmax))
    }
  }
  #hcheck used to check whether value of CGR has surpassed control limit
  hcheck <- FALSE
  stopind <- FALSE
  stopctime <- NA
  #Calculate maxoverk for every construction time ctimes by applying maxoverj on all ctimes.
  fin <- sapply(ctimes, function(x){ if(isFALSE(hcheck)){ maxoverj(x, h = h)} else{return(c(0,0,1))}})
  if(displaypb){
    close(pb)
  }
  Gt <- as.data.frame(t(fin))
  if(!is.na(stopctime)){
    Gt <- Gt[1:stopctime,]
  }
  Gt[,2] <- exp(Gt[,2])
  Gt[,3] <- helperstimes[Gt[,3]]
  if(!is.na(stopctime)){
    Gt <- cbind(ctimes[1:stopctime], Gt)
  } else{
    Gt <- cbind(ctimes, Gt)
  }
  colnames(Gt) <- c("time", "value", "exp_theta_t", "S_nu")

  #return list of relevant values
  return(list(Gt = Gt,
              stopind = stopind))
}
