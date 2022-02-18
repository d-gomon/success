#' Continuous time Generalized Rapid response CUSUM (CGR-CUSUM) helper - matrix
#' formulation of the problem - version 2
#'
#' @description This function calculates the value of the CGR-CUSUM using a
#' matrix formulation of the problem - this can require a lot of available RAM.
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
#' @importFrom parallel detectCores
#' @importFrom parallel stopCluster
#' @importFrom parallel parSapply
#' @importFrom parallel makeCluster
#' @importFrom pbapply pbmapply
#' @importFrom pbapply pbsapply
#' @importFrom pbapply pbapply
#' @importFrom parallel clusterExport
#' @importFrom pbapply pboptions
#'
#' @noRd
#' @keywords internal
#'
#' @author Daniel Gomon
#'
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
#' cgr2 <- cgr_helper_mat_2(data = tdat, ctimes = unique(tdat$entrytime + tdat$survtime),
#'                          coxphmod = tcoxmod, cbaseh = tcbaseh, displaypb = TRUE)
#' }

cgr_helper_mat_down <- function(data, ctimes, h, coxphmod, cbaseh, ncores, displaypb = FALSE, dependencies){
  #!
  #DEVELOPMENT VERSION - DOWNWARDS CGR-CUSUM
  #!
  #This function consists of the following functions:
  #maxoverk: calculate CGI value at specific ctime, considering only patients with entrytime S_i >= k
  #maxoverj: Calculate CGR value at specific ctime by applying maxoverk on all relevant starting times (helperstimes)
  #fin = apply maxoverj to all required construction times to determine the values of CGR(t) for all ctimes
  #TRIPLE APPLY CONSTRUCTION.

  #Lambdamat is determined using an mapply on the columns of a matrix which is structured as:
  #      ctime_1  ctime_2 ctime_3  ctime_4 ................ ctime_p
  # i=1 (                                                            )
  # i=2 (     accumulated hazard of subject 2 at time in column      )
  # i=3 (                           subject 3                        )
  # ... (                             ...                            )
  # i=p (                                                            )
  #
  # Then we determine all patients which are relevant for constructing at a certain timepoint
  # by indexing over the data (checking their Stimes) and determine their active
  # contribution by summing over the column at the desired ctime.
  # The number of failures at that timepoint is determined the same way
  ctimes <- ctimes[which(ctimes >= min(data$entrytime))]

  #Remove R check warnings
  entrytime <- NULL

  if(ncores > 1){
    #Create a cluster for parallel computing and error check
    real_cores <- detectCores()
    if(ncores > real_cores){
      warning(paste0("More cores requested (", ncores ,") than detected (", real_cores,") by R. \n Proceed at own risk."))
    }
    cl <- makeCluster(ncores)
  } else{
    cl <- NULL
  }

  #Calculate Lambda_i for a single person at times - can be applied to matrix columns
  Lambdafun_noparallel <- function(entrytime, otime, times, cbaseh){
    lamval <- sapply(times, FUN = function(x) {
      if(x >= entrytime){
        cbaseh(min(x, otime) - entrytime)
      } else{0}
    })
  }
  #Same as lamdafun_noparallel, but can be applied to matrix rows using apply
  Lambdafun <- function(y, times, cbaseh){
    lamval <- sapply(times, FUN = function(x) {
      if(x >= y[1]){
        cbaseh(min(x, y[2]) - y[1])
      } else{0}
    })
  }

  #---------NEW FCT BODY-------------
  #Calculate subject specific failure risk (exp(\beta Z))
  riskdat <- calc_risk(data, coxphmod)
  #ctimes are already pre-sorted in cgrcusum
  #Create matrix containing Lambda_i(t) per patient. Rows represent patient i
  if(displaypb){
    message("Step 1/2: Determining hazard contributions.")
    pboptions(type = "timer")
  } else{
    pboptions(type = "none")
  }


  #MULTI-CORE:
  #Provide dependencies to cluster (such as functions required for clusters)
  if(!missing(dependencies)){
    clusterExport(cl,dependencies, envir=environment())
  }
  #Create matrix containing Lambda_i(t) per patient. Rows represent subject i, columns represent ctimes
  #the columns are exactly ctimes
  #tryCatch: if MULTI-CORE doesn't work, go back to SINGLE-CORE
  lambdamat <- tryCatch({t(pbapply(X = data[, c("entrytime", "otime")], MARGIN = 1,
                                   FUN = Lambdafun, times = ctimes, cbaseh = cbaseh, cl = cl))},
                        error = function(cond){
                          asd <- cond
                          message(cond)
                          message(paste0("\nPlease specify above missing functions/variables as a list to cgr_cusum under"
                                         ," the argument 'dependencies'"))
                          message("\nStep 1/2 will continue WITHOUT PARALLELIZATION - MAY BE SLOW!")
                          t(pbmapply(FUN = Lambdafun_noparallel, entrytime = data$entrytime,
                                     otime = data$otime, MoreArgs = list( times = ctimes, cbaseh = cbaseh)))
                        }
  )


  #Risk-adjust the cumulative intensity calculated at times
  lambdamat <- lambdamat * as.vector(riskdat)

  #Determine times from which to construct the CGR
  helperstimes <- sort(unique(data$entrytime))
  #helperfailtimes are not used in lower sided CGR-CUSUMS, because any
  #contribution to Lambda(t) will drift the chart downwards... :(
  #helperfailtimes <- numeric(length(helperstimes))
  #Determine failure times, as we only have to use helperstimes with failure times
  #smaller than the construction time.
  #for(i in seq_along(helperstimes)){
  #   failtemptime <- min(subset(data, entrytime == helperstimes[i])$otime)
  #   helperfailtimes[i] <- failtemptime
  #   if(failtemptime == helperstimes[i]){
  #     if(i > 1){
  #       helperfailtimes[i-1] <- failtemptime
  #     }
  #   }
  # }


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
      thetat <- min(0, thetat)
    } else {thetat <- 0}
    #Determine value of CGI-CUSUM using only patients with S_i > helperstimes[i]
    CGIvalue <- - thetat * NDT + (exp(thetat)- 1) * AT
    #Return both the value of CGI and the MLE (to be used later)
    return(c(CGIvalue, thetat))
  }

  #Function to calculate CGR value at one ctime by ways of maximizing over all CGI(t) values.
  #Achieved by applying the maxoverk function and determining maxima.
  maxoverj <- function(y){
    #don't use helperfailtimes in downward CUSUM
    #a <- sapply(helperstimes[which(helperstimes <= y & helperfailtimes <= y)], function(x) maxoverk(helperstime = x,  ctime = y, ctimes = ctimes, data = data, lambdamat = lambdamat))
    #We first apply the maxoverk function to determine CGI values
    #starting from all relevant helper S_i times (patient entry times)
    #We only need to calculate the CGI value when the the starting time of patients
    a <- sapply(helperstimes[which(helperstimes <= y)],
                function(x) maxoverk(helperstime = x,  ctime = y, ctimes = ctimes, data = data, lambdamat = lambdamat))
    #If there are no values to be calculated, return trivial values
    if(length(a) == 0){
      return(c(0,0,1))
    }else{
      #First row is value of chart, second row associated value of theta
      #Determine which entry is the smallest (largest CGI value)
      tidmin <- which.min(a[1,])
      #Determine the corresponding value of CGI(t)
      atemp <- a[,tidmin]
      #Returns c(chartval, thetaval) at maximum of CGI(t)
      return(c(atemp, tidmin))
    }
  }

  maxoverj_h <- function(y, h){
    #We first apply the maxoverk function to determine CGI values
    #starting from all relevant helper S_i times (patient entry times)
    #We only need to calculate the CGI value when the the starting time of patients
    #is before our construction time
    #For when I fix helperfailtimes problem.
    #a <- sapply(helperstimes[which(helperstimes <= y & helperfailtimes <= y)], function(x) maxoverk(helperstime = x,  ctime = y, ctimes = ctimes, data = data, lambdamat = lambdamat))
    a <- sapply(helperstimes[which(helperstimes <= y)],
                function(x) maxoverk(helperstime = x,  ctime = y, ctimes = ctimes, data = data, lambdamat = lambdamat))
    #First row is value of chart, second row associated value of theta
    #Determine which entry is the largest (largest CGI value)
    tidmin <- which.min(a[1,])
    #Determine the corresponding value of CGI(t)
    atemp <- a[,tidmin]
    if(abs(atemp[1]) >= abs(h)){
      hcheck <<- TRUE
      stopind <<- TRUE
      stopctime <<- match(y, ctimes)
    }
    #Returns c(chartval, thetaval) at maximum of CGI(t)
    return(c(atemp, tidmin))
  }

  #hcheck used to check whether value of CGR has surpassed control limit
  hcheck <- FALSE
  stopind <- FALSE
  stopctime <- NA

  #Calculate maxoverk for every construction time ctimes
  #If progress bar, then pbsapply (with possible multi-core)
  if(displaypb){
    message("Step 2/2: Determining chart values.")
  }
  if(!missing(h)){
    #Calculate maxoverk for every construction time ctimes by applying maxoverj on all ctimes.
    if(ncores > 1){
      stopCluster(cl)
      cl <- NULL
    }
    fin <- pbsapply(ctimes, function(x){ if(isFALSE(hcheck)){ maxoverj_h(x, h = h)} else{return(c(0,0,1))}}, cl = cl)
  } else{
    fin <- pbsapply(ctimes, maxoverj, cl = cl)
    if(ncores > 1){
      stopCluster(cl)
    }
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
  Gt <- rbind(c(min(data$entrytime), 0, 1, 0),Gt)
  colnames(Gt) <- c("time", "value", "exp_theta_t", "S_nu")

  #return list of relevant values
  return(Gt)
}
