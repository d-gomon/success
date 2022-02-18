


gbk_cusum <- function(data, coxphmod, cbaseh, ctimes, h, stoptime,
                      C, pb = FALSE, ncores = 1, cmethod = "memory",
                      dependencies, detection = "upper", assist){

  if(!missing(assist)){
    list2env(assist, envir = environment())
  }
  call = match.call()
  #-------------------------------DATA CHECKS-----------------------------#
  #Basic data checks (global for BK, CGR and Bernoulli)
  if(missing(data)){
    stop("Please provide data to construct chart.")
  } else{
    data <- check_data(data)
  }
  #  compriskcheck <- "cause" %in% colnames(data)
  #  if(compriskcheck){
  #    message("Competing risks specified.")
  #  }
  # determine chronological failure times
  data$otime <- data$entrytime + data$survtime
  if(!missing(C)){
    tempidx <- which(data$otime < data$entrytime + C)
    data[tempidx,]$otime <- data$entrytime + C
    data[tempidx,]$censorid <- rep(length(tempidx), 0)
  }
  # determine the default construction times (all failtimes), if none specified
  if(missing(ctimes)){
    ctimes <- unique(data$otime)
  }
  if(missing(stoptime)){
    stoptime <- max(data$otime[is.finite(data$otime)])
  }
  checkcbase <- FALSE
  if(missing(coxphmod)){
    coxphmod <- NULL
  } else if(inherits(coxphmod, "coxph") & missing(cbaseh)){
    checkcbase <- TRUE
    message("Missing cumulative baseline hazard. Determining using provided Cox PH model.")
    cbaseh <- extract_hazard(coxphmod)$cbaseh
    #Old Loess method:
    # cbaselo <- loess(cbase_temp$hazard~cbase_temp$time)
    # cbaseh <- function(x) predict(cbaselo, x)
  } else if(is.list(coxphmod)){
    if(all(c("formula", "coefficients") %in% names(coxphmod))){
      checkcoxlist <- TRUE
    } else if(!is.null(coxphmod)){
      stop("coxphmod does not contain $formula and/or $coefficients.")
    }
  } else if(is.null(coxphmod)){
    coxphmod <- NULL
  } else{ stop("coxphmod is not a list or survival object.")}
  if(missing(cbaseh)){
    if(!checkcbase){
      stop("Please specify cbaseh (function) or coxphmod as Survival object.")
    }
  }
  if(!missing(h)){
    if(length(h) > 1){
      stop("Please specify only 1 control limit.")
    }
  }

  #----------------------------FUNCTION BODY--------------------------#

  #Only determine value of the chart at at relevant times
  ctimes <- sort(ctimes[which(ctimes <= stoptime)])


  #If method = CPU, iteratively calculate the value at each time
  Gt <- matrix(0, nrow = 2*length(ctimes)+1, ncol = 4)
  startval <- 0
  stopind <- FALSE

  minval <- 0
  #progress bar
  if(pb){
    pb2 <- txtProgressBar(min= 0, max = length(ctimes), style = 3)
  }
  #Calculate value of chart at each construction time separately
  for(j in seq_along(ctimes)){
    if(pb){
      setTxtProgressBar(pb2, value = j)
    }
    #Calculate value of CGR at time ctimes[j]
    temgbk <- gbk_helper(data = data, ctime = ctimes[j],
                         coxphmod = coxphmod, cbaseh = cbaseh)
    if(temgbk$val_down <= minval){
      minval <- temgbk$val_down
    }
    Gt[j+1,] <- c(ctimes[j], temgbk$val_down - minval, temgbk$val_up - minval, exp(temgbk$theta))
    if (!missing(h)){if(temgbk$val_up >= h) {stopind = TRUE; break}}
  }
  colnames(Gt) <- c("time", "val_down", "val_up", "exp_theta_t")
  if(pb){
    close(pb2)
  }
  gbk <- list(GBK = Gt,
              stopind = stopind,
              call = call)


  if(!missing(h)){
    if(detection == "upper"){
      cgr$h <- abs(h)
    } else{
      cgr$h <- -abs(h)
    }
  }
  class(gbk) <- "gbkcusum"
  gbk
}
