#' @noRd
#' @keywords internal
#'
#' @author Daniel Gomon
#'
#' @param data data frame to check for compliance
#'

check_data <- function(data){
  # check data.frame and column names
  if(!is.data.frame(data)){
    warning("Provided data is not a data frame, attempting to convert.",
            immediate. = TRUE)
    data <- as.data.frame(data)
  }
  if(!"entrytime" %in% colnames(data)){
    stop("Entry time missing for subjects. Please specify them as named column
        'entrytime' in your data frame.")
  }
  if(!"survtime" %in% colnames(data)){
    stop("Survival time missing for subjects. Please specify them as named
          column 'survtime' in your data frame.")
  }
  if(!"censorid" %in% colnames(data)){
    warning("No censoring mechanism specified. Assuming data is uncensored.")
    data$censorid <- rep(1, nrow(data))
  }
  if(any(is.na(data))){
    stop("Please make sure 'data' has no missing values.")
  }
  #  compriskcheck <- "cause" %in% colnames(data)
  #  if(compriskcheck){
  #    message("Competing risks specified.")
  #  }
  return(data)
}
