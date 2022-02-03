#' @title Calculate the Cox Proportional hazards relative risk associated
#' with the covariates of subjects
#'
#' @description  This function can be used to calculate the change in the
#' relative risk of a subject pertaining to their covariates
#'  under a specified Cox proportional hazards model.
#'
#' @details The subject-specific relative risk is
#'  \eqn{e^{\beta  Z_i}}{exp(\beta * Z_i)},
#' where \eqn{\beta}{\beta} is a vector of regression coefficients
#' and \eqn{Z_i}{Z_i} a vector of covariates for subject i.
#'
#'
#' @param data A data frame containing the covariates to be used
#'  for risk-adjustment as named columns.
#' @param coxphmod (optional): A cox proportional hazards model generated using
#' \code{\link[survival:coxph]{coxph()}} or a list containing:
#' \describe{
#'   \item{\code{formula}:}{a \code{\link[stats:formula]{formula()}} in the form \code{~ covariates};}
#'   \item{\code{coefficients}:}{a named vector specifying risk adjustment coefficients
#'   for covariates. Names must be the same as in \code{formula} and colnames of \code{data}.}
#' }
#'
#'
#'
#' @return A vector of \code{nrow(data)} specifying the increased/decreased
#' risk of failure for each subject.
#'
#' @importFrom stats model.matrix
#' @importFrom stats terms
#' @importFrom stats update
#' @export
#'
#'
#' @author Daniel Gomon
#' @examples
#' #Small example data
#' crdat <- data.frame(age = rnorm(10, 40, 5), BMI = rnorm(10, 24, 3))
#' #Example risk-adjustment list (can also specify coxphmod)
#' crlist <- list(formula = as.formula("~age + BMI"), coefficients = c("age"= 0.02, "BMI"= 0.009))
#' #Calculate the increase or decrease of the relative risk for the subjects
#' #in crdat.
#' calc_risk(crdat, crlist)


calc_risk <- function(data, coxphmod = NULL){
  #Calculate risk for dataset data according to Cox PH model using specified model
  #coxphmod must either be a COXPH model or a list containing $formula and named vector $coefficients




  #Data checks
  if(missing(data)){
    return(1)
  } else if(is.null(data)){
    return(1)
  } else if(nrow(data) == 0){
    return(1)
  }
  #Coxphmod checks - if valid, determine risk
  if(missing(coxphmod) | is.null(coxphmod)){
    #When coxphmod is NULL, just return unadjusted probabilities
    return(rep(1, nrow(data)))
  } else if(length(labels(terms(coxphmod$formula))) == 0){
    return(rep(1, nrow(data)))
  } else{
    #Extract right side of model matrix from formula
    mmatrix <- model.matrix(update(coxphmod$formula, NULL ~ .), data)[,-1] #removes the intercept
    coeffs <- coxphmod$coefficients[colnames(mmatrix)]
    if(nrow(data) == 1){
      coeffs <- coxphmod$coefficients[names(mmatrix)]
    }
    #Return exp(beta * Z_i)
    return(c(exp(mmatrix %*% coeffs)))
  }
}
