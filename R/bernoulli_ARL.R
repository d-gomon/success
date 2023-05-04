#' Average Run Length for Bernoulli CUSUM
#'
#' @description This function allows to estimate the Average Run Length (ARL)
#' of the risk-adjusted Bernoulli CUSUM (see \code{\link[success:bernoulli_cusum]{bernoulli_cusum()}})
#' through a Markov Chain Approach (Brook & Evans(1972) & Steiner et al. (2000)).
#' The function requires the specification of one of the following combinations of parameters
#' as arguments to the function:
#' \itemize{
#' \item \code{glmmod} & \code{theta}
#' \item \code{p0} & \code{theta}
#' \item \code{p0} & \code{p1}
#' } Average run length of lower-sided Bernoulli CUSUM charts can be determined
#' by specifying \code{theta} < 0.
#'
#' @param h Control limit for the Bernoulli CUSUM
#' @param t Number of state spaced used to discretize the outcome space.
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
#'  If \eqn{\theta >= 0}{\theta >= 0}, the average run length for the upper one-sided
#'  Bernoulli CUSUM will be determined. If \eqn{\theta < 0}{\theta < 0},
#' the average run length for the lower one-sided CUSUM will be determined.
#' Note that \deqn{p_1 = \frac{p_0 e^\theta}{1-p_0 +p_0 e^\theta}.}{p1 = (p0 * e^\theta)/(1-p0+p0 * e^\theta).}
#' @param p0 The baseline failure probability at \code{entrytime + followup} for individuals.
#' @param p1 The alternative hypothesis failure probability at \code{entrytime + followup} for individuals.
#' @param followup (optional) The value of the follow-up time used to determine event time. Will be used to
#' calculate Average Run Length in time units instead of number of observed outcomes.
#' @param smooth_prob Should the probability distribution of failure under the null distribution be smoothed?
#' Useful for small samples. Can only be TRUE when \code{glmmod} is supplied. Default = FALSE.
#'
#'
#'
#' @return A list containing:
#' \itemize{
#' \item \code{ARL_0}: A numeric value indicating the average run length in number of outcomes and
#' time (if followup specified) when starting from state E_0.
#' \item \code{ARL}: A \code{data.frame} containing the average run length (#outcomes/time)
#' depending on the state in which the process starts (E_0, E_1, ..., E_{t-1})
#' \describe{
#'   \item{\code{t_start}:}{State in which the CUSUM process starts;}
#'   \item{\code{#outcomes}:}{ARL from starting state in #outcomes;}
#'   \item{\code{Time}:}{(only if followup specified) ARL in time units;}
#' }
#' \item \code{R}: A transition probability \code{matrix} containing the transition
#' probabilities between states \eqn{E_0, \ldots, E_{t-1}}{E_0, ..., E_{t-1}}.
#' \eqn{R_{i,j}}{R_{i,j}} is the transition probability from state i to state j.
#' }
#'
#' @references Brook, D., & Evans, D. A. (1972). An Approach to the Probability
#' Distribution of Cusum Run Length. Biometrika, 59(3), 539–549.
#' \doi{10.2307/2334805}
#'
#' Steiner, S. H., Cook, R. J., Farewell, V. T., & Treasure, T. (2000).
#' Monitoring surgical performance using risk-adjusted cumulative sum charts.
#' Biostatistics, 1(4), 441–452. \doi{10.1093/biostatistics/1.4.441}
#'
#'
#' @details The average run length of a CUSUM chart \eqn{S_n}{S_n} is given by
#' \eqn{\mathbb{E}[\tau_n],}{E[\tau_n],} where \eqn{\tau_n}{\tau_n} is defined as:
#' \deqn{\tau_n = \inf\{n \geq 0: S_n \geq h\}.}{\tau_n = inf(n >= 0: S_n >= h).}
#' This function determines the average run length using the Markov Chain approach described
#' in Brook & Evans (1972), using the risk-adjustment correction proposed in
#' Steiner et al. (2000). The idea is to discretize the domain (0, h) into $t-1$
#' state spaces \eqn{E_0}{E_0} of width \eqn{w/2}{w/2} and \eqn{E_1, \ldots, E_{t-1}}{E_1, ..., E_{t-1}} of width \eqn{w}{w}, such that
#' \eqn{E_t}{E_t} is an absorbing state.  This is done using the following steps:
#' \itemize{
#' \item \eqn{w}{w} is determined using the relationship \eqn{\frac{2h}{2t-1}}{2h/(2t-1)}.
#' \item Transition probabilities between the states are determined and
#' 'transition matrix' \eqn{R}{R} is constructed.
#' \item The equation \eqn{(\bm{I}-\bm{R}) \bm{ARL} = \bm{1}}{(I-R) ARL = 1} is
#' solved to find the ARL starting from each of the states.
#' }
#'
#'
#'
#' @export
#' @importFrom Matrix solve
#'
#' @author Daniel Gomon
#' @family average run length
#' @seealso \code{\link[success]{bernoulli_cusum}}, \code{\link[success]{bernoulli_control_limit}}
#'
#'
#' @examples
#' #We consider patient outcomes 100 days after their entry into the study.
#' followup <- 100
#' #Determine a risk-adjustment model using a generalized linear model.
#' #Outcome (failure within 100 days) is regressed on the available covariates:
#' exprfitber <- as.formula("(survtime <= followup) & (censorid == 1)~ age + sex + BMI")
#' glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
#' ARL <- bernoulli_ARL(h = 2.5, t = 50, glmmod = glmmodber, theta = log(2))
#'
#'




bernoulli_ARL <- function(h, t, glmmod, theta, theta_true, p0, p1, followup, smooth_prob = FALSE){
  #------------Variable checks------------------

  #Input checks
  if(!all(is.numeric(h) & length(h) == 1)){
    stop("Parameter 'h' must be a single numeric variable.")
  }
  if(!all(is.numeric(t), t%%1 == 0, t > 1, length(t) == 1)){
    stop("Parameter 't' must be an integer larger than 1.")
  }
  if(!all(is.logical(smooth_prob), length(smooth_prob) == 1)){
    stop("Parameter 'smooth_prob' must be a single logical value.")
  }
  if(!missing(followup)){
    if(!all(is.numeric(followup), followup > 0, length(followup) == 1)){
      stop("Parameter 'followup' must be a positive numeric variable.")
    }
  }
  if(!missing(theta_true)){
    if(!all(is.numeric(theta_true) & theta != 0)){
      stop("Parameter 'theta' must be a numeric variable not equal to 0.")
    }
  } else{
    theta_true = NULL
  }
  if(!missing(glmmod)){
    if(!inherits(glmmod, "glm")){
      stop("Parameter 'glmmod' must be of class 'glm'.")
    }
  }
  if(!missing(theta)){
    if(!all(is.numeric(theta) & theta != 0)){
      stop("Parameter 'theta' must be a numeric variable not equal to 0.")
    }
  }
  if(!missing(p0)){
    if(!all(is.numeric(p0), p0 >= 0, p0 <= 1)){
      stop("Parameter 'p0' must be a positive probability between 0 and 1. (numeric)")
    }
  }
  if(!missing(p1)){
    if(!all(is.numeric(p1), p1 >= 0, p1 <= 1)){
      stop("Parameter 'p1' must be a positive probability between 0 and 1. (numeric)")
    }
  }

  if(h < 0){
    h <- sign(h)*h
  }

  #Correct combination of parameters specified?
  if(!((!missing(p0) & !missing(p1)) | (!missing(p0) & !missing(theta)) | (!missing(glmmod) & !missing(theta)))){
    stop("Please specify any of the following parameter combinations: 'glmmod' & 'theta' or 'p0' & 'theta' or 'p0' & 'p1'.")
  } else if(!missing(p0) & !missing(p1)){
    theta <- log((p1)*(1-p0)/((p0)*(1-p1)))
  } else if(!missing(p0) & !missing(theta)){
    p1 <- p0*exp(theta)/(1-p0 + exp(theta)*p0)
  } else if(!missing(glmmod) & !missing(theta)){
    #Do nothing (for now at least)
  }
  if(isTRUE(smooth_prob) & missing(glmmod)){
    stop("Probability distribution can only be smoothed when 'glmmod' is specified.")
  }

  #Calculate the value w used for discretizing the state space
  w <- 2*h/(2*t -1)
  #Initiatlize transition matrix
  R <- matrix(0, nrow = t, ncol = t)

  #If we do not want to smooth probability distribution
  if(isFALSE(smooth_prob)){
    if(!missing(glmmod)){
      null_probs <- sort(glmmod$fitted.values)
      if(theta > 0){
        #If we have glmmod, we can do a risk-adjusted calculation of the ARL.
        Wncdf_discrete <- function(x, null_probs, theta, theta_true = NULL){
          #Function to calculate cdf of W_n in discrete case
          if(!is.null(theta_true)){
            adjnull_probs <- null_probs*exp(theta_true)/(1 - null_probs + exp(theta_true)*null_probs)
          } else{
            adjnull_probs <- null_probs
          }
          N <- length(null_probs)
          if(x <= 0){
            #ind <- which(null_probs >= (1-exp(-x))/(1-exp(theta)))
            #Finding the probabilities which are larger than cut-off value (see above line)
            ind <- Rfast::binary_search(null_probs, (1-exp(-x))/(1-exp(theta)), index = TRUE)
            if(ind == N + 1){
              ind <- NULL
            } else{
              ind <- ind:N
            }
            #\sum_{p \geq (1-e^{-K})/(1-e^\theta)} (1-p) P(p_n = p)
            #or
            #\sum_{p \geq (1-e^{-K})/(1-e^\theta)} (1-(p*e^{theta_true}/(1-p + e^{theta_true}*p))) P(p_n = p)
            return(sum((1-adjnull_probs[ind]))*1/N)
          } else{
            #ind <- which(null_probs >= (1-exp(-(x-theta)))/(1-exp(theta)))
            ind <- Rfast::binary_search(null_probs, (1-exp(-(x-theta)))/(1-exp(theta)), index = TRUE)
            if(ind == N + 1){
              ind <- NULL
            } else{
              ind <- ind:N
            }
            #P(W_n \leq 0) + \sum_{p \geq (1-e^{-(K-\theta)})/(1-e^\theta)} p* P(p_n = p)
            return(Wncdf_discrete(0, null_probs = null_probs, theta = theta, theta_true = theta_true) + sum(adjnull_probs[ind])*1/N)
          }
        }
      } else if(theta < 0){
        #If we have glmmod, we can do a risk-adjusted calculation of the ARL.
        Wncdf_discrete <- function(x, null_probs, theta, theta_true = NULL){
          #Function to calculate cdf of W_n in discrete case
          N <- length(null_probs)
          if(!is.null(theta_true)){
            adjnull_probs <- null_probs*exp(theta_true)/(1 - null_probs + exp(theta_true)*null_probs)
          } else{
            adjnull_probs <- null_probs
          }
          if(x <= 0){
            #ind <- which(null_probs <= (1-exp(-x-theta))/(1-exp(theta)))
            #Finding the probabilities which are larger than cut-off value (see above line)
            ind <- Rfast::binary_search(null_probs, (1-exp(-(x-theta)))/(1-exp(theta)), index = TRUE)
            if(ind == 1){
              ind <- NULL
            } else{
              ind <- 1:(ind-1)
            }
            #\sum_{p \leq (1-e^{-K})/(1-e^\theta)} (1-p) P(p_n = p)
            #or
            #\sum_{p \leq (1-e^{-K})/(1-e^\theta)} (1-(p*e^{theta_true}/(1-p + e^{theta_true}*p))) P(p_n = p)
            return(sum(adjnull_probs[ind])*1/N)
          } else{
            #ind <- which(null_probs <= (1-exp(-(x-theta)))/(1-exp(theta)))
            ind <- Rfast::binary_search(null_probs, (1-exp(-(x)))/(1-exp(theta)), index = TRUE)
            if(ind == 1){
              ind <- NULL
            } else{
              ind <- 1:(ind-1)
            }
            #P(W_n \leq 0) + \sum_{p \geq (1-e^{-(K-\theta)})/(1-e^\theta)} p* P(p_n = p)
            return(Wncdf_discrete(0, null_probs = null_probs, theta = theta, theta_true = theta_true) + sum((1-adjnull_probs[ind]))*1/N)
          }
        }
      }
    }else{
      #if no glmmod, we don't perform risk-adjustment.
      #IMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANT
      #We use null_probs to represent p0
      #IMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANT
      if(theta > 0){
        Wncdf_discrete <- function(x, null_probs, theta, theta_true = NULL){
          #Function to calculate cdf of W_n in discrete unadjusted
          if(!is.null(theta_true)){
            adjnull_probs <- null_probs*exp(theta_true)/(1 - null_probs + exp(theta_true)*null_probs)
          } else{
            adjnull_probs <- null_probs
          }
          if(x <= 0){
            #Read line below as if(p0 >= ....)
            if(null_probs >= (1-exp(-x))/(1-exp(theta))){
              return(1-adjnull_probs)
            } else{
              return(0)
            }
          } else{
            if(null_probs >= (1-exp(-(x-theta)))/(1-exp(theta))){
              return(1)
            }else{
              return(1-adjnull_probs)
            }
          }
        }
      } else if(theta < 0){
        Wncdf_discrete <- function(x, null_probs, theta, theta_true = NULL){
          if(!is.null(theta_true)){
            adjnull_probs <- null_probs*exp(theta_true)/(1 - null_probs + exp(theta_true)*null_probs)
          } else{
            adjnull_probs <- null_probs
          }
          #Function to calculate cdf of W_n in discrete unadjusted
          if(x <= 0){
            #Read line below as if(p0 <= ....)
            if(null_probs <= (1-exp(-(x-theta)))/(1-exp(theta))){
              return(adjnull_probs)
            } else{
              return(0)
            }
          } else{
            if(null_probs <= (1-exp(-(x)))/(1-exp(theta))){
              return(1)
            }else{
              return(adjnull_probs)
            }
          }
        }
      }
    }
    #Using Wncdf_discrete we can calculate the transition probabilities
    trans_prob <- function(state1, state2, w, t, null_probs, theta, theta_true = NULL){
      if(state2 == 0){
        #P(W_n \leq -i2 + 1/2*w)
        Wncdf_discrete((0.5-state1) * w, null_probs = null_probs, theta = theta, theta_true = theta_true)
      } else if(state2 == t){
        1- Wncdf_discrete((t-state1 - 0.5) * w, null_probs = null_probs, theta = theta, theta_true = theta_true)
      }else{
        Wncdf_discrete((state2-state1 + 0.5) * w, null_probs = null_probs, theta = theta, theta_true = theta_true) - Wncdf_discrete((state2-state1 - 0.5) * w, null_probs = null_probs, theta = theta, theta_true = theta_true)
      }
    }
    #We only need to calculate some specific transition probabilities:
    if(missing(glmmod)){
      #To avoid code duplication, define null_probs as p0 here.
      null_probs = p0
    }


    #First we calculate the probabilities in the first row and first column of R
    #For this we create the vector (p_{t-1, 0}, p_{t-2, 0}, ..., p_{0,0}, ..., p_{0, t-2}, p_{0, t-1})
    probstmin10 <- sapply((t-1):0, function(x) trans_prob(x, 0, w = w, t = t, null_probs = null_probs, theta = theta, theta_true = theta_true))
    probs1tmin1 <- sapply(1:(t-1), function(x) trans_prob(0, x, w = w, t = t, null_probs = null_probs, theta = theta, theta_true = theta_true))
    probs0 <- c(probstmin10, probs1tmin1)

    #Then we calculate the probabilities on the diagonal of the rest of the matrix:
    #These are given by the vector (p_{t-1, 1}, ..., p_{1, 1}, p_{1, 2}, ..., p_{1, t-1})
    probstmin11 <- sapply((t-1):1, function(x) trans_prob(x, 1, w = w, t = t, null_probs = null_probs, theta = theta, theta_true = theta_true))
    probs2tmin1 <- sapply(2:(t-1), function(x) trans_prob(1, x, w = w, t = t, null_probs = null_probs, theta = theta, theta_true = theta_true))
    probs1 <- c(probstmin11, probs2tmin1)

    #Iterate through matrix to fill values. For details, see Overleaf file (picture displaying the matrix).
    for(i in 1:t){
      for(j in 1:t){
        #If we are in the first row or column, we use probs0 vector
        if(i == 1| j == 1){
          R[i,j] <- probs0[t + (j-i)]
        } else{
          #In any other row/column we use the probs1 vector
          R[i,j] <- probs1[(t-1) + (j-i)]
        }
      }
    }
  } else if(isTRUE(smooth_prob)){
    #smooth_prob is only usefull if glmmod is supplied.
    if(missing(glmmod)){
      stop("Probabilities can only be smoothed when 'glmmod' is supplied.")
    }
    density_estimate <- density(glmmod$fitted.values)

  }

  #Now we need to decompose the "transition matrix" R to obtain the ARL vector mu
  #Each entry mu[i] represents the ARL when starting from state i-1
  #First entry is therefore the average run length when starting from 0.
  mu = round(Matrix::solve(diag(1, nrow = t, ncol = t) - R, matrix(1, nrow = t, ncol = 1)))
  rownames(mu) <- 0:(t-1)

  #Post-processing
  if(!missing(followup)){
    ARL <- cbind(mu, mu*followup)
    colnames(ARL) <- c("#outcomes", "Time")
    ARL_0 = c(mu[1,])
  } else{
    ARL <- mu
    colnames(ARL) <- c("#outcomes")
    ARL_0 = mu[1]
    names(ARL_0) <- "#outcomes"
  }
  ARL <- cbind(0:(t-1), ARL)
  colnames(ARL)[1] <- "t_start"

  print(ARL_0)
  return(invisible(list(ARL_0 = ARL_0,
              ARL = ARL,
              R = R)))
}






