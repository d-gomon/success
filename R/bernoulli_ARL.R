#' Average Run Length for Bernoulli CUSUM
#'
#' @description This function allows to estimate the Average Run Length (ARL)
#' of the risk-adjusted Bernoulli CUSUM (see \code{\link[success:bernoulli_cusum]{bernoulli_cusum()}})
#' through a Markov Chain Approach (Brook & Evans(1972) & Steiner et al. (2000)) or
#' an Integral Equation Approach (Kemp (1971)).
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
#' @param n_grid Number of state spaces used to discretize the outcome space (when \code{method = "MC"})
#' or number of grid points used for trapezoidal integration (when \code{method = "IntEq"}).
#' Increasing this number improves accuracy, but can also significantly increase computation time.
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
#' @param theta_true The true log odds ratio \eqn{\theta}{\theta}, describing the
#' true increase in failure rate from the null-hypothesis. Default = log(1), indicating
#' no increase in failure rate.
#' @param p0 The baseline failure probability at \code{entrytime + followup} for individuals.
#' @param p1 The alternative hypothesis failure probability at \code{entrytime + followup} for individuals.
#' @param method The method used to obtain the average run length. Either "MC" for Markov Chain
#' or "IntEq" for integral equation. Default = "MC".
#' @param smooth_prob Should the probability distribution of failure under the null distribution be smoothed?
#' Useful for small samples. Can only be TRUE when \code{glmmod} is supplied. Default = FALSE.
#'
#'
#'
#' @return A list containing:
#' \itemize{
#' \item \code{ARL_0}: A numeric value indicating the average run length in
#' number of outcomes when starting from state E_0.
#' \item \code{ARL}: A \code{data.frame} containing the average run length (#outcomes)
#' depending on the state in which the process starts (E_0, E_1, ..., E_{n_grid-1})
#' \describe{
#'   \item{\code{start_val}:}{Starting value of the CUSUM, corresponding to the
#'    discretized state spaces E_{i};}
#'   \item{\code{#outcomes}:}{ARL for the CUSUM with
#'   initial value \code{start_val};}
#' }
#' \item \code{R}: A transition probability \code{matrix} containing the transition
#' probabilities between states \eqn{E_0, \ldots, E_{t-1}}{E_0, ..., E_{n_grid-1}}.
#' \eqn{R_{i,j}}{R_{i,j}} is the transition probability from state i to state j.
#' \item \code{h}: Value of the control limit.
#' } The value of \code{ARL_0} will be printed to the console.
#'
#' @references Brook, D., & Evans, D. A. (1972). An Approach to the Probability
#' Distribution of Cusum Run Length. Biometrika, 59(3), 539–549.
#' \doi{10.2307/2334805}
#'
#' Steiner, S. H., Cook, R. J., Farewell, V. T., & Treasure, T. (2000).
#' Monitoring surgical performance using risk-adjusted cumulative sum charts.
#' Biostatistics, 1(4), 441–452. \doi{10.1093/biostatistics/1.4.441}
#'
#' Kemp, K. W. (1971). Formal Expressions which Can Be Applied to Cusum Charts.
#' Journal of the Royal Statistical Society. Series B (Methodological), 33(3),
#' 331–360. \doi{10.1111/j.2517-6161.1971.tb01521.x}
#'
#'
#' @details The average run length of a CUSUM chart \eqn{S_n}{S_n} is given by
#' \eqn{\mathbb{E}[\tau_n],}{E[\tau_n],} where \eqn{\tau_n}{\tau_n} is defined as:
#' \deqn{\tau_n = \inf\{n \geq 0: S_n \geq h\}.}{\tau_n = inf(n >= 0: S_n >= h).}
#'
#' When \code{method = "MC"}, the average run length will be determined by
#' the Markov Chain approach described
#' in Brook & Evans (1972), using the risk-adjustment correction proposed in
#' Steiner et al. (2000). The idea is to discretize the domain (0, h) into \eqn{n_{grid} -1}{n_grid-1}
#' state spaces, with \eqn{E_0}{E_0} of width \eqn{w/2}{w/2}
#' and \eqn{E_1, \ldots, E_{n_{grid}-1}}{E_1, ..., E_{n_grid-1}} of width \eqn{w}{w}, such that
#' \eqn{E_{n_{grid}}}{E_{n_grid}} is an absorbing state.  This is done using the following steps:
#' \itemize{
#' \item \eqn{w}{w} is determined using the relationship \eqn{\frac{2h}{2t-1}}{2h/(2t-1)}.
#' \item Transition probabilities between the states are determined and
#' 'transition matrix' \eqn{R}{R} is constructed.
#' \item The equation \eqn{(\bm{I}-\bm{R}) \bm{ARL} = \bm{1}}{(I-R) ARL = 1} is
#' solved to find the ARL starting from each of the states.
#' }
#'
#' When \code{method = "IntEq"}, the average run length will be determined by
#' the integral equation approach described in Kemp (1971), using the risk-adjustment
#' correction proposed in Steiner et al. (2000). The idea is to exploit the
#' connection between the run length of a Sequential Probability Ratio Test (SPRT)
#' and that of a CUSUM. If N is the run length of a SPRT, P(0) the probability of
#' a SPRT terminating on the lower boundary of zero and R the run length of
#' a CUSUM, then: \deqn{\mathbb{E}[R] = \frac{\mathbb{E}[N]}{1 - P(0)}.}{E[R] = E[N]/(1-P(0)).}
#' \eqn{\mathbb{E}[N]}{E[N]} and \eqn{P(0)}{P(0)} are completely determined by
#' \deqn{G_n(z) = \int_0^h F(z-w) dG_{n-1}(w)}{G_n(z) = \int_0^h F(z-w) dG_{n-1}(w)}
#' with \eqn{F(x)}{F(x)} the cdf of the singletons \eqn{W_n}{Wn}. The integral can be
#' approximated using the generalized trapezoidal quadrature rule:
#' \deqn{G_n(z) = \sum_{i=0}^{n_{grid}-1} \frac{F(z-x_{i+1}) + F(z-x_{i})}{2} \left(G_{n-1}(x_{i+1}) - G_{n-1}(x_{i})  \right)}{Gn(z) = \sum_{i=0}^{n_grid-1} (F(z-x_{i+1}) + F(z-x_{i}))/2 * (G_{n-1}(x_{i+1}) - G_{n-1}(x_{i}))}
#'
#'
#'
#' @export
#'
#' @author Daniel Gomon
#' @family average run length
#' @seealso \code{\link[success]{bernoulli_cusum}}, \code{\link[success]{bernoulli_control_limit}}
#'
#'
#' @examples
#' #Determine a risk-adjustment model using a generalized linear model.
#' #Outcome (failure within 100 days) is regressed on the available covariates:
#' glmmodber <- glm((survtime <= 100) & (censorid == 1)~ age + sex + BMI,
#'                   data = surgerydat, family = binomial(link = "logit"))
#' #Determine the Average Run Length in number of outcomes for
#' #control limit h = 2.5 with (0, h) divided into n_grid = 200 segments
#' ARL <- bernoulli_ARL(h = 2.5, n_grid = 200, glmmod = glmmodber, theta = log(2))
#' #Calculate ARL, but now exploiting connection between SPRT and CUSUM:
#' #n_grid now decides the accuracy of the Trapezoidal rule for integral approximation
#' ARLIntEq <- bernoulli_ARL(h = 2.5, n_grid = 200, glmmod = glmmodber,
#' theta = log(2), method = "IntEq")
#'
#' \donttest{
#' #We can compare our ARL with that determined using the VLAD package
#' #See \url{https://cran.r-project.org/package=vlad}
#' if(require("vlad")){
#'    fi <- as.numeric(table(glmmodber$fitted.values)/length(glmmodber$fitted.values))
#'    pi1 <- sort(unique(glmmodber$fitted.values))
#'    pmix1 <- data.frame(fi, pi1, pi1)
#'    vlad_ARL <- round(vlad::racusum_arl_mc(pmix = pmix1, RA = 2, RQ = 1, h = 2.5, scaling = 200))
#' }
#' }




bernoulli_ARL <- function(h, n_grid, glmmod, theta, theta_true, p0, p1,
                          method = c("MC", "IntEq"), smooth_prob = FALSE){
  #------------Variable checks------------------
  method <- match.arg(method)
  #Input checks
  if(!all(is.numeric(h) & length(h) == 1)){
    stop("Parameter 'h' must be a single numeric variable.")
  }
  if(!all(is.numeric(n_grid), n_grid%%1 == 0, n_grid > 1, length(n_grid) == 1)){
    stop("Parameter 'n_grid' must be an integer larger than 1.")
  }
  if(!all(is.logical(smooth_prob), length(smooth_prob) == 1)){
    stop("Parameter 'smooth_prob' must be a single logical value.")
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
  if(!((method == "MC") | (method == "IntEq"))){
    stop("Parameter 'method' must be either 'MC' or 'IntEq'.")
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

  #############################CALCULATE CDF of W_n############################################
  Wncdf <- calc_Wncdf(glmmod = glmmod, theta = theta, theta_true = theta_true, p0 = p0, smooth_prob = smooth_prob)

  #Determine Average Run Length depending on method
  if(method == "MC"){
    R <- calc_MC_trans_matrix(h = h, n_grid = n_grid, Wncdf = Wncdf, glmmod = glmmod, p0 = p0, theta = theta, theta_true = theta_true)
    return(c(bernoulli_ARL_MC(n_grid = n_grid, R = R, h = h),
                h = h))
  } else if(method == "IntEq"){
    return(c(bernoulli_ARL_IntEq(h = h, n_grid = n_grid, Wncdf = Wncdf, glmmod = glmmod, p0 = p0, theta = theta, theta_true = theta_true),
                h = h))
  }
}



#' Cumulative distribution function (cdf) of Run Length for Bernoulli CUSUM
#'
#' @description Calculate the cdf of the Run Length of the Bernoulli CUSUM,
#' starting from initial value between 0 and \code{h}.
#'
#'
#' @inheritParams bernoulli_ARL
#' @param x Quantile at which to evaluate the cdf.
#' @param exact Should the cdf be determined exactly (TRUE), or approximately
#' (FALSE)? The approximation works well for large \code{x}, and can cut computation
#' time significantly. Default = TRUE.
#'
#'
#' @details Let \eqn{K}{K} denote the run length of the Bernoulli CUSUM with control limit \code{h}, then
#' this function can be used to evaluate \eqn{\mathbb{P}(K \leq x)}{P(K <= x)}.
#'
#' When \code{method = "MC"}, the formula on page 543 of Brook & Evans (1972)
#' is used if \code{exact = TRUE}. When \code{exact = FALSE}, formula (3.9) on
#' page 545 is used instead, approximating the transition matrix using its
#' Jordan canonical form. This can save computation time considerably, but is
#' not appropriate for small values of \code{x}.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{Fr_0}: A numeric value indicating the probability of the run
#' length being smaller than \code{x}.
#' \item \code{Fr}: A \code{data.frame} containing the cumulative distribution function of the run length
#' depending on the state in which the process starts (E_0, E_1, ..., E_{n_grid-1})
#' \describe{
#'   \item{\code{start_val}:}{Starting value of the CUSUM, corresponding to the
#'    discretized state spaces E_{i};}
#'   \item{\code{P(K <= x)}:}{Value of the cdf at \code{x} for the CUSUM with
#'   initial value \code{start_val};}
#' }
#' \item \code{R}: A transition probability \code{matrix} containing the transition
#' probabilities between states \eqn{E_0, \ldots, E_{t-1}}{E_0, ..., E_{n_grid-1}}.
#' \eqn{R_{i,j}}{R_{i,j}} is the transition probability from state i to state j.
#' } The value of \code{ARL_0} will be printed to the console.
#'
#'
#' @export
#'
#' @references Brook, D., & Evans, D. A. (1972). An Approach to the Probability
#' Distribution of Cusum Run Length. Biometrika, 59(3), 539–549.
#' \doi{10.2307/2334805}
#'
#' Steiner, S. H., Cook, R. J., Farewell, V. T., & Treasure, T. (2000).
#' Monitoring surgical performance using risk-adjusted cumulative sum charts.
#' Biostatistics, 1(4), 441–452. \doi{10.1093/biostatistics/1.4.441}
#'
#' Kemp, K. W. (1971). Formal Expressions which Can Be Applied to Cusum Charts.
#' Journal of the Royal Statistical Society. Series B (Methodological), 33(3),
#' 331–360. \doi{10.1111/j.2517-6161.1971.tb01521.x}
#'
#' @examples
#' #Determine a risk-adjustment model using a generalized linear model.
#' #Outcome (failure within 100 days) is regressed on the available covariates:
#' glmmodber <- glm((survtime <= 100) & (censorid == 1)~ age + sex + BMI,
#'                   data = surgerydat, family = binomial(link = "logit"))
#' #Determine probability of run length being less than 600
#' prob600 <- bernoulli_RL_cdf(h = 2.5, x = 600, n_grid = 200, glmmod = glmmodber, theta = log(2))


bernoulli_RL_cdf <- function(h, x, n_grid, glmmod, theta, theta_true, p0, p1,
                             method = c("MC", "IntEq"), smooth_prob = FALSE, exact = TRUE){
  #------------Variable checks------------------
  method <- match.arg(method)
  #Input checks
  if(!all(is.numeric(h) & length(h) == 1)){
    stop("Parameter 'h' must be a single numeric variable.")
  }
  if(!all(is.numeric(n_grid), n_grid%%1 == 0, n_grid > 1, length(n_grid) == 1)){
    stop("Parameter 'n_grid' must be an integer larger than 1.")
  }
  if(!all(is.logical(smooth_prob), length(smooth_prob) == 1)){
    stop("Parameter 'smooth_prob' must be a single logical value.")
  }
  if(!all(is.logical(exact), length(exact) == 1)){
    stop("Parameter 'exact' must be a single logical value.")
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
  if(!((method == "MC") | (method == "IntEq"))){
    stop("Parameter 'method' must be either 'MC' or 'IntEq'.")
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

  #############################CALCULATE CDF of W_n############################################
  Wncdf <- calc_Wncdf(glmmod = glmmod, theta = theta, theta_true = theta_true, p0 = p0, smooth_prob = smooth_prob)

  #Determine Average Run Length depending on method
  if(method == "MC"){
    R <- calc_MC_trans_matrix(h = h, n_grid = n_grid, Wncdf = Wncdf, glmmod = glmmod, p0 = p0, theta = theta, theta_true = theta_true)
    return(bernoulli_cdf_MC(n_grid = n_grid, R = R, r = x, h = h, exact = exact))
  } else if(method == "IntEq"){
    stop("Not implemented yet")
  }
}



#############################UTIL FUNCTIONS BELOW###############################
##################################INTERNAL######################################


#' Calculate cdf of singletons W_n for CUSUM
#'
#' @description Internal function to calculate cdf of singletons \eqn{W_n}{Wn}
#' of the Bernoulli CUSUM chart. The cdf is used to create the transition matrix
#' when Markov Chain methodology is used or to determine the integral equation/probabilities
#' of a Wald test when integral equation or Kemp's methodology is used.
#'
#' @inheritParams bernoulli_ARL
#'
#' @importFrom Rfast binary_search
#'
#' @keywords internal
calc_Wncdf <- function(glmmod, theta, theta_true, p0, smooth_prob = FALSE){
  #If we do not want to smooth probability distribution
  if(isFALSE(smooth_prob)){
    if(!missing(glmmod)){
      #IMPORTANT:
      #null_probs = sort(glmmod$fitted.values)
      if(theta > 0){
        #If we have glmmod, we can do a risk-adjusted calculation of the ARL.
        Wncdf <- function(x, null_probs, theta, theta_true = NULL){
          #Note that null_probs must be sort(glmmod$fitted.values)!!!!
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
            return(Wncdf(0, null_probs = null_probs, theta = theta, theta_true = theta_true) + sum(adjnull_probs[ind])*1/N)
          }
        }
      } else if(theta < 0){
        #If we have glmmod, we can do a risk-adjusted calculation of the ARL.
        Wncdf <- function(x, null_probs, theta, theta_true = NULL){
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
            return(Wncdf(0, null_probs = null_probs, theta = theta, theta_true = theta_true) + sum((1-adjnull_probs[ind]))*1/N)
          }
        }
      }
    }else{
      #if no glmmod, we don't perform risk-adjustment.
      #IMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANT
      #We use null_probs to represent p0
      #IMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANTIMPORTANT
      if(theta > 0){
        Wncdf <- function(x, null_probs, theta, theta_true = NULL){
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
        Wncdf <- function(x, null_probs, theta, theta_true = NULL){
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
  } else if(isTRUE(smooth_prob)){
    #smooth_prob is only usefull if glmmod is supplied.
    if(missing(glmmod)){
      stop("Probabilities can only be smoothed when 'glmmod' is supplied.")
    }
    #density_estimate <- density(glmmod$fitted.values)
    stop("smooth_prob not incorporated yet")
    #Parametric fit of p_0 distribution
    #Hein: user specifies distribution for linear predictor -> Normal(mu, sigma^2)
    #Then integrate from logit scale.
  }
  return(Wncdf)
}



#' Transition probability matrix for Bernoulli CUSUM
#'
#' @description Calculates the transition probability matrix for the Bernoulli
#' CUSUM described in Brook & Evans (1972).
#'
#' @inheritParams bernoulli_ARL
#'
#' @keywords internal
#'
#'

calc_MC_trans_matrix <- function(h, n_grid, Wncdf, glmmod, p0, theta, theta_true){
  #####################CALCULATE TRANSITION PROBABILITIES#######################################
  #Calculate the value w used for discretizing the state space
  w <- 2*h/(2*n_grid -1)
  #Initialize transition matrix
  R <- matrix(0, nrow = n_grid, ncol = n_grid)
  #Using Wncdf we can calculate the transition probabilities
  trans_prob <- function(state1, state2, w, t, null_probs, theta, theta_true = NULL){
    if(state2 == 0){
      #P(W_n \leq -i2 + 1/2*w)
      Wncdf((0.5-state1) * w, null_probs = null_probs, theta = theta, theta_true = theta_true)
    } else if(state2 == t){
      1- Wncdf((t-state1 - 0.5) * w, null_probs = null_probs, theta = theta, theta_true = theta_true)
    }else{
      Wncdf((state2-state1 + 0.5) * w, null_probs = null_probs, theta = theta, theta_true = theta_true) - Wncdf((state2-state1 - 0.5) * w, null_probs = null_probs, theta = theta, theta_true = theta_true)
    }
  }
  #We only need to calculate some specific transition probabilities:
  if(missing(glmmod)){
    #To avoid code duplication, define null_probs as p0 here.
    null_probs = p0
  } else{
    null_probs = sort(glmmod$fitted.values)
  }


  #First we calculate the probabilities in the first row and first column of R
  #For this we create the vector (p_{t-1, 0}, p_{t-2, 0}, ..., p_{0,0}, ..., p_{0, t-2}, p_{0, t-1})
  probstmin10 <- sapply((n_grid-1):0, function(x) trans_prob(x, 0, w = w, t = n_grid, null_probs = null_probs, theta = theta, theta_true = theta_true))
  probs1tmin1 <- sapply(1:(n_grid-1), function(x) trans_prob(0, x, w = w, t = n_grid, null_probs = null_probs, theta = theta, theta_true = theta_true))
  probs0 <- c(probstmin10, probs1tmin1)

  #Then we calculate the probabilities on the diagonal of the rest of the matrix:
  #These are given by the vector (p_{t-1, 1}, ..., p_{1, 1}, p_{1, 2}, ..., p_{1, t-1})
  probstmin11 <- sapply((n_grid-1):1, function(x) trans_prob(x, 1, w = w, t = n_grid, null_probs = null_probs, theta = theta, theta_true = theta_true))
  probs2tmin1 <- sapply(2:(n_grid-1), function(x) trans_prob(1, x, w = w, t = n_grid, null_probs = null_probs, theta = theta, theta_true = theta_true))
  probs1 <- c(probstmin11, probs2tmin1)

  #Iterate through matrix to fill values. For details, see Overleaf file (picture displaying the matrix).
  #We only need to calculate the first row and column values
  #And the values on the diagonals of the matrix without first row/column
  #Then we fill the matrix by checking on which diagonal we are (value of j-i)
  for(i in 1:n_grid){
    for(j in 1:n_grid){
      #If we are in the first row or column, we use probs0 vector
      if(i == 1| j == 1){
        R[i,j] <- probs0[n_grid + (j-i)]
      } else{
        #In any other row/column we use the probs1 vector
        R[i,j] <- probs1[(n_grid-1) + (j-i)]
      }
    }
  }
  return(R)
}




#' Average run length for Bernoulli CUSUM using Markov Chain methodology
#'
#' @description Internal function that discretizes grid and solves
#' matrix equation involving transition matrix for Markov Chain methodology
#'
#' @inheritParams bernoulli_ARL
#' @param R Transition probability matrix obtained from \code{calc_MC_trans_matrix}
#'
#' @importFrom Matrix solve
#'
#' @keywords internal
#'

bernoulli_ARL_MC <- function(n_grid, R, h){
  #############################DECOMPOSE TRANSITION MATRIX FOR ARL#######################
  #Now we need to decompose the "transition matrix" R to obtain the ARL vector mu
  #Each entry mu[i] represents the ARL when starting from state i-1
  #First entry is therefore the average run length when starting from 0.
  mu = round(Matrix::solve(diag(1, nrow = n_grid, ncol = n_grid) - R, matrix(1, nrow = n_grid, ncol = 1)))
  rownames(mu) <- 0:(n_grid-1)

  #Calculate value of w (grid discretization size)
  w <- 2*h/(2*n_grid -1)

  #Post-processing
  ARL_0 = mu[1]
  names(ARL_0) <- "#outcomes"
  ARL <- cbind(seq(0, h, w), mu)
  colnames(ARL) <- c("start_val", "#outcomes")
  rownames(ARL) <- 0:(n_grid -1)

  print(ARL_0)
  return(invisible(list(ARL_0 = ARL_0,
                        ARL = ARL,
                        R = R)))
}

#' Average run length for Bernoulli CUSUM using Markov Chain methodology
#'
#' @description Internal function that discretizes grid and solves
#' matrix equation involving transition matrix for Markov Chain methodology
#'
#' @inheritParams bernoulli_ARL
#' @param R Transition probability matrix obtained from \code{calc_MC_trans_matrix}
#'
#' @importFrom matrixcalc matrix.power
#'
#' @keywords internal
#'

bernoulli_cdf_MC <- function(n_grid, R, r, h, exact = TRUE){
  #Calculate value of w (grid discretization)
  w <- 2*h/(2*n_grid -1)
  #############################DECOMPOSE TRANSITION MATRIX FOR cdf#######################
  #According to Brook & Evans, the cdf is given by F(X_0 <= r, X_1 <= r, ..., X_{h-1} <= r) =  (I-R^r)1
  #First entry is therefore the cdf starting from 0
  if(isTRUE(exact)){
    Rr <- matrixcalc::matrix.power(R, r-1)
    Fr <- 1 - (Rr %*% matrix(1, nrow = n_grid, ncol = 1))
  } else{
    ######################Right-EV#############################
    #R = VLV^{-1}, so that RV = VL   -> Left eigenvector
    #Eigen decompose R
    eigenR <- eigen(R)
    #Check which eigenvalues are Real
    RealEV <- which(abs(Im(eigenR$values)) < 1e-8)
    #Determine largest real eigenvalue
    RealEVR <- Re(eigenR$values[RealEV])[1]
    #Determine associated eigenvector
    RealEVec <- Re(eigenR$vectors[,RealEV[1]])
    RealEVec <- sign(RealEVec[1]) * RealEVec

    ######################Left-EV#############################
    #R^T = VLV^{-1}, so that R^TV = VL  -> Left eigenvector
    LefteigenR <- eigen(t(R))
    #Check which eigenvalues are Real
    LeftRealEV <- which(abs(Im(LefteigenR$values)) < 1e-8)
    #Determine largest real eigenvalue
    LeftRealEVR <- Re(LefteigenR$values[LeftRealEV])[1]
    #Determine associated eigenvector
    LeftRealEVec <- Re(LefteigenR$vectors[,LeftRealEV[1]])
    LeftRealEVec <- sign(LeftRealEVec[1]) * LeftRealEVec

    Fr <- 1 - ((RealEVR^(r-1)) * (sum(LeftRealEVec))/(sum(RealEVec * LeftRealEVec))) * RealEVec
  }

  Fr_0 <- Fr[1]
  Fr <- cbind(seq(0, h, w), Fr)
  colnames(Fr) <- c("start_val", "P(K <= x)")
  rownames(Fr) <- 0:(n_grid -1)

  if(is.unsorted(Fr[, 2])){
    warning("Probabilities are not increasing with starting value, resulting probabilities may not be accurate.
            Consider increasing 'x' and/or 'ngrid'.")
  } else if(any(sign(Fr[,2]) < 0)){
    warning("Probability below 0 calculated, resulting probabilities may not be accurate.
            Consider increasing 'x' and/or 'ngrid'.")
  }
  print(Fr_0)
  return(invisible(list(Fr_0 = Fr_0,
                        Fr = Fr,
                        R = R)))
}


#' Average run length for Bernoulli CUSUM using Integral Equation methodology
#'
#' @description Internal function that calculates the ARL using the connection
#' between the ARL of a Wald SPRT and a CUSUM.
#'
#' @inheritParams bernoulli_ARL
#' @param Wncdf A function returning the values of the (risk-adjusted) cumulative
#' distribution function (cdf) for the singletons Wn.
#'
#' @keywords internal
#'

bernoulli_ARL_IntEq <- function(h, n_grid, Wncdf, glmmod, theta, theta_true, p0, tol = 1e-6){
  #TO-DO: Initialize cdf of W_n as in the function above.
  #Wncdf <- function
  #For now we check for normal:
  stieltjes_trapezoid <- function(fvals, gvals, n){
    #Stieltjes trapezoidal integration approximation
    #Inputs are:
    #fvals (Wncdf evaluated at relevant points z - [0,h])
    #gvals (G_{n-1}(w) evaluated at the grid of [0,h])
    #We want to calculate sum_{i=1}^{n_grid -1} (f(x_{i=1}) + f(x_i))/2 * (g(x_{i+1}) - g(x_{i}))
    #The terms (f(x_{i=1}) + f(x_i)) and (g(x_{i+1}) - g(x_{i})) can be calculated by performing
    #a single vector operation below. We then retrieve the sum and multiply by 0.5
    0.5*sum((fvals[-(n+1)] + fvals[-1]) * (gvals[-1] - gvals[-(n+1)]))
    #OLD - SLOW IMPLEMENTATION
    #Calculate the summation terms (f(x_i)+f(x_{i+1}))/2 * (G(x_i+1) - G(x_i))
    #approxvals <- sapply(1:n, function(x) (fvals[x] + fvals[x+1])/2 * (gvals[x+1] - gvals[x]))
    #Return sum as approximator of integral
    #sum(approxvals)
  }

  if(missing(glmmod)){
    #To avoid code duplication, define null_probs as p0 here.
    null_probs = p0
  } else{
    null_probs = sort(glmmod$fitted.values)
  }
  #Initiate 3 column matrix for storing G_N(0) in the first column,
  #G_N(h) in the second and G_N(Inf) in the third
  G_storage <- matrix(NA, nrow = 1, ncol = 3)
  #Initialize grid for Stieltjes integration
  #The grid is given by: [0, h/n_grid, 2*h/n_grid, ..., h]
  d <- (h-0)/n_grid
  xj <- (0:n_grid)*d
  #Initialize extended grid for cdf memoization
  xj_extended <- -h + (0:(2*n_grid))*d
  #Pre-calculate cdf on required values:
  #We calculate Wncdf on [-h, -h + h/ngrid, ..., 0, h/ngrid, ..., h]
  #Because we only ever integrate the cdf over these values!
  Fvals <- sapply(xj_extended, function(x) Wncdf(x, null_probs = null_probs, theta = theta, theta_true = theta_true))
  #Pre-calculate G_1(x) on [0,h]
  G_current <- Fvals[(n_grid+1):(2*n_grid + 1)]
  #Produce first output (G_1(0), G_1(h), G_1(Inf))
  #Note that G_1(Inf) = 1
  G_storage[1, ] <- c(G_current[1], G_current[n_grid + 1], 1)
  i <- 1
  #We continue calculating values until the tolerance is reached
  #Tolerance is defined as G_n(h) - G_n(0), the probability of a Wald test not stopping in N steps.
  #When this is sufficiently small, we barely lose any accuracy
  while((G_storage[i,2] - G_storage[i,1]) > tol){
    #Calculate values of G_{i}(z) over the grid
    #Remember that we have to integrate
    #G_n(z) = \int_0^h Wncdf(z-w) dG_{n-1}(w)
    #Let 0 = x_0 < x_1 < ..., x_{n_grid} = h be an equally spaced grid of [0,h]
    #We approximate the integral using the trapezoidal rule:
    #G_n(z) = \sum_{i = 1}^{ngrid-1} (Wncdf(z - x_{i+1}) + Wncdf(z - x_{i}))/2 * (G_{n-1}(x_{i+1}) - G_{n-1}(x_i))
    #Seeing as we are only interested in G_n(z) with z \in [0, h] we only
    #need to evaluate Wncdf on the grid [-h, h]. We use memoization to pre-calculate
    #Wncdf on a grid -h = s_0 < s_1 < ... < s_{2n_grid + 1} = h
    #This grid is given by xj_extended above. And the memoization is given by Fvals.
    G_new <- sapply(1:(n_grid + 1), function(x) stieltjes_trapezoid(Fvals[(n_grid + x):x], G_current, n = n_grid))
    G_inf <- stieltjes_trapezoid(rep(1, n_grid+1), G_current, n = n_grid)
    G_storage <- rbind(G_storage, c(G_new[1], G_new[n_grid + 1], G_inf))
    G_current <- G_new
    i <- i + 1
  }
  #E[R; 0] = (sum_{N=1}^Inf (G_N(Inf) - G_N(h) + G_N(0))*N)/(1-sum_{N=1}^Inf G_N(0))
  #Equation (11) in Kemp(1971)
  ARL <- round(sum((rowSums(G_storage) - G_storage[,2]*2) * 1:nrow(G_storage))/(1-sum(G_storage[,1])))
  names(ARL) <- "#outcomes"
  print(ARL)
  return(invisible(list(ARL_0 = ARL,
                        G = G_storage)))
}





