#
# #------------CHECKING CBASEH generation from coxphmod----------------
#
# exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
# tcoxmod <- coxph(exprfit, surgerydat)
#
# tyu <- basehaz(tcoxmod, centered = FALSE)
#
#
# cbase_temp <- basehaz(tcoxmod, centered = FALSE)
# cbase_loess <- loess(cbase_temp$hazard~cbase_temp$time)
# cbaseh <- function(t) predict(cbase_loess, t)
#
#
# plot(cbase_temp$time, cbase_temp$hazard, type = "l")
# plot(cbase_temp$time, exp(-cbase_temp$hazard), type = "l")
# xseq <- seq(1, 600, 1)
# plot(xseq, cbaseh(xseq))
#
# cbaseh_coxphmod <- function(coxphmod){
#   cbase_temp <- basehaz(coxphmod, centered = FALSE)
#   cbasefct <- function(t){
#     approxfun(x = cbase_temp$time, y = cbase_temp$hazard, )
#   }
# }
#
# cbase_temp <- basehaz(tcoxmod, centered = FALSE)
# cbaseh_coxphmod <- approxfun(x = cbase_temp$time, y = cbase_temp$hazard)
#
# RootLinearInterpolant <- function (x, y, y0 = 0) {
#   if (is.unsorted(x)) {
#     ind <- order(x)
#     x <- x[ind]; y <- y[ind]
#   }
#   z <- y - y0
#   ## which piecewise linear segment crosses zero?
#   k <- which(z[-1] * z[-length(z)] < 0)
#   ## analytically root finding
#   xk <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
#   xk
# }
#
# inv_cbaseh_temp <- function(y){
#   if(y <= max(cbase_temp$hazard)){
#     return(RootLinearInterpolant(x = cbase_temp$time, y = cbase_temp$hazard, y0 = y))
#   } else{
#     return(max(cbase_temp$time))
#   }
# }
#
# inv_cbaseh <- Vectorize(inv_cbaseh_temp)
#
# survtimes3 <- gen_surv_times(invchaz = inv_cbaseh, data = 10, coxphmod = coxphmod)
#
#
#
#
#
#
#
#
#
# #------------------------TEST RUNTIME OF CONTROL LIMIT----------------------
#
#
# library(tictoc)
#
# tic("INV_CBASEH")
# a <- cgr_control_limit(time = 500, alpha = 0.1, cbaseh = function(t) chaz_exp(t, lambda = 0.02),
#                        inv_cbaseh = function(t) inv_chaz_exp(t, lambda = 0.02), psi = 0.5, n_sim = 10)
# toc()
# #INV_CBASEH: 25.53 sec elapsed
#
# tic("coxmod + resampling")
# b <- cgr_control_limit(time = 500, alpha = 0.1, coxphmod = tcoxmod, psi = 0.5,
#                        baseline_data = surgerydat, pb = TRUE)
# toc()
# #coxmod + resampling: 77.02 sec elapsed
#
# tic("coxmod without resampling")
# b2 <- cgr_control_limit(time = 500, alpha = 0.1, coxphmod = tcoxmod, psi = 0.5, n_sim = 10)
# toc()
# #coxmod without resampling: 24.27 sec elapsed
#
# tic("cbaseh without resampling")
# c <- cgr_control_limit(time = 500, alpha = 0.1, cbaseh = function(t) chaz_exp(t, lambda = 0.02), psi = 0.5, n_sim = 10)
# toc()
# #cbaseh without resampling: 27.34 sec elapsed
#
#
# tic("BK with coxphmod + resampling")
# cbk <- bk_control_limit(time = 6*365, alpha = 0.1, coxphmod = tcoxmod, psi = 2,
#                         n_sim = 10, pb = TRUE, theta = log(2),
#                         baseline_data = surgerydat)
# toc()
#
# #Resampling seems to take a lot of time
#
#
# #------------LROI like test-------------
# tic("LROI + resampling + small hospital")
# b <- cgr_control_limit2(time = 6*365, alpha = 0.1, coxphmod = tcoxmod, psi = 0.5,
#                        baseline_data = surgerydat, pb = TRUE, chartpb = TRUE, ncores = 4)
# toc()
#
# #About 1h per 20 hospitals - 3min/hospital!!!! SLOW AF
# #With 2 cores: 1:30 per hospital. 30 min total So it scales with number of cores!
#
# tic("LROI + NO resampling + small hospital")
# b <- cgr_control_limit(time = 6*365, alpha = 0.1, coxphmod = tcoxmod, psi = 0.5,
#                        pb = TRUE, chartpb = TRUE, ncores = 4)
# toc()
#
# #About 30minutes. 1.5min/hospital.
# #2 cores: About 15 minutes 0.75min/hospital
#
#
# #-------------TEST WHETHER multi-core requires dependencies when using linear interpolant
#
# library(survival)
# exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
# tcoxmod <- coxph(exprfit, data= surgerydat)
# cbaseh1_temp <- basehaz(tcoxmod, centered = FALSE)
# cbaseh1 <- approxfun(x = cbaseh1_temp$time, y = cbaseh1_temp$hazard)
#
#
# aty <- cgr_cusum(data = subset(surgerydat, hosp_num < 15), cbaseh = cbaseh1, ncores = 3, pb = TRUE)
#
#
#
# #Testing generate_units
#
# b <- cgr_control_limit2(time = 500, alpha = 0.1, coxphmod = tcoxmod, psi = 0.5, n_sim = 5,
#                        baseline_data = surgerydat, pb = TRUE, chartpb = TRUE)
# c <- cgr_control_limit(time = 500, alpha = 0.1, coxphmod = tcoxmod, psi = 0.5, n_sim = 5,
#                         baseline_data = surgerydat, pb = TRUE, chartpb = TRUE)
#
# y <- generate_units(time = 500, psi = 0.5, coxphmod = tcoxmod, baseline_data = surgerydat)
#
#
#
#
# #Test behaviour of missingness and NULL
#
# tft <- function(t = NULL){
#   if(missing(t)){
#     print("asd")
#   }
#   if(is.null(t)){
#     print("hallo")
#   }
# }
#
#
#
#
# tft2 <- function(par1 = NULL, par2){
#   trt <- function(par1, par2){
#     if(missing(par1)){
#       print("par1 miss")
#     }
#     if(missing(par2)){
#       print("par2 miss")
#     }
#   }
#   trt(par1 = par1, par2 = par2)
# }
#
#
#
#
#
#
#
#
#
#
#
# inv_tcbaseh_temp <- function(y, lower = 0, upper = 9e12){
#   tryCatch({return(unname(unlist(uniroot((function (x) tcbaseh(x) - y),
#                                          lower = 0, upper = 9e12)$root)))},
#            error = function(cond){ return(upper)})
# }
#
#
#
# exprfit2 <- as.formula("Surv(survtime, censorid) ~ 1")
#
# tcoxmod2 <- coxph(exprfit2, surgerydat)
#
# exprfit3 <- as.formula("Surv(survtime, censorid) ~ NULL")
#
# tcoxmod3 <- coxph(exprfit3, surgerydat)
#
# calc_risk(data = surgerydat, coxphmod = tcoxmod2)
#
#
# exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
# tcoxmod <- coxph(exprfit, data= surgerydat)
#
# a <- calc_risk(data = surgerydat, coxphmod = tcoxmod)
#
# b <- unname(predict(tcoxmod, newdata = surgerydat, type = "risk", reference = "zero"))
# c <- unname(predict(tcoxmod3, newdata = surgerydat, type = "risk", reference = "zero"))
#
# library(tictoc)
# tic("predict")
# a <- replicate(1000,unname(predict(tcoxmod, newdata = surgerydat, type = "risk", reference = "zero")))
# toc()
#
# tic("calc_risk")
# b<- replicate(1000,calc_risk(data = surgerydat, coxphmod = tcoxmod))
# toc()
#
#
# followup <- 100
# #Determine a risk-adjustment model using a generalized linear model.
# #Outcome (failure within 100 days) is regressed on the available covariates:
# exprfitber <- as.formula("(survtime <= followup) & (censorid == 1)~ age + sex + BMI")
# glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
# tic("bernoulli")
# a <- bernoulli_control_limit(time = 500, alpha = 0.1, followup = followup,
#  psi = 0.5, n_sim = 100, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat)
# toc()
#





a <- cgr_control_limit(time = 50, alpha = 0.05, psi = 2, n_sim = 10, cbaseh = function(t) chaz_exp(t, lambda = 0.02), inv_cbaseh = function(t) inv_chaz_exp(t, lambda = 0.02))$h
b <- cgr_control_limit(time = 50, alpha = 0.05, psi = 2, n_sim = 10, cbaseh = function(t) chaz_exp(t, lambda = 0.02), inv_cbaseh = function(t) inv_chaz_exp(t, lambda = 0.02), maxtheta = Inf)$h
c <- cgr_control_limit(time = 50, alpha = 0.05, psi = 2, n_sim = 10, cbaseh = function(t) chaz_exp(t, lambda = 0.02), inv_cbaseh = function(t) inv_chaz_exp(t, lambda = 0.02), maxtheta = log(4))$h
d <- bk_control_limit(time = 50, alpha = 0.05, psi = 2, n_sim = 10, cbaseh = function(t) chaz_exp(t, lambda = 0.02), inv_cbaseh = function(t) inv_chaz_exp(t, lambda = 0.02), theta = log(2))$h
test_that("Decreasing maxthetat reduces control limit",{
  expect_lt(a, b)
  expect_lt(c, a)
})



