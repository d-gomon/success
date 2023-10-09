require(survival)

test_that("input checks", {
  tdat <- subset(surgerydat, unit == 1)
  tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data = surgerydat)
  expect_error(cgr_cusum(coxphmod = tcoxmod, cbaseh = tcbaseh), "Please provide data to construct chart.")
  expect_error(cgr_cusum(data = tdat),
               "Please specify cbaseh (function) or coxphmod as Survival object.", fixed = TRUE)
  expect_error(cgr_cusum(data = tdat, cbaseh = tcbaseh,  h = c(3,4, 5)),
               "Please specify only 1 control limit.", fixed = TRUE)
  expect_error(cgr_cusum(data = tdat, cbaseh = tcbaseh, maxtheta = -0.3), "Parameter 'maxtheta' must be larger than 0.
      Expect at most a doubling (exp(theta) = 2) of cumulative hazard? theta = log(2)
      Expect at most a halving (exp(theta) = 0.5) of cumulative hazard rate? theta = log(2) and detection = 'lower'", fixed = TRUE)
  expect_error(cgr_cusum(data = tdat, cbaseh = tcbaseh, ctimes = seq(-10, -5, 1)),
               "Cannot construct chart before subjects enter into study (max(ctimes) <= min(data$entrytime)). Please re-asses the argument 'ctimes'.", fixed = TRUE)
  expect_error(cgr_cusum(data = tdat, cbaseh = tcbaseh, coxphmod = list(asd = c(3, 4), tyu = c(5,6))),
               "coxphmod does not contain $formula and/or $coefficients.", fixed = TRUE)
})


test_that("cores don't influence results",{
  skip()
  tdat <- subset(surgerydat, unit == 1 & entrytime < 100)
  tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  cgr1c <- cgr_cusum(data = tdat, coxphmod = tcoxmod)
  cgr3c <- cgr_cusum(data = tdat, coxphmod = tcoxmod, ncores = detectCores())
  expect_equal(cgr1c$CGR, cgr3c$CGR)
}
)

test_that("Specifying C reduces ctimes", {
  tdat <- subset(surgerydat, unit == 1 & entrytime < 100)
  tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  cgrwithoutC <- cgr_cusum(data = tdat, coxphmod = tcoxmod)
  cgrwithC <- cgr_cusum(data = tdat, coxphmod = tcoxmod, C = 50)
  expect_lt(nrow(cgrwithC$CGR), nrow(cgrwithoutC$CGR))
})


test_that("maxtheta is working as intended",{
  tdat <- subset(surgerydat, unit == 1 & entrytime < 100)
  tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  cgr1c <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh)
  maxtt <- log(6)
  cgr1clower <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, detection = "lower", maxtheta = maxtt)

  expect_equal(all(cgr1c$CGR$exp_theta_t <= exp(maxtt)), TRUE)
  expect_equal(all(cgr1clower$CGR$exp_theta_t >= exp(-maxtt)), TRUE)
})


test_that("CPU also works + progress bar", {
  tdat <- subset(surgerydat, unit == 1 & entrytime < 100)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  expect_no_error(cgr_cusum(data = tdat, coxphmod = tcoxmod, cmethod = "CPU"))
  expect_output(cgr_cusum(data = tdat, coxphmod = tcoxmod, cmethod = "CPU", pb = TRUE))
})


test_that("specifying control limit works", {
  tdat <- subset(surgerydat, unit == 1 & entrytime < 100)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  cgr <- cgr_cusum(data = tdat, coxphmod = tcoxmod, h = 4)
  expect_true(sum(cgr$CGR$value >= 4) <= 1)
  cgrlower <- cgr_cusum(data = tdat, coxphmod = tcoxmod, h = 2, detection = "lower")
  expect_true(sum(cgrlower$CGR$value <= -2) <= 1)
})


test_that("parameter assist works as expected", {
  #Specifying all parameters
  pars <- parameter_assist(baseline_data = surgerydat,
                           data = subset(surgerydat, unit == 1),
                           formula = formula("survtime ~ age + sex + BMI"))
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data = surgerydat)
  cgr <- cgr_cusum(assist = pars)
  cgr2 <- cgr_cusum(data = subset(surgerydat, unit == 1),
                  coxphmod = tcoxmod)
  expect_equal(cgr$CGR, cgr2$CGR)
})



#Known issues:

#When patients with survtime = 0 are present, the ML estimate cannot be calculated
#at the failure time of the patient.
#If no patients are present in dataset at that point in time, then at the time
#of entry of the next patient into the data set S_2, we can calculate CGR(t) at
#S_2 + \epsilon to obtain an arbitrarily large value of CGR(t). This is because
#thetat = log(NDT/AT) = log(1/epsilon) -> infty when epsilon -> 0.
#Demonstration:
#cgr <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, ctimes = seq(5, 30, 1))
#cgr2 <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, ctimes = seq(5, 30, 0.1))
#cgr3 <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, ctimes = seq(5, 30, 0.01))
#plot the CGR-CUSUM and note how the value can increase infinitely
#plot(cgr)
#plot(cgr2)
#plot(cgr3)
