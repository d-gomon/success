require(survival)

test_that("input checks", {
  tdat <- subset(surgerydat, unit == 1)
  tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data = surgerydat)
  expect_error(bk_cusum(theta = log(2),
                        coxphmod = tcoxmod, cbaseh = tcbaseh), "Please provide data to construct chart.")
  expect_error(bk_cusum(data = tdat, cbaseh = tcbaseh),
               "Please specify a value for theta (ln(expected hazard ratio)).", fixed = TRUE)
  expect_error(bk_cusum(data = tdat, theta = log(2)),
               "Please specify cbaseh (function) or coxphmod as Survival object.", fixed = TRUE)
  expect_error(bk_cusum(data = tdat, cbaseh = tcbaseh, theta = log(2), twosided = TRUE, h = c(3,4)),
               "When specifying 2 control limits the two values should have reverse signs.", fixed = TRUE)
  expect_error(bk_cusum(data = tdat, cbaseh = tcbaseh, theta = log(2), twosided = TRUE, h = c(3,4, 5)),
               "Please provide 1 or 2 values for the control limit.", fixed = TRUE)
  expect_error(bk_cusum(data = tdat, cbaseh = tcbaseh, theta = log(2),  h = c(3,4, 5)),
               "Please provide only 1 value for the control limit", fixed = TRUE)
})



test_that("output checks", {
  tdat <- subset(surgerydat, unit == 1)
  tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data = surgerydat)
  bkcus <- bk_cusum(data = tdat, theta = log(2),
                        coxphmod = tcoxmod, cbaseh = tcbaseh)
  bkcus2 <- bk_cusum(data = tdat, theta = log(2),
                     coxphmod = tcoxmod, cbaseh = tcbaseh, h = 3)
  bkcus3 <- bk_cusum(data = tdat, theta = log(2),
                     coxphmod = tcoxmod, cbaseh = tcbaseh, ctimes = seq(6, 100, 1))
  bkcus4 <- bk_cusum(data = tdat, theta = log(2),
                     coxphmod = tcoxmod, cbaseh = tcbaseh, C = 100)
  bksmaller3 <- which(bkcus$BK$value < 3)
  expect_equal(bkcus$BK$value[bksmaller3], bkcus2$BK$value[bksmaller3])
})

test_that("Automatic cbaseh determination", {
  tdat <- subset(surgerydat, unit == 1)
  tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data = surgerydat)
  bkcus <- bk_cusum(data = tdat, theta = log(2),
                    coxphmod = tcoxmod)
  bkcus2 <- bk_cusum(data = tdat, theta = log(2),
                     coxphmod = tcoxmod, cbaseh = extract_hazard(tcoxmod)$cbaseh)
  expect_equal(bkcus$BK, bkcus2$BK)
})

test_that("Lower sided CUSUM", {
  tdat <- subset(surgerydat, unit == 1)
  tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data = surgerydat)
  bkcus <- bk_cusum(data = tdat, theta = -log(2),
                    coxphmod = tcoxmod)
  expect_true(all(bkcus$BK$value <= 0))
})


test_that("Two-sided vs one-sided", {
  tdat <- subset(surgerydat, unit == 1)
  tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data = surgerydat)
  bkcus <- bk_cusum(data = tdat, theta = log(2),
                    coxphmod = tcoxmod)
  bkcus2 <- bk_cusum(data = tdat, theta = log(2),
                     coxphmod = tcoxmod, twosided = TRUE)
  expect_equal(bkcus2$BK$val_up, bkcus$BK$value)
  expect_output(bk_cusum(data = tdat, theta = log(2),
                         coxphmod = tcoxmod, pb = TRUE))
  bkcus3 <- bk_cusum(data = tdat, theta = log(2),
                     coxphmod = tcoxmod, twosided = TRUE, h = 3)
  bkcus4 <- bk_cusum(data = tdat, theta = log(2),
                     coxphmod = tcoxmod, twosided = TRUE, h = c(-3, 3))
  expect_equal(bkcus3$BK, bkcus4$BK)
})


test_that("parameter assist works as expected", {
  #Specifying all parameters
  pars <- parameter_assist(baseline_data = surgerydat,
                           data = subset(surgerydat, unit == 1),
                           formula = formula("survtime ~ age + sex + BMI"), followup = 100)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data = surgerydat)
  bk <- bk_cusum(assist = pars)
  bk2 <- bk_cusum(data = subset(surgerydat, unit == 1),
                  coxphmod = tcoxmod, theta = log(2))
  expect_equal(bk$BK, bk2$BK)
})

