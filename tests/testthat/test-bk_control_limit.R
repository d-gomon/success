require(survival)

test_that("input checks", {
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  expect_error(bk_control_limit(time = -500, alpha = 0.1, theta = log(2),
                                coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat),
               "Argument time must be a single positive numeric value.")
  expect_error(bk_control_limit(time = 500, alpha = -0.1, theta = log(2),
                                coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat),
               "Argument alpha must be a single numeric value between 0 and 1.")
  expect_error(bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
                                coxphmod = tcoxmod, psi = -0.5, n_sim = 10, baseline_data = surgerydat),
               "Argument psi must be a single numeric value larger than 0.")
  expect_error(bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
                                coxphmod = tcoxmod, psi = 0.5, n_sim = 10.5, baseline_data = surgerydat),
               "Argument n_sim must be a single integer value larger than 0.")

})





test_that("Output checks", {
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  bkcontrol <- bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
                                coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat)
  expect_length(bkcontrol$charts, 10)
  bkcontrol2 <- bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
                                 coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat, seed = 1041996)
  expect_equal(bkcontrol2$h, bkcontrol$h)
  expect_output(bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
                                 coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat, pb = TRUE))
  #Control limit should be lower than 0 when lower-sided CUSUM
  expect_lt(bk_control_limit(time = 500, alpha = 0.1, theta = -log(2),
                             coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat)$h, 0)
})



test_that("Parameter_assist works as expected", {
  skip_on_cran()
  pars <- parameter_assist(baseline_data = surgerydat,
                           data = subset(surgerydat, unit == 1),
                           formula = formula("survtime ~ age + sex + BMI"), time = 500)
  bk_assist <- bk_control_limit(assist = pars)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  bk_noassist <- bk_control_limit(time = 500, alpha = 0.05, theta = log(2),
                                  coxphmod = tcoxmod, psi = arrival_rate(subset(surgerydat, unit == 1)), n_sim = 200, baseline_data = surgerydat)
  expect_equal(bk_assist$h, bk_noassist$h)
})


test_that("Theory", {
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  psilow <- bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
                             coxphmod = tcoxmod, psi = 0.3, n_sim = 10, baseline_data = surgerydat)
  psihigh <- bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
                              coxphmod = tcoxmod, psi = 0.7, n_sim = 10, baseline_data = surgerydat)
  #Increasing psi should increase control limit
  expect_lt(psilow$h, psihigh$h)
  ##############################################
  thetalow <- bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
                               coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat)
  thetahigh <- bk_control_limit(time = 500, alpha = 0.1, theta = log(6),
                                coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat)
  #Increasing theta should increase control limit
  expect_lt(thetalow$h, thetahigh$h)
  ##############################################
  alphalow <- bk_control_limit(time = 500, alpha = 0.1, theta = log(2),
                               coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat)
  alphahigh <- bk_control_limit(time = 500, alpha = 0.4, theta = log(2),
                                coxphmod = tcoxmod, psi = 0.5, n_sim = 10, baseline_data = surgerydat)
  #Increasing alpha decreases control limit
  expect_gt(alphalow$h, alphahigh$h)
  ##############################################
})



