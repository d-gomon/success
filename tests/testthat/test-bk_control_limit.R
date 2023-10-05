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



