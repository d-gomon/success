test_that("input checks", {
  followup <- 100
  exprfitber <- (survtime <= followup) & (censorid == 1)~ age + sex + BMI
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  expect_error(bernoulli_control_limit(time = -500, alpha = 0.1, followup = followup,
                               psi = 0.5, n_sim = 10, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat),
               "Argument time must be a single positive numeric value.")
  expect_error(bernoulli_control_limit(time = 500, alpha = -0.1, followup = followup,
                                       psi = 0.5, n_sim = 10, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat),
               "Argument alpha must be a single numeric value between 0 and 1.")
  expect_error(bernoulli_control_limit(time = 500, alpha = 0.1, followup = "asd",
                                       psi = 0.5, n_sim = 10, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat),
               "Argument followup must be a single numeric value larger than 0.")
  expect_error(bernoulli_control_limit(time = 500, alpha = 0.1, followup = followup,
                                       psi = -0.5, n_sim = 10, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat),
               "Argument psi must be a single numeric value larger than 0.")
  expect_error(bernoulli_control_limit(time = 500, alpha = 0.1, followup = followup,
                                       psi = 0.5, n_sim = 10.5, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat),
               "Argument n_sim must be a single integer value larger than 0.")
  expect_error(bernoulli_control_limit(time = 99, alpha = 0.1, followup = followup,
                                       psi = 0.5, n_sim = 10, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat),
               "Argument followup must be greater than time, otherwise no events will be observed.")

})


test_that("Output checks", {
  followup <- 100
  exprfitber <- (survtime <= followup) & (censorid == 1)~ age + sex + BMI
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  bercontrol <- bernoulli_control_limit(time = 500, alpha = 0.1, followup = followup,
                          psi = 1, n_sim = 10, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat)
  expect_length(bercontrol$charts, 10)
  bercontrol2 <- bernoulli_control_limit(time = 500, alpha = 0.1, followup = followup,
                                         psi = 1, n_sim = 10, theta = log(2), glmmod = glmmodber, baseline_data = surgerydat, seed = 1041996)
  expect_equal(bercontrol2$h, bercontrol$h)
  expect_output(bernoulli_control_limit(time = 500, alpha = 0.1, followup = followup,
                                        psi = 1, n_sim = 10, theta = log(2), glmmod = glmmodber,
                                        baseline_data = surgerydat, seed = 1041996, pb = TRUE))
  #Control limit should be lower than 0 when lower-sided CUSUM
  expect_lt(bernoulli_control_limit(time = 500, alpha = 0.1, followup = followup,
                                    psi = 1, n_sim = 10, theta = -log(2), glmmod = glmmodber,
                                    baseline_data = surgerydat, seed = 1041996)$h, 0)
})
