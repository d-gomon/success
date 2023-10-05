test_that("Input checks", {
  followup <- 100
  exprfitber <- (survtime <= followup) & (censorid == 1)~ age + sex + BMI
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  expect_error(bernoulli_cusum(glmmod = glmmodber,
                            followup = followup, theta = log(2)), "Please provide data to construct chart.")
  expect_error(bernoulli_cusum(data = surgerydat, glmmod = glmmodber,
                               followup = "asd", theta = log(2)),"Argument followup must be a single numeric value larger than 0.")
  expect_error(bernoulli_cusum(data = surgerydat, glmmod = glmmodber,
                               followup = followup, theta = log(2), twosided = TRUE, h = c(3,4)),
               "When specifying 2 control limits the two values should have reverse signs.")
  expect_error(bernoulli_cusum(data = surgerydat, glmmod = glmmodber,
                               followup = followup, theta = log(2), h = c(3,4)), "Please provide only 1 value for the control limit")
}
)




test_that("No data to construct cusum", {
  followup <- 0.2
  exprfitber <- (survtime <= followup) & (censorid == 1)~ age + sex + BMI
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  expect_warning(bernoulli_cusum(data = subset(surgerydat, unit == 1 & survtime > 1), glmmod = glmmodber,
                               followup = followup, theta = log(2), stoptime = 2), "No failures observed in specified time frame.
Decrease 'followup' or consider a larger time frame for construction.
Returning trivial chart.")
})



test_that("Two-sided bernoulli CUSUM vs one-sided", {
  followup <- 100
  exprfitber <- (survtime <= followup) & (censorid == 1)~ age + sex + BMI
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
                            followup = followup, theta = log(2))
  bercus2 <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
                             followup = followup, theta = log(2), twosided = TRUE)
  expect_equal(bercus$CUSUM$value, bercus2$CUSUM$val_up)
  bercus3 <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
                             followup = followup, theta = log(2), twosided = TRUE, h = 3)
  expect_equal(bercus3$h, c(-3, 3))
})



test_that("Control limit stops the chart + checks", {
  followup <- 100
  exprfitber <- (survtime <= followup) & (censorid == 1)~ age + sex + BMI
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
                            followup = followup, theta = log(2), h = 3)
  expect_true(sum(bercus$CUSUM$value >= 3) <= 1)
  bercus_lower <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
                                  followup = followup, theta = -log(2), h = 3)
  expect_equal(bercus$h, 3)
  expect_equal(bercus_lower$h, -3)
  bercus_two <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
                                followup = followup, theta = -log(2), h = c(3, -3), twosided = TRUE)
  expect_true(sum(bercus_two$CUSUM$val_up > 3) + sum(bercus_two$CUSUM$val_down < -3) <= 1)
})


test_that("Constructing with pre-specified p0 and theta/p1 + equality with theta specification", {
  followup <- 100
  exprfitber <- (survtime <= followup) & (censorid == 1)~ age + sex + BMI
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 1), p0 = 0.5, p1 = (0.5*2)/(1-0.5 + 0.5*2),
                            followup = followup, h = 3)
  bercus2 <- bernoulli_cusum(data = subset(surgerydat, unit == 1), p0 = 0.5, theta = log(2),
                             followup = followup, h = 3)
  expect_equal(bercus$CUSUM, bercus2$CUSUM)
})


