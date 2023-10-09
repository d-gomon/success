test_that("Check summary output", {
  exprfitfunnel <- (survtime <= 100) & (censorid == 1)~ age + sex + BMI
  glmmodfun <- glm(exprfitfunnel, data = surgerydat, family = binomial(link = "logit"))
  suppressWarnings(funnel <- funnel_plot(data = surgerydat, ctime = 3*365, glmmod = glmmodfun, followup = 100))
  expect_equal(funnel$data, summary(funnel))
})
