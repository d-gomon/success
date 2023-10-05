
test_that("Input checks",{
  exprfitfunnel <- (survtime <= 100) & (censorid == 1) ~ age + sex + BMI
  glmmodfun <- glm(exprfitfunnel, data = surgerydat, family = binomial(link = "logit"))
  expect_error(funnel_plot(ctime = 3*365, glmmod = glmmodfun, followup = 100),
               "Please provide data to construct chart.")
  expect_error(funnel_plot(data = surgerydat, ctime = 3*365, glmmod = glmmodfun, followup = 100, predlim = "asd"),
               "Argument predlim must be numeric vector with values between 0 and 1.")
  expect_error(funnel_plot(data = surgerydat, ctime = 3*365, glmmod = glmmodfun, followup = "asd"),
               "Argument followup must be a single numeric value larger than 0.")
  #Expect warning when p0 is not specified.
  expect_warning(funnel_plot(data = surgerydat, ctime = 3*365, glmmod = glmmodfun, followup = 100),
                 NULL)
})


test_that("Output checks", {
  exprfitfunnel <- (survtime <= 100) & (censorid == 1) ~ age + sex + BMI
  glmmodfun <- glm(exprfitfunnel, data = surgerydat, family = binomial(link = "logit"))
  suppressWarnings(funnel <- funnel_plot(data = surgerydat, ctime = 3*365, glmmod = glmmodfun, followup = 100))
  #Expect each unit in data to be present in funnel plot
  expect_equal(length(unique(surgerydat$unit)), length(unique(funnel$data$unit)))
  #Expect lower prediction interval bounds to be smaller than upper
  expect_true(all(funnel$plotdata$lower < funnel$plotdata$upper))
}
)

test_that("Internal checks", {
  pars <- parameter_assist(baseline_data = surgerydat,
                           data = subset(surgerydat, unit == 1),
                           formula = formula("survtime ~ age + sex + BMI"), followup = 100)
  assist_funnel <- funnel_plot(assist = pars)
  exprfitfunnel <- (survtime <= 100) & (censorid == 1) ~ age + sex + BMI
  suppressWarnings(funnel <- funnel_plot(surgerydat,
                        glmmod = glm(exprfitfunnel, data = surgerydat, family = binomial(link = "logit")),
                        followup = 100))
  expect_equal(assist_funnel$data, funnel$data)
})


