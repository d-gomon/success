test_that("input checks", {
  expect_equal(calc_risk(NULL), 1)
  expect_equal(calc_risk(), 1)
  expect_equal(calc_risk(data.frame(asd = c())), 1)
  expect_equal(calc_risk(data.frame(asd = rnorm(10)), coxphmod = NULL), rep(1, 10))
})


test_that("Manual Risk-adjustment", {
  crdat <- data.frame(age = rnorm(10, 40, 5), BMI = rnorm(10, 24, 3))
  crlist <- list(formula = as.formula("~age + BMI"), coefficients = c("age"= 0.02, "BMI"= 0.009))
  out <- calc_risk(crdat, crlist)
  expect_type(out, "double")
  expect_length(out, nrow(crdat))
})




test_that("survival package vs manual RA", {
  crdat <- data.frame(age = rnorm(10, 40, 5), BMI = rnorm(10, 24, 3))
  phmod <- coxph(Surv(survtime, censorid) ~ age +  BMI, data = surgerydat)
  expect_length(calc_risk(crdat, phmod), nrow(crdat))
  ph_manual <- list(formula = ~ age + BMI,
                    coefficients = phmod$coefficients)
  expect_equal(unname(calc_risk(crdat, phmod)), calc_risk(crdat, ph_manual))
})
