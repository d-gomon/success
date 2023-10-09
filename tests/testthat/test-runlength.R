test_that("Input checks", {
  #Bernoulli CUSUM
  exprfitber <- (survtime <= 100) & (censorid == 1) ~ age + sex + BMI
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 14), glmmod = glmmodber,
                            followup = 100, theta = log(2))
  expect_error(runlength(bercus), NULL)
  expect_error(runlength(bercus, h = "asd"), NULL)
  expect_error(runlength.bercusum(h = 3), NULL)
  expect_error(runlength.bercusum(bercus, h = c(3,4)), NULL)


  #BK-CUSUM
  require(survival)
  tdat <- subset(surgerydat, unit == 1)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  bk <- bk_cusum(data = tdat, theta = log(2), coxphmod = tcoxmod)
  expect_error(runlength(bk), NULL)
  expect_error(runlength(bk, h = "asd"), NULL)
  expect_error(runlength.bkcusum(h = 3), NULL)
  expect_error(runlength.bkcusum(bk, h = c(3, 4)), NULL)

  #CGR-CUSUM
  tdat <- subset(surgerydat, unit == 1)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  cgr <- cgr_cusum(data = tdat, coxphmod = tcoxmod)
  expect_error(runlength(cgr), NULL)
  expect_error(runlength(cgr, h = "asd"), NULL)
  expect_error(runlength.cgrcusum(h = 3), NULL)
  expect_error(runlength.cgrcusum(cgr, h = c(3, 4)), NULL)
  expect_type(runlength(cgr, h = 2), "double")
  expect_equal(runlength(cgr, h = 100), Inf)
})


test_that("Output checks", {
  #Two-sided cusum control limit (Ber + BK - two control limits, CGR has no two sided)
  #

  #Bernoulli CUSUM
  exprfitber <- (survtime <= 100) & (censorid == 1) ~ age + sex + BMI
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 14), glmmod = glmmodber,
                            followup = 100, theta = log(2))
  expect_type(runlength(bercus, h = 1), "double")
  expect_equal(runlength(bercus, h = 100), Inf)
  bercus2 <- bernoulli_cusum(data = subset(surgerydat, unit == 14), glmmod = glmmodber,
                             followup = 100, theta = log(2), twosided = TRUE)
  expect_type(runlength(bercus2, h = 1), "double")
  expect_equal(runlength(bercus2, h = 100), Inf)
  expect_type(runlength(bercus2, h = c(-1, 1)), "double")
  expect_error(runlength(bercus2, h = c(1, 3)), NULL)
  expect_type(runlength(bercus2, h = -1), "double")
  expect_error(runlength(bercus2, h = c(-1, 1, 3)), NULL)

  #BK-CUSUM
  require(survival)
  tdat <- subset(surgerydat, unit == 1)
  exprfit <- Surv(survtime, censorid) ~ age + sex + BMI
  tcoxmod <- coxph(exprfit, data= surgerydat)
  bk <- bk_cusum(data = tdat, theta = log(2), coxphmod = tcoxmod)
  expect_type(runlength(bk, h = 1), "double")
  expect_equal(runlength(bk, h = 100), Inf)
  bk2 <- bk_cusum(data = subset(surgerydat, unit == 14), coxphmod = tcoxmod,
                             theta = log(2), twosided = TRUE)
  expect_type(runlength(bk2, h = 1), "double")
  expect_equal(runlength(bk2, h = 100), Inf)
  expect_type(runlength(bk2, h = c(-1, 1)), "double")
  expect_error(runlength(bk2, h = c(1, 3)), NULL)
  expect_type(runlength(bk2, h = -1), "double")
  expect_error(runlength(bk2, h = c(-1, 1, 3)), NULL)

})
