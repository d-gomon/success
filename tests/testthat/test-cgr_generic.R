#Are plots created without errors for the examples in the package?

test_that(
  "Bernoulli CUSUM plot", {
  followup <- 100
  exprfitber <- as.formula("(survtime <= followup) & (censorid == 1)~ age + sex + BMI")
  glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
  bercus <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
                            followup = followup, theta = log(2))
  expect_no_error(plot(bercus))
  expect_no_error(plot(bercus, h = 3))
  bercus$h <- 3
  expect_no_error(plot(bercus))
  bercus2sided <- bernoulli_cusum(data = subset(surgerydat, unit == 1), glmmod = glmmodber,
                                  followup = followup, theta = log(2), twosided = TRUE)
  expect_no_error(plot(bercus2sided))
  expect_no_error(plot(bercus2sided, h = c(-3, 4)))
  bercus2sided$h <- c(-3, 4)
  expect_no_error(plot(bercus2sided))
  }
)

test_that(
  "BK-CUSUM plot", {
    require(survival)
    tdat <- subset(surgerydat, unit == 1)
    tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
    exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
    tcoxmod <- coxph(exprfit, data= surgerydat)
    bk <- bk_cusum(data = tdat, theta = log(2), coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE)
    expect_no_error(plot(bk))
    expect_no_error(plot(bk, h = 3))
    bk$h <- 3
    expect_no_error(plot(bk))
    bk2sided <- bk_cusum(data = tdat, theta = log(2), coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, twosided = TRUE)
    expect_no_error(plot(bk2sided))
    expect_no_error(plot(bk2sided, h = c(-3, 4)))
    bk2sided$h <- c(-3, 4)
    expect_no_error(plot(bk2sided))
  }
)

test_that(
  "CGR-CUSUM plot", {
    require(survival)
    tdat <- subset(surgerydat, unit == 1 & entrytime < 365)
    tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
    exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
    tcoxmod <- coxph(exprfit, data= surgerydat)
    cgr <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE)
    expect_no_error(plot(cgr))
    expect_no_error(plot(cgr, h = 3))
    cgr$h <- 4
    expect_no_error(plot(cgr))
  }
)

test_that(
  "funnel plot", {
    exprfitfunnel <- as.formula("(survtime <= 100) & (censorid == 1)~ age + sex + BMI")
    glmmodfun <- glm(exprfitfunnel, data = surgerydat, family = binomial(link = "logit"))
    suppressWarnings(funnel <- funnel_plot(data = surgerydat, ctime = 3*365, glmmod = glmmodfun, followup = 100))
    expect_no_error(plot(funnel))
    expect_no_error(plot(funnel, percentage = FALSE))
  }
)
