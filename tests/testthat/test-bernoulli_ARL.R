# ####################################Tests############################

test_that("Inputs", {
  expect_error(bernoulli_ARL(h = 2, t = 0, p0 = 0.5, theta = log(2)),
               "Parameter 't' must be an integer larger than 1.")
  expect_error(bernoulli_ARL(h = 2, t = 3.3, p0 = 0.5, theta = log(2)),
               "Parameter 't' must be an integer larger than 1.")
  expect_error(bernoulli_ARL(h = c(2,3), t = 4, p0 = 0.5, theta = log(2)),
               "Parameter 'h' must be a single numeric variable.")
  expect_error(bernoulli_ARL(h = 2, t = 4, p0 = 0.5, theta = log(2), smooth_prob = 3),
               "Parameter 'smooth_prob' must be a single logical value.")
  expect_error(bernoulli_ARL(h = 2, t = 4, p0 = 0.5, theta = log(2), followup = "death"),
               "Parameter 'followup' must be a positive numeric variable.")
  expect_error(bernoulli_ARL(h = 2, t = 4, p0 = 0.5, theta = log(2), glmmod = 0.3),
               "Parameter 'glmmod' must be of class 'glm'.")
  expect_error(bernoulli_ARL(h = 2, t = 4, p0 = 0.5, theta = 0),
               "Parameter 'theta' must be a numeric variable not equal to 0.")
  expect_error(bernoulli_ARL(h = 2, t = 4, p0 = 1.2, p1 = 0.8),
               "Parameter 'p0' must be a positive probability between 0 and 1. (numeric)", fixed = TRUE)
  expect_error(bernoulli_ARL(h = 2, t = 4, p0 = 0.5, p1 = 1.2),
               "Parameter 'p1' must be a positive probability between 0 and 1. (numeric)", fixed = TRUE)
})




test_that("Theory", {
  glmmodber <- glm((survtime <= 100) & (censorid == 1)~ age + sex + BMI, data = surgerydat, family = binomial(link = "logit"))
  ARL <- bernoulli_ARL(h = 2.5, t = 100, glmmod = glmmodber, theta = log(2))
  t <- 100
  #Check whether theta is well defined
  expect_equal(bernoulli_ARL(h = 2, t = t, p0 = 0.5, theta = log(2)), bernoulli_ARL(h = 2, t = t, p0 = 0.5, p1 = 2/3))
  #Create inverted glmmod
  glmmod2 <- glmmodber
  glmmod2$fitted.values <- 1 - glmmod2$fitted.values
  #Replacing theta with -theta we should get the same answer when the specified probabilities are inverted (1-p).
  expect_equal(ARL, bernoulli_ARL(h = 2.5, t = 100, glmmod = glmmod2, theta = -log(2)))
  expect_equal(bernoulli_ARL(h = 2.5, t = 100, glmmod = glmmodber, theta = -log(2)), bernoulli_ARL(h = 2.5, t = 100, glmmod = glmmod2, theta = log(2)))
  #Same for probabilities specified:
  expect_equal(bernoulli_ARL(h = 2.5, t = 100, p0 = 0.3, theta = log(2)), bernoulli_ARL(h = 2.5, t = 100, p0 = 0.7, theta = -log(2)))
  #Equal when theta_true = 0
  expect_equal(ARL, bernoulli_ARL(h = 2.5, t = 100, glmmod = glmmodber, theta = log(2), theta_true = log(1)))
  #ARL_0 Smaller when theta_true > 0
  expect_gt(ARL$ARL_0, bernoulli_ARL(h = 2.5, t = 100, glmmod = glmmodber, theta = log(2), theta_true = log(1.5))$ARL_0)
})


test_that("Outputs",{
  glmmodber <- glm((survtime <= 100) & (censorid == 1)~ age + sex + BMI, data = surgerydat, family = binomial(link = "logit"))
  ARL <- bernoulli_ARL(h = 2.5, t = 100, glmmod = glmmodber, theta = log(2))
  t <- 100
  expect_true(all(dim(ARL$R) == t))
})






# ####################################Profiling Code#############################
# library(profvis)
# #We consider patient outcomes 100 days after their entry into the study.
# followup <- 100
# #Determine a risk-adjustment model using a generalized linear model.
# #Outcome (failure within 100 days) is regressed on the available covariates:
# exprfitber <- as.formula("(survtime <= followup) & (censorid == 1)~ age + sex + BMI")
# glmmodber <- glm(exprfitber, data = surgerydat, family = binomial(link = "logit"))
# profvis({bernoulli_ARL(h = 2.5, t = 200, glmmod = glmmodber, theta = log(2))})
# profvis({bernoulli_ARL(h = 2.5, t = 200, p0 = 0.5, theta = log(2))})
#
# ####################################Speed testing#############################
# library(tictoc)
# tic("t = 100")
# ARL <- bernoulli_ARL(h = 2.5, t = 100, glmmod = glmmodber, theta = log(2))
# toc()
# tic("t = 300")
# ARL <- bernoulli_ARL(h = 2.5, t = 300, glmmod = glmmodber, theta = log(2))
# toc()
# tic("t = 1000")
# ARL <- bernoulli_ARL(h = 2.5, t = 1000, glmmod = glmmodber, theta = log(2))
# toc()