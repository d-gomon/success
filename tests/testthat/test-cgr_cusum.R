require(survival)
#preliminaries
tdat <- subset(surgerydat, unit == 1 & entrytime < 100)
tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)
exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
tcoxmod <- coxph(exprfit, data= surgerydat)




# Test if no active individuals is handled appropriately
# Should be a straight line between time 5 and time 10
tdat <- data.frame(entrytime = c(1,  10), survtime = c(4, 5), censorid = rep(1, 2))
cgr_inactive <- cgr_cusum(tdat, cbaseh = tcbaseh, ctimes = seq(0, 15, 1))

test_that("Straight line when no individuals",{
  expect_equal(cgr_inactive$CGR[7, 2], cgr_inactive$CGR[12, 2])
})


#Known issues:

#When patients with survtime = 0 are present, the ML estimate cannot be calculated
#at the failure time of the patient.
#If no patients are present in dataset at that point in time, then at the time
#of entry of the next patient into the data set S_2, we can calculate CGR(t) at
#S_2 + \epsilon to obtain an arbitrarily large value of CGR(t). This is because
#thetat = log(NDT/AT) = log(1/epsilon) -> infty when epsilon -> 0.
#Demonstration:
#cgr <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, ctimes = seq(5, 30, 1))
#cgr2 <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, ctimes = seq(5, 30, 0.1))
#cgr3 <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, ctimes = seq(5, 30, 0.01))
#plot the CGR-CUSUM and note how the value can increase infinitely
#plot(cgr)
#plot(cgr2)
#plot(cgr3)


cgr1c <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh)
#cgr3c <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, ncores = 3, dependencies = list("chaz_exp"))

maxtt <- log(6)
cgr1clower <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, detection = "lower", maxtheta = maxtt)

#test_that("cores don't influence results",{
#          skip_on_cran()
#          expect_equal(cgr1c$CGR, cgr3c$CGR)})



test_that("maxtheta is working as intended",{
  expect_equal(all(cgr1c$CGR$exp_theta_t <= exp(maxtt)), TRUE)
  expect_equal(all(cgr1clower$CGR$exp_theta_t >= exp(-maxtt)), TRUE)
})



# Profvis time checks
tdat <- subset(surgerydat, unit <= 15 & entrytime < 100)
exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
tcoxmod <- coxph(exprfit, data= surgerydat)
cppcus <- cgr_cusum(data = tdat, coxphmod = tcoxmod)
normcus <- cgr_cusum(data = tdat, coxphmod = tcoxmod)
library(profvis)

profvis({
  cgr_cusum(data = tdat, coxphmod = tcoxmod)
})






