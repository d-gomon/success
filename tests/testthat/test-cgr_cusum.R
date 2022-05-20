#Known issues:

#When patients with survtime = 0 are present, the ML estimate cannot be calculated
#at the failure time of the patient.
#If no patients are present in dataset at that point in time, then at the time
#of entry of the next patient into the data set S_2, we can calculate CGR(t) at
#S_2 + \epsilon to obtain an arbitrarily large value of CGR(t). This is because
#thetat = log(NDT/AT) = log(1/epsilon) -> infty when epsilon -> 0.
#Demonstration:
require(survival)
#Select only the data of the first year of the first hospital in the surgerydat data set
tdat <- subset(surgerydat, unit == 1 & entrytime < 365)

#We know that the cumulative baseline hazard in the data set is
#Exponential(0.01). If you don't know the cumulative baseline, we suggest
#leaving the cbaseh argument empty and determining a coxphmod (see help)
tcbaseh <- function(t) chaz_exp(t, lambda = 0.01)

#Determine a risk-adjustment model using a Cox proportional hazards model.
#Outcome (survival) is regressed on the available covariates:
exprfit <- as.formula("Surv(survtime, censorid) ~ age + sex + BMI")
tcoxmod <- coxph(exprfit, data= surgerydat)

#Determine the values of the chart
cgr <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, ctimes = seq(5, 30, 1))
cgr2 <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, ctimes = seq(5, 30, 0.1))
cgr3 <- cgr_cusum(data = tdat, coxphmod = tcoxmod, cbaseh = tcbaseh, pb = TRUE, ctimes = seq(5, 30, 0.01))
#plot the CGR-CUSUM (exact hazard)
plot(cgr)
plot(cgr2)
plot(cgr3)

TRUE
