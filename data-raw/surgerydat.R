# code to prepare `surgerydat` dataset goes here
#Multiple hospitals perform the same surgery (procedure) - prolong life on
#avg by 1 year
#Patients arrive according to poisson process
#Failure - patient experiences an unfavourable outcome (a while) after surgery
#Covariates - age, sex, BMI
#Failures according to exponential lambda = 0.01 (expect to fail after 100 days)
#time in days


set.seed(01041996)
surgerydat <- data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, c("entrytime", "survtime", "censorid", "hosp_num", "expmu", "psival", "age", "sex", "BMI"))))
#censor some observations?
#i = hospital number
#j = value of psi
#k = exp(mu) >1 or not, sample from normal with high

#So we have 15 hospitals with exp(mu = 1), 15 hospitals with exp(mu = 2)
#They are divided into 3 groups of 5 hospitals with psi = 0.5, 1, 1.5
psivals <- c(0.5, 1, 1.5)
expmu <- c(1, 2)
s <- 1
for(k in 1:2){
  for(j in 1:3){
    for(i in 1:5){
        arrivtimes <- round(gen_arriv_times(psi = psivals[j], t = 400))
        n <- length(arrivtimes)
        age <- round(pmin(110, pmax(10, rnorm(n, mean = 60, 20))))
        BMI <- round(pmin(50, pmax(10, rnorm(n, mean = 25, 5))), 2)
        sex <- as.factor(sample(c("male", "female"), size = n, replace = TRUE, prob = c(0.45, 0.55)))
        censorid <- sample(c(0, 1), size = n, replace = TRUE, prob = c(0.02, 0.98))
        tdat <- data.frame(entrytime = arrivtimes, age = age, BMI = BMI, sex = sex)
        coxmodt <- list(formula = formula("~ age + BMI + sex"), coefficients = c(age = 0.003, BMI = 0.02, sexmale = 0.2))
        survtime <- round(gen_surv_times(invchaz = function(t) inv_chaz_exp(t, lambda = 0.01), mu = log(expmu[k]), data = tdat,
                                    coxphmod = coxmodt))
        tdat2 <- data.frame(entrytime = arrivtimes, survtime = survtime, censorid = censorid, hosp_num = rep(s, n), expmu = rep(expmu[k], n), psival = rep(psivals[j], n), age = age, sex = sex, BMI = BMI)
        surgerydat <- rbind(surgerydat, tdat2)
        s <- s+1
    }
  }
}

usethis::use_data(surgerydat, overwrite = TRUE)
