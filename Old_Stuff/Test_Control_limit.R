op <- subset(surgerydat, hosp_num == 1 & survtime > 0)

rcox <- coxph(exprfit, data = op)
po <- basehaz(rcox, centered = FALSE)

library(tictoc)

tic("INV_CBASEH")
a <- cgr_control_limit(time = 500, alpha = 0.1, cbaseh = function(t) chaz_exp(t, lambda = 0.02),
                       inv_cbaseh = function(t) inv_chaz_exp(t, lambda = 0.02), psi = 0.5, n_sim = 10)
toc()

tic("coxmod + resampling")
b <- cgr_control_limit(time = 500, alpha = 0.1, coxphmod = tcoxmod, psi = 0.5, n_sim = 10,
                       data = subset(surgerydat, hosp_num == 1))
toc()

tic("coxmod without resampling")
b2 <- cgr_control_limit(time = 500, alpha = 0.1, coxphmod = tcoxmod, psi = 0.5, n_sim = 10)
toc()

tic("cbaseh without resampling")
c <- cgr_control_limit(time = 500, alpha = 0.1, cbaseh = function(t) chaz_exp(t, lambda = 0.02), psi = 0.5, n_sim = 10)
toc()
