#We create a list of in-control and out-of-control instances, both upwards
#and downwards. Then we grey out the in-control instances, and highlight the out-of-control ones
#Create pictures which shows this and make sticker

load_all()
set.seed(01041996)
time <- 100
psi <- 2
n_sim <- 20
cbaseh <- function(t) chaz_exp(t, lambda = 0.02)
inv_cbaseh <- function(t) inv_chaz_exp(t, lambda = 0.02, mu = log(1))
coxphmod <- NULL
baseline_data <- NULL
interval <- c(0, 9e12)
df_temp <- generate_units(time = time, psi = psi, n_sim = n_sim, cbaseh = cbaseh,
                          inv_cbaseh = inv_cbaseh)

df_temp2 <- generate_units(time = time, psi = psi, n_sim = n_sim, cbaseh = cbaseh,
                          inv_cbaseh = inv_cbaseh, mu = log(1.6))

plot(bk_cusum(data = subset(df_temp, unit == 1), cbaseh = function(t) chaz_exp(t, lambda = 0.02),theta = -log(2), stoptime = 100))
plot(bk_cusum(data = subset(df_temp, unit == 1), cbaseh = function(t) chaz_exp(t, lambda = 0.02),theta = log(2), stoptime = 100))


plot(bk_cusum(data = subset(df_temp2, unit == 1), cbaseh = function(t) chaz_exp(t, lambda = 0.02),theta = -log(2), stoptime = 100))
plot(bk_cusum(data = subset(df_temp2, unit == 1), cbaseh = function(t) chaz_exp(t, lambda = 0.02),theta = log(2), stoptime = 100))
