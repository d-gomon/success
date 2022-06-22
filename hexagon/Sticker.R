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
ic_hosp <- generate_units(time = time, psi = psi, n_sim = n_sim, cbaseh = cbaseh,
                          inv_cbaseh = inv_cbaseh)

oc_hosp_up <- generate_units(time = time, psi = psi, n_sim = n_sim, cbaseh = cbaseh,
                          inv_cbaseh = inv_cbaseh, mu = log(3))

oc_hosp_down <- generate_units(time = time, psi = psi, n_sim = n_sim, cbaseh = cbaseh,
                             inv_cbaseh = inv_cbaseh, mu = -log(3))

bk_oc_up <- bk_cusum(data = subset(oc_hosp_up, unit == 1), cbaseh = function(t) chaz_exp(t, lambda = 0.02),theta = log(2), stoptime = 50)
bk_oc_down <- bk_cusum(data = subset(oc_hosp_down, unit == 2), cbaseh = function(t) chaz_exp(t, lambda = 0.02),theta = -log(2), stoptime = 50)

#cgr_oc_up <- cgr_cusum(data = subset(oc_hosp_up, unit == 1), cbaseh = function(t) chaz_exp(t, lambda = 0.02), stoptime = 50)
#cgr_oc_down <- cgr_cusum(data = subset(oc_hosp_down, unit == 2), cbaseh = function(t) chaz_exp(t, lambda = 0.02), stoptime = 50, detection ="lower")

c <- ggplot()
for(i in 1:n_sim){
  bktemp <- bk_cusum(data = subset(ic_hosp, unit == i), cbaseh = function(t) chaz_exp(t, lambda = 0.02),theta = log(2), stoptime = 50)
  bktemp2 <- bk_cusum(data = subset(ic_hosp, unit == i), cbaseh = function(t) chaz_exp(t, lambda = 0.02),theta = -log(2), stoptime = 50)
  c <- c +
    geom_line(data = bktemp$BK, mapping= aes(x = time, y = value),
                     col = "grey", size = 0.6, lty = 1) +
    geom_line(data = bktemp2$BK, mapping= aes(x = time, y = value),
              col = "grey", size = 0.6, lty = 1)
}
c <- c +
  geom_line(data = bk_oc_up$BK, mapping= aes(x = time, y = value),
            col = "red", size = 1, lty = 1) +
  geom_line(data = bk_oc_down$BK, mapping= aes(x = time, y = value),
            col = "green", size = 1, lty = 1) +
  geom_hline(yintercept = 12, col = "red", lty = 2, size = 1.3) +
  geom_hline(yintercept = -12, col = "green", lty = 2, size = 1.3) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),  #remove y axis ticks
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
  )+ scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  geom_hline(yintercept = 0, size = 1)


ggsave(filename = "success.png", plot = c,
       path = dirname(rstudioapi::getActiveDocumentContext()$path),
       height = 700, width = 1000, units = "px")



rm(list = ls())

library(showtext)
## Loading Google fonts (fonts.google.com)
#font_add_google("Source Code Pro", "myfont")
font_add_google("PT mono", "myfont")

library(hexSticker)

#col = 'aquamarine2'
#col = 'greenyellow'
#col = 'deepskyblue4'
col = 'dodgerblue4'
asd <- getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sticker("success.png", dpi = 1000,
        s_x=1, s_y=1.08, s_width = 1, s_height = 1,
        package = "success", p_size = 40, p_color = col,
        p_family = 'myfont',
        p_x = 1, p_y = 0.30,
        h_fill = 'white', h_color = "red",
        asp = 0.9, # default asp = 1 sformatta input figure!!!
        filename="success_logo.png",
        white_around_sticker = TRUE)
setwd(asd)


