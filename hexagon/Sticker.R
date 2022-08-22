library(ggplot2)
ic_chart <- data.frame(time = c(0, 2.5, 4, 6),
                       value = c(0, 2.5/(cot((30/360)*2*pi)), 0.7, 6))
oc_chart <- data.frame(time = c(0, 2, 3, 6), value = c(0, -2, -0.5, -4))

#Tuning parameters
t_size <- 1.2


#initiate plot
c <- ggplot()

#Add chart lines
c <- c +
  geom_line(data = ic_chart, mapping= aes(x = time, y = value),
            col = "red", size = t_size, lty = 1)+
  geom_line(data = oc_chart, mapping= aes(x = time, y = value),
            col = "green", size = t_size, lty = 1)

#Add control limits
c <- c + geom_hline(yintercept = 7, lty = 5, col = "red", size = t_size) +
  geom_hline(yintercept = -6, col = "green", lty = 5, size = t_size)

#Add arrows
c <- c + geom_segment(aes(x = ic_chart$time[3], y = ic_chart$value[3], xend = ic_chart$time[4], yend = ic_chart$value[4]),
                      arrow = arrow(length = unit(0.5, "cm")), col = "red", size = t_size)
c <- c + geom_segment(aes(x = oc_chart$time[3], y = oc_chart$value[3], xend = oc_chart$time[4], yend = oc_chart$value[4]),
                      arrow = arrow(length = unit(0.5, "cm")), col = "green", size = t_size)


#Remove stuff from theme
c <- c +
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
  scale_y_continuous(expand=c(0,0))


ggsave(filename = "success2.png", plot = c,
       path = dirname(rstudioapi::getActiveDocumentContext()$path),
       height = 2000, width = 1728, units = "px")



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
sticker("success_gimp.png", dpi = 1000,
        s_x=1, s_y=1, s_width = 1, s_height = 1,
        package = "", p_size = 1, p_color = "red",
        p_family = 'myfont',
        p_x = 1, p_y = 0.30,
        h_fill = 'white', h_color = "red",
        asp = 0.9, # default asp = 1 sformatta input figure!!!
        filename="success_logo2.png",
        white_around_sticker = FALSE)
setwd(asd)

asd <- getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sticker("success_gimp_v3.png", dpi = 1000,
        s_x=1, s_y=1, s_width = 1, s_height = 1,
        package = "", p_size = 1, p_color = "red",
        p_family = 'myfont',
        p_x = 1, p_y = 0.30,
        h_fill = 'white', h_color = "red",
        asp = 0.9, # default asp = 1 sformatta input figure!!!
        filename="success_logo3.png",
        white_around_sticker = FALSE)
setwd(asd)

asd <- getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sticker("TyronIdee.png", dpi = 1000,
        s_x=1, s_y=1, s_width = 1, s_height = 1,
        package = "", p_size = 1, p_color = "red",
        p_family = 'myfont',
        p_x = 1, p_y = 0.30,
        h_fill = 'black', h_color = "black",
        asp = 0.9, # default asp = 1 sformatta input figure!!!
        filename="success_logo_Tyron.png",
        white_around_sticker = FALSE)
setwd(asd)


