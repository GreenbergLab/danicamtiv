library(data.table)
library(ggplot2)
library(cowplot)

data <- fread("/home/brent/lasertrapr/project_dani-pCa-7.5/summary/2024-03-29_project_dani-pCa-7.5_all-measured-events.csv")
data[, time_off := (((time_off_1 + time_off_2)/2)/20000)]

colorz <- c("#666666", "#e7298a")

ggplot(data = data)+
  stat_ecdf(aes(time_off, color = conditions), pad = F)+
  scale_color_manual(values = colorz)+
  coord_cartesian(c(0, 15))+
  xlab("Time Between Events (s)")+
  ylab("Cumulative Probability")+
  theme_cowplot()


ggsave("~/Downloads/pCa-7.5-dani-time-off.png", bg = "white")

ggplot(data = data)+
  geom_boxplot(aes(conditions, time_off, fill = conditions), show.legend = F)+
  scale_y_log10()+
  ylab("Time Between Events (s)")+
  scale_fill_manual(values = colorz)+
  theme_cowplot()



kruskal.test(time_off ~ conditions, data = data)
