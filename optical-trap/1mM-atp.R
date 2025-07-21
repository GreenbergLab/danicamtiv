library(data.table)
library(ggplot2)
library(cowplot)

dat <- fread("optical-trap/lasertrapr/project_dani-1mM-atp/summary/2024-06-03_project_dani-1mM-atp_all-measured-events.csv")

colorz <- c("#666666", "#e7298a")

dat$conditions <- factor(dat$conditions, levels = c("control", "10uM-dani_1mM-atp"))

ggstep <-
ggplot(dat)+
  stat_ecdf(aes(displacement_nm, color = conditions), pad = F, linewidth = 0.8)+
  scale_color_manual(values = colorz)+
  ## coord_cartesian(xlim = c(-15, 20))+
  ylab("Probability")+
  xlab("Total Displacement (nm)")+
  ggtitle("")+
  theme_cowplot(11)+
  theme(
  legend.position = "none"
  )


ggton <-
ggplot(dat)+
  stat_ecdf(aes(time_on_ms/1000, color = conditions), pad = F, linewidth = 0.8)+
  scale_color_manual(values = colorz)+
  coord_cartesian(xlim = c(0, 0.25))+
  ylab("Probability")+
  xlab("Attachment Duration (s)")+
  ggtitle("")+
  theme_cowplot(11)+
  theme(
  legend.position = "none"
  )

png("supplemental-figures/supp-fig-4-1mM-atp-trap.png", width = 6, height = 2.5, units = "in", res = 300)
plot_grid(ggstep, ggton, labels = LETTERS)
dev.off()


pdf("supplemental-figures/supp-fig-4-1mM-atp-trap.pdf", width = 6, height = 2.5)
plot_grid(ggstep, ggton, labels = LETTERS)
dev.off()


dat[, .(mean = mean(displacement_nm),
        sd = sd(displacement_nm)),
    by = conditions]
