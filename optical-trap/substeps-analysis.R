library(data.table)
library(ggplot2)
library(viridis)
library(cowplot)

substeps_files <- list.files("~/lasertrapr/project_dani-single-molecule",
                             pattern = "substeps.csv",
                             full.names = TRUE,
                             recursive = TRUE)

substep_data <- rbindlist(lapply(substeps_files, fread))
substep_data$conditions <- factor(substep_data$conditions, levels = c("0uM-dani_10uM-atp", "1uM-dani_10uM-atp", "10uM-dani_10uM-atp", "100uM-dani_10uM-atp"))
substep_data$conditions2 <- sub("_10uM-atp","", substep_data$conditions)
substep_data$conditions2 <- sub("10uM-atp_","", substep_data$conditions2)

substep_data$conditions2 <- factor(substep_data$conditions2, levels = c("0uM-dani",
                                                     "1uM-dani",
                                                     "10uM-dani",
                                                     "100uM-dani"))

substep_data <- substep_data[, .(prior_unbound_position_nm = mean(prior_unbound_position_nm),
                                 bead_position_substep_1_nm = mean(bead_position_substep_1_nm),
                                 substep_1_nm = mean(substep_1_nm),
                                 bead_position_substep_2_nm = mean(bead_position_substep_2_nm),
                                 substep_2_nm = mean(substep_2_nm),
                                 after_unbound_position_nm = mean(after_unbound_position_nm)),
                             by = c("project", "conditions", "conditions2", "date", "obs", "event_id")
                             ]

substep_data[, total_step := substep_1_nm + substep_2_nm ]
substep_data[, total_step2 := bead_position_substep_2_nm - after_unbound_position_nm ]
substep_data[, substep_2_nm2 := total_step2 - substep_1_nm ]


ggplot(substep_data)+
  geom_histogram(aes(x = total_step,
                     y = after_stat(density),
                     fill = conditions),
                 color = "black",
                 binwidth = 2,
                 size = 0.2)+
  facet_wrap(~conditions)

ggplot(substep_data)+
  stat_ecdf(aes(x = total_step, color = conditions), linewidth = 1)+
  scale_color_manual(values = c("#666666", RColorBrewer::brewer.pal(8, "Dark2")[c(3, 4, 5, 6)]))+
  facet_wrap(~conditions)

shapiro.test(substep_data[conditions=="0uM-dani_10uM-atp"]$total_step)
shapiro.test(substep_data[conditions=="1uM-dani_10uM-atp"]$total_step)
shapiro.test(substep_data[conditions=="10uM-dani_10uM-atp"]$total_step)
shapiro.test(substep_data[conditions=="100uM-dani_10uM-atp"]$total_step)

ggplot(substep_data)+
  geom_histogram(aes(x = substep_1_nm,
                     y = after_stat(density),
                     fill = conditions),
                 color = "black",
                 binwidth = 2,
                 size = 0.2)+
  facet_wrap(~conditions, ncol = 1)

ggplot(substep_data)+
  geom_histogram(aes(x = substep_2_nm,
                     y = after_stat(density),
                     fill = conditions),
                 color = "black",
                 binwidth = 1,
                 size = 0.2)+
  facet_wrap(~conditions, ncol = 1)

shapiro.test(substep_data[conditions=="0uM-dani_10uM-atp"]$substep_1_nm)
shapiro.test(substep_data[conditions=="1uM-dani_10uM-atp"]$substep_1_nm)
shapiro.test(substep_data[conditions=="10uM-dani_10uM-atp"]$substep_1_nm)
shapiro.test(substep_data[conditions=="100uM-dani_10uM-atp"]$substep_1_nm)


shapiro.test(substep_data[conditions=="0uM-dani_10uM-atp"]$substep_2_nm)
shapiro.test(substep_data[conditions=="1uM-dani_10uM-atp"]$substep_2_nm)
shapiro.test(substep_data[conditions=="10uM-dani_10uM-atp"]$substep_2_nm)
shapiro.test(substep_data[conditions=="100uM-dani_10uM-atp"]$substep_2_nm)


ggplot(substep_data, aes(x = substep_1_nm, y = substep_2_nm))+
  ## geom_tile(aes(x = substep_1_nm, y = substep_2_nm, fill = after_stat(density)))+
  stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  scale_fill_viridis(discrete  = FALSE, option = "mako", direction = -1)+
  facet_wrap(~conditions2, ncol = 1)+
  xlab("Substep 1 (nm)")+
  ylab("Substep 2 (nm)")+
  coord_cartesian(xlim = c(-20, 20), y = c(-20, 20))+
  theme_cowplot()+
  theme(
   strip.background = element_rect(fill = "white")
  )





substep_data[, .(s1_avg = mean(substep_1_nm),
                s2_avg = mean(substep_2_nm),
                total_step_avg = mean(total_step)),
            by = conditions2]
