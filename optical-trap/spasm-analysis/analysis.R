library(ggplot2)
library(ggtext)
library(cowplot)
library(tidyr)
library(minpack.lm)
library(data.table)
library(magick)
devtools::load_all("~/git/spasmtools")

colorz <- c("#666666", "#e7298a")


spasm_files <- list("0 &micro;M" = "optical-trap/spasm-analysis/control/spasm-included.xlsx",
                    "10 &micro;M" = "optical-trap/spasm-analysis/10uM-dani/spasm-included-manual.xlsx"
                    )

spasm_data <- spasm_read(spasm_files)

## spasm_data <- spasm_data |> separate(id, c("isoform", "construct"), sep = "-", remove = FALSE)

con <- spasm_ensemble_average(spasm_file = "optical-trap/spasm-analysis/control/combinedEnsembleAxesData.xlsx",
                              forward_time_filter = 0.3,
                              reverse_time_filter = 0.5,
                              step_estimation_method = "fit",
                              color = colorz[1],
                              textsize = 8,
                              ylim = c(-2, 9),
                              x_shift = 0.9,
                              title = "0 &micro;M"
                              )


dani <- spasm_ensemble_average(spasm_file = "optical-trap/spasm-analysis/10uM-dani/manual-combinedEnsembleAxesData.xlsx",
                              forward_time_filter = 0.3,
                              reverse_time_filter = 0.5,
                              step_estimation_method = "fit",
                              color = colorz[2],
                              textsize = 8,
                              ylim = c(-2, 9),
                              x_shift = 0.9,
                              title = "10 &micro;M"
                              )



con_ea_f <- con$ea_f
con_ea_f$id <- "0 &micro;M"

dani_ea_f <- dani$ea_f
dani_ea_f$id <- "10 &micro;M"

ea_dat <- rbind(con_ea_f, dani_ea_f)

con_f_pred <- con$f_pred
dani_f_pred <- dani$f_pred

ggea <-
ggplot(data = dplyr::filter(ea_dat, seconds > -0.02))+
    geom_line(
              aes(x = seconds,
                  y = nanometers,
                  color = id),
              ## alpha = 0.8,
              linewidth = 0.6,
              show.legend = FALSE)+
    annotate("text",
             x = 0.19,
             y = min(ea_dat$nanometers),
             label = "0.15 s",
             vjust = 1.25,
             hjust = 0.5,
             color = "black",
             size = 12/.pt)+
    geom_line(data = con_f_pred,
              aes(x = seconds,
                  y = y_fit),
              linewidth = 0.4)+
    geom_line(data = dani_f_pred,
              aes(x = seconds,
                  y = y_fit),
              linewidth = 0.4)+
  draw_line(c(0.1, 0.26), c(-0.5))+
  draw_line(c(0.2, 0.2), c(4.5, 3), arrow = arrow(length=unit(0.20,"cm")))+
    annotate("text",
             x = 0.25,
             y = 3.5,
             label = "50%",
             vjust = 0,
             hjust = 0.5,
             color = "black",
             size = 12/.pt)+
  coord_cartesian(ylim = c(-1, 5.75))+
  scale_y_continuous(breaks = 0:5)+
  scale_color_manual(values = colorz)+
  ylab("Displacement (nm)")+
  theme_cowplot(11)+
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())



ton <- spasm_fit_ton(spasm_data = spasm_data, tmin = 0.01, colorz = colorz)

ggton <-
ton$plot+
coord_cartesian(xlim = c(0, 1))+
ggtitle("")+
ylab("Probability")+
xlab("Attachment Duration (s)")+
theme_cowplot(11)+
  theme(legend.position = "none")



## con_to_write <-dplyr::filter(ton$boot_data, id == "0 &micro;M")$boot_df[[1]]

## fwrite(con_to_write, file = "optical-trap/spasm-analysis/control/ton-boot-control.csv", col.names = F)


## dani_to_write <-dplyr::filter(ton$boot_data, id == "10 &micro;M")$boot_df[[1]]
## fwrite(dani_to_write, file = "optical-trap/spasm-analysis/10uM-dani/ton-boot-dani.csv", col.names = F)


con_trace <- readRDS("~/lasertrapr/project_dani-single-molecule/summary/figures/0uM-dani_10uM-atp_2024-04-05_obs-07_50.157-52.6062.rds")+
  coord_cartesian(xlim = c(-0.02, 2.25))+
  draw_line(x = c(0, 0.2), y = -70)+
  annotate("text", x = 0.1, y = -75, label = "0.2 s", size = 3)+
  draw_line(x = -0.05, y = c(-30, -10))+
  annotate("text", x = -0.1, y = -20, label = "20 nm", angle = 90, size = 3)+
  ggtitle("0 &micro;M danicamtiv")+
  theme_void(11)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold"))





con_trace_data <- con_trace$data

con_trace <-
  ggplot(con_trace_data)+
  geom_point(aes(new_time_index/20000, y = processed_bead_1, color = population),
             size = 0.25,
             alpha = 0.5,
             shape = 16)+
  coord_cartesian(xlim = c(-0.02, 2.25))+
  draw_line(x = c(0, 0.2), y = -70)+
  annotate("text", x = 0.1, y = -75, label = "0.2 s", size = 3)+
  draw_line(x = -0.05, y = c(-30, -10))+
  annotate("text", x = -0.1, y = -20, label = "20 nm", angle = 90, size = 3)+
  scale_color_manual(values = c("black", rep("#666666", 100)))+
  ggtitle("0 &micro;M danicamtiv")+
  theme_void(10)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold"))




dani_trace <- readRDS("~/lasertrapr/project_dani-single-molecule/summary/figures/10uM-dani_10uM-atp_2023-11-06_obs-14_28.2474-30.0897.rds")+
  coord_cartesian(xlim = c(-0.04, 1.8))+
  draw_line(x = c(0, 0.2), y = -80)+
  annotate("text", x = 0.1, y = -85, label = "0.2 s", size = 3)+
  draw_line(x = -0.05, y = c(-37, -17))+
  annotate("text", x = -0.1, y = -27, label = "20 nm", angle = 90, size = 3)+
  ggtitle("10 &micro;M danicamtiv")+
  theme_void(11)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold"))
dani_trace

dani_trace_data <- dani_trace$data

dani_trace <-
  ggplot(dani_trace_data)+
  geom_point(aes(new_time_index/20000, y = processed_bead_1, color = population),
             size = 0.25,
             alpha = 0.5,
             shape = 16)+
  coord_cartesian(xlim = c(-0.04, 1.8))+
  draw_line(x = c(0, 0.2), y = -80)+
  annotate("text", x = 0.1, y = -85, label = "0.2 s", size = 3)+
  draw_line(x = -0.05, y = c(-37, -17))+
  annotate("text", x = -0.1, y = -27, label = "20 nm", angle = 90, size = 3)+
  ggtitle("10 &micro;M danicamtiv")+
  scale_color_manual(values = c("black", rep("#e7298a", 100)))+
  theme_void(10)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold"))

fig_top <- plot_grid(con_trace, dani_trace, labels = c("A", "B"))


spasm_summary <-
  spasm_data |>
  dplyr::group_by(id)|>
  dplyr::summarize("Avg Total Step" = mean(Displacements),
                   "SD Total Step" = sd(Displacements, na.rm = T),
                   "Avg Substep 1" = mean(`Substep 1`),
                   "SD Substep 1" = sd(`Substep 1`, na.rm = T),
                   "Avg Substep 2" = mean(`Substep 2`),
                   "SD Substep 2" = sd(`Substep 2`, na.rm = T),
                   n = n()
                   )

rnorm_list <- vector("list")
  for(i in 1:nrow(spasm_summary)){
    rnorm_list[[i]] <-
      data.table(total = rnorm(10000,
                                       mean = spasm_summary$`Avg Total Step`[[i]],
                                       sd = spasm_summary$`SD Total Step`[[i]]),
                 id = spasm_summary$id[[i]]
                 )
  }
rnorm_df <- do.call("rbind", rnorm_list)

ggstep1 <-
ggplot(spasm_data)+
  stat_ecdf(aes(x = Displacements, color = id), linewidth = 1, alpha = 0.7)+
  stat_ecdf(data = rnorm_df, aes(x = total, color = id), linewidth = 0.5, linetype = "dashed")+
  scale_color_manual(values = colorz)+
  coord_cartesian(xlim = c(-15, 20))+
  ylab("Probability")+
  xlab("Total Displacement (nm)")+
  ggtitle("")+
  theme_cowplot(11)+
  theme(
  legend.position = "none"
  )

ggstep2 <- ggstep1+coord_cartesian(xlim = c(-10, 15))+theme_cowplot(11)+theme(legend.position = "none")

ggstep3 <- ggdraw(ggstep2)+draw_plot(ggstep1, 0.6, 0.15, 0.4, 0.4)

## d0_facet <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = Displacements, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   facet_wrap(~isoform)+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Total Step")+
##   theme_cowplot()

## d1 <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = `Substep 1`, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Substep 1")+
##   theme_cowplot()+
##   theme(
##   legend.position = "none"
##   )



## d1_facet <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = `Substep 1`, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   facet_wrap(~isoform)+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Substep 1")+
##   theme_cowplot()

## d2 <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = `Substep 2`, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Substep 2")+
##   theme_cowplot()

## d2_facet <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = `Substep 2`, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   facet_wrap(~isoform)+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Substep 2")+
##   theme_cowplot()

bottom <- plot_grid(ggea, ggstep1, ggton, nrow = 1, labels = c("C", "D", "E"))
gg_trap <- image_read_svg("figures/trap-cartoon.svg", height = 10000, width = 10000)
gg_trap <- ggplot()+draw_image(gg_trap)+theme_nothing()
traces_right <- plot_grid(con_trace, dani_trace, nrow = 2)
top <- plot_grid(gg_trap, traces_right, rel_widths = c(1, 1.2), labels = c("A", "B"))

png("figures/for-the-boss2/figure-3.png", width = 6.5, height = 4.5, units = "in", res = 1000)
plot_grid(top, bottom, nrow = 2, rel_heights = c(1.25, 1))
dev.off()


pdf("figures/for-the-boss2/figure-3.pdf", width = 6.5, height = 4.5)
plot_grid(top, bottom, nrow = 2, rel_heights = c(1.25, 1))
dev.off()






ton_df <- dplyr::select(ton$boot_data, id, html_label)

spasm_summary <-
  spasm_data |>
  dplyr::group_by(id)|>
  dplyr::summarize("Avg Total Step" = mean(Displacements),
                   "SD Total Step" = sd(Displacements, na.rm = T),
                   "Avg Substep 1" = mean(`Substep 1`),
                   "SD Substep 1" = sd(`Substep 1`, na.rm = T),
                   "Avg Substep 2" = mean(`Substep 2`),
                   "SD Substep 2" = sd(`Substep 2`, na.rm = T),
                   n = n()
                   )


dat <- left_join(spasm_summary, ton_df)


write.csv(dat, "summary-table.csv")






con$mod_f
