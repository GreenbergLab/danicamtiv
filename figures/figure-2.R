library(ggplot2)
library(ggtext)
library(cowplot)
library(tidyr)
library(minpack.lm)
library(data.table)
library(magick)
library(ggbeeswarm)

devtools::load_all("~/git/spasmtools")

theme_set(theme_cowplot(10))

colorz <- c("#666666",
            "#f7b8d8",
            "#ef72b2",
            "#e7298a")


spasm_files <- list("0 &micro;M" = "optical-trap/spasm-analysis/control/spasm-included.xlsx",
                    "1 &micro;M" = "optical-trap/spasm-analysis/1uM-dani/spasm_1uM-dani_10uM-atp.xlsx",
                    "10 &micro;M" = "optical-trap/spasm-analysis/10uM-dani/spasm-included-manual.xlsx",
                    "100 &micro;M" = "optical-trap/spasm-analysis/100uM-dani/spasm_100uM-dani_10uM-atp.xlsx"
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
                              color = colorz[3],
                              textsize = 8,
                              ylim = c(-2, 9),
                              x_shift = 0.9,
                              title = "10 &micro;M"
                              )



dani_1uM <- spasm_ensemble_average(spasm_file = "optical-trap/spasm-analysis/1uM-dani/1uM-dani_10uM-atp_combinedEnsembleAxesData.xlsx",
                              forward_time_filter = 0.3,
                              reverse_time_filter = 0.5,
                              step_estimation_method = "fit",
                              color = colorz[2],
                              textsize = 8,
                              ylim = c(-2, 9),
                              x_shift = 0.9,
                              title = "1 &micro;M"
                              )


dani_100uM <- spasm_ensemble_average(spasm_file = "optical-trap/spasm-analysis/100uM-dani/100uM-dani_10uM-atp_combinedEnsembleAxesData.xlsx",
                              forward_time_filter = 0.3,
                              reverse_time_filter = 0.5,
                              step_estimation_method = "fit",
                              color = colorz[4],
                              textsize = 8,
                              ylim = c(-2, 9),
                              x_shift = 0.9,
                              title = "100 &micro;M"
                              )


con_f_final_nm <- mean(tail(con$ea_f$nanometers), 200)
con_f_final_nm
mean(tail(dani_1uM$ea_f$nanometers), 200)/con_f_final_nm
mean(tail(dani$ea_f$nanometers), 200) / con_f_final_nm
mean(tail(dani_100uM$ea_f$nanometers), 200) /con_f_final_nm


dat_labels <- data.frame(y = c(mean(tail(con$ea_f$nanometers, 200)),
                                   mean(tail(dani_1uM$ea_f$nanometers), 200) ,
                                   mean(tail(dani$ea_f$nanometers), 200)+0.06,
                                   mean(tail(dani_100uM$ea_f$nanometers), 200)-0.085),
                         x = 0.355,
                         id = c("0 &micro;M", "1 &micro;M", "10 &micro;M", "100 &micro;M"))


con_ea_f <- con$ea_f
con_ea_f$id <- "0 &micro;M"

dani_ea_f <- dani$ea_f
dani_ea_f$id <- "10 &micro;M"


dani_1uM_ea_f <- dani_1uM$ea_f
dani_1uM_ea_f$id <- "1 &micro;M"


dani_100uM_ea_f <- dani_100uM$ea_f
dani_100uM_ea_f$id <- "100 &micro;M"

ea_dat <- rbind(con_ea_f, dani_ea_f, dani_1uM_ea_f, dani_100uM_ea_f)

con_f_pred <- con$f_pred
dani_f_pred <- dani$f_pred
dani_1uM_f_pred <- dani_1uM$f_pred
dani_100uM_f_pred <- dani_100uM$f_pred

(
ggea <-
ggplot(data = dplyr::filter(ea_dat, seconds > -0.02 & seconds < 0.25))+
    geom_line(
              aes(x = seconds,
                  y = nanometers,
                  color = id),
              ## alpha = 0.8,
              linewidth = 0.6,
              show.legend = FALSE)+
  geom_richtext(data = dat_labels,
            aes(x = x,
                y = y,
                color = id,
                label = id),
            label.color = NA,
            fill = "transparent",
            hjust = 1,
            size = 8/.pt)+
    annotate("text",
             x = 0.175,
             y = min(ea_dat$nanometers),
             label = "0.15 s",
             vjust = 1,
             hjust = 0.5,
             color = "black",
             size = 10/.pt)+
    ## geom_line(data = con_f_pred,
    ##           aes(x = seconds,
    ##               y = y_fit),
    ##           linewidth = 0.4)+
    ## geom_line(data = dani_f_pred,
    ##           aes(x = seconds,
    ##               y = y_fit),
    ##           linewidth = 0.4)+
  draw_line(c(0.1, 0.25), c(-0.85))+
  ## draw_line(c(0.2, 0.2), c(4.5, 3), arrow = arrow(length=unit(0.20,"cm")))+
  ##   annotate("text",
  ##            x = 0.25,
  ##            y = 3.5,
  ##            label = "50%",
  ##            vjust = 0,
  ##            hjust = 0.5,
  ##            color = "black",
  ##            size = 12/.pt)+
  coord_cartesian(ylim = c(-1.5, 5.75))+
  scale_y_continuous(breaks = 0:5)+
  scale_color_manual(values = c(colorz[1], colorz[2], colorz[3], colorz[4]))+
  ylab("Displacement (nm)")+
  ggtitle("Working Stroke")+
  theme_cowplot(10)+
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_markdown(size = 10),
        plot.title = element_markdown(size = 9, hjust = 0.5),
        legend.position = "none")
)


## ggsave("supplemental-figures/ensemble-average-0-1-10-100-uM-dani.png", bg = "white")




ton <- spasm_fit_ton(spasm_data = spasm_data, tmin = 0.01, colorz = c(colorz[1],
                                                                      colorz[2],
                                                                      colorz[3],
                                                                      colorz[4]))

label_ton <- paste0("<span style='color:#666666'>0 &micro;M = ", ton$boot_data$html_label[[1]],"  </span><br>",
                    "<span style='color:#f7b8d8'> 1 &micro;M = ", ton$boot_data$html_label[[2]],"</span><br>",
                    "<span style='color:#ef72b2'> 10 &micro;M = ", sub("24", "24.0", ton$boot_data$html_label[[3]]),"  </span><br>",
                    "<span style='color:#e7298a'> 100 &micro;M = ", ton$boot_data$html_label[[4]],"  </span><br>" )


(
  ggton <-
ton$plot+
  annotate("richtext",
           x = 0.6,
               y = 0.5,
               label = label_ton,
            label.color = "transparent",
           label.size = 1,
            fill = "transparent",
            hjust = 0.5,
           vjust = 0.5,
            size = 8/.pt)+
coord_cartesian(xlim = c(0, 1))+
## scale_x_continuous(expand = expansion(c(NA, 0), c(NA, 0)))+
ggtitle("")+
ylab("Probability")+
xlab("Time (s)")+
ggtitle("Detachment Rate")+
theme_cowplot(10)+
  theme(
    legend.position = "none",
    plot.title = element_markdown(size = 9, hjust = 0.5),
    ## plot.margin = margin(r = 30, l = 30, t=5, b = 5),
   axis.title =  element_markdown(size = 10)
    )
  )

## ggsave("supplemental-figures/all-concs-trap-ton.png", bg = "white")

## con_to_write <-dplyr::filter(ton$boot_data, id == "0 &micro;M")$boot_df[[1]]

## fwrite(con_to_write, file = "optical-trap/spasm-analysis/control/ton-boot-control.csv", col.names = F)


## dani_to_write <-dplyr::filter(ton$boot_data, id == "10 &micro;M")$boot_df[[1]]
## fwrite(dani_to_write, file = "optical-trap/spasm-analysis/10uM-dani/ton-boot-dani.csv", col.names = F)


con_trace <- readRDS("optical-trap/lasertrapr/project_dani-single-molecule/summary/figures/0uM-dani_10uM-atp_2024-04-05_obs-07_50.157-52.6062.rds")+
  coord_cartesian(xlim = c(-0.021, 2.25))+
  draw_line(x = c(0, 0.2), y = -70)+
  annotate("text", x = 0.1, y = -75, label = "0.2 s", size = 8/.pt)+
  draw_line(x = -0.05, y = c(-30, -10))+
  annotate("text", x = -0.1, y = -20, label = "20 nm", angle = 90, size = 8/.pt)+
  ggtitle("0 &micro;M danicamtiv")+
  theme_void(10)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold", size = 10))





con_trace_data <- con_trace$data

con_trace <-
  ggplot(con_trace_data)+
  geom_point(aes(new_time_index/20000, y = processed_bead_1, color = population),
             size = 0.25,
             alpha = 0.5,
             shape = 16)+
  coord_cartesian(xlim = c(-0.023, 2.25))+
  draw_line(x = c(0, 0.2), y = -70)+
  annotate("text", x = 0.1, y = -75, label = "0.2 s", size = 8/.pt)+
  draw_line(x = -0.05, y = c(-30, -10))+
  annotate("text", x = -0.1, y = -20, label = "20 nm", angle = 90, size = 8/.pt)+
  scale_color_manual(values = c("black", rep("#666666", 100)))+
  ggtitle("0 &micro;M danicamtiv")+
  theme_void(10)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold", size = 10))




dani_trace <- readRDS("optical-trap/lasertrapr/project_dani-single-molecule/summary/figures/10uM-dani_10uM-atp_2023-11-06_obs-14_28.2474-30.0897.rds")+
  coord_cartesian(xlim = c(-0.041, 1.8))+
  draw_line(x = c(0, 0.2), y = -80)+
  annotate("text", x = 0.1, y = -85, label = "0.2 s", size = 8/.pt)+
  draw_line(x = -0.05, y = c(-37, -17))+
  annotate("text", x = -0.1, y = -27, label = "20 nm", angle = 90, size = 8/.pt)+
  ggtitle("10 &micro;M danicamtiv")+
  theme_void(10)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold", size = 10))
dani_trace

dani_trace_data <- dani_trace$data

dani_trace <-
  ggplot(dani_trace_data)+
  geom_point(aes(new_time_index/20000, y = processed_bead_1, color = population),
             size = 0.25,
             alpha = 0.5,
             shape = 16)+
  coord_cartesian(xlim = c(-0.043, 1.8))+
  draw_line(x = c(0, 0.2), y = -80)+
  annotate("text", x = 0.1, y = -85, label = "0.2 s", size = 8/.pt)+
  draw_line(x = -0.05, y = c(-37, -17))+
  annotate("text", x = -0.1, y = -27, label = "20 nm", angle = 90, size = 8/.pt)+
  ggtitle("10 &micro;M danicamtiv")+
  scale_color_manual(values = c("black", rep("#e7298a", 100)))+
  theme_void(10)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold", size = 10))

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


#only load svg one or it will crash R sessions
gg_trap <- image_read_svg("figures/trap-cartoon.svg", height = 10000, width = 10000)
gg_trap <- ggplot()+draw_image(gg_trap)+theme_nothing()
traces_right <- plot_grid(con_trace, dani_trace, nrow = 2)


## (top <- plot_grid(gg_trap, traces_right, ggton,
##                   rel_widths = c(1, 1.2, 0.75),
##                   labels = c("A", "B", "C"),
##                   nrow = 1))



cor_dat <- fread("dani-dose-motility-speed.csv")
step_cor <- cor.test(cor_dat$speed_nm_s, cor_dat$step_nm)
lmod2 <- lm(step_nm ~ speed_nm_s, data = cor_dat)
pline2 <- predict(lmod2, newdata = data.frame(speed_nm_s = seq(100, 300, by = 1)))
pdat2 <- data.frame(x = seq(100, 300, by = 1),
                   y = pline2)

(
gg_step_speed_cor <-
  ggplot(cor_dat, aes(x = speed_nm_s, y = step_nm))+
  geom_line(data = pdat2, aes(x, y))+
  geom_point(aes(color = id), size = 7)+
  geom_richtext(aes(label = num),
                size = 6/.pt,
                fill = NA,
                label.color = NA, # remove background and outline
                label.padding = grid::unit(rep(0, 4), "pt"))+ # remove padding )+
  annotate("richtext", x = Inf, y = -Inf, label = paste0("<br> r = ", round(step_cor$estimate, 2), "<br>",
                                                         "p = ", round(step_cor$p.value, 2)),
           hjust = 1,
           vjust = 0,
           label.color = NA,
           size = 8/.pt)+
  scale_y_continuous(breaks = 0:6)+
  coord_cartesian(xlim = c(0, 300), ylim = c(0, 5.5))+
  scale_color_manual(values = colorz)+
  ylab("Working Stroke (nm)")+
  xlab("Motility Speed (nm/s)")+
  theme_cowplot(10)+
  theme(
    axis.title = element_markdown(size = 10),
    legend.position = "none")
)


read_mtrackj <- function(folder, id_split, type = "Tracks"){
### helper function to read in motility track data
  read_mot <- function(x){
    d <- fread(x)
    d$path <- x
    d
  }
### pCa curves
  paths <- list.files(folder,
                      pattern = paste0("MTrackJ: ", type, ".csv"),
                      recursive = TRUE,
                      full.names = FALSE)
#
  paths <- file.path(folder, paths)
#
  data <- rbindlist(lapply(paths, read_mot))
### split the file path into new columns to get out conditions metadata
  data[, (id_split) := tstrsplit(path, "/", fixed = TRUE)]
}

mot_data <- read_mtrackj(folder = "motility/unreg-dani-drug-curve",
                     id_split = c("mot", "folder", "date","unreg", "id", "video", "filename"),
                     type = "Tracks")

mot_data$id <- sub("uM-dani", " <br> &micro;M", mot_data$id)
mot_data$id <- factor(mot_data$id, levels = c("0 <br> &micro;M", "1 <br> &micro;M", "10 <br> &micro;M", "100 <br> &micro;M"))

gg_unreg_dani_curve <-
ggplot()+
  geom_quasirandom(data = mot_data, aes(id, `Mean v [nm/sec]`, color = id),size = 1, width = 0.25)+
  coord_cartesian(ylim = c(-0.1, NA))+
  ylab("Speed (nm/s)")+
  xlab("[danicamtiv]")+
  ggtitle("Motility")+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  theme_cowplot(10)+
  theme(
    plot.title = element_markdown(size = 9, hjust = 0.5),
    legend.position = "none",
   axis.text.x = element_markdown(),
   axis.title =  element_markdown(size = 10)
   )

## mot_titration <- readRDS("figures/gg_mot_unreg_dani_titration.rds")+
##   ggtitle("Motility")+
##   scale_color_manual(values = colorz)+
##   theme_cowplot(11)+
##   theme(plot.title = element_markdown(size=9, hjust=0.5),
##         axis.text = element_markdown(size = 9),
##         axis.title = element_markdown(10),
##         legend.position = "none")
## trap_cor <- readRDS("figures/trap-mot-cor.rds")

top <- plot_grid(gg_trap, traces_right, rel_widths = c(1, 1.35), labels = c("A", "B"))

(middle <- plot_grid(ggea, ggton, nrow = 1, labels = c("C", "D"), rel_widths = c(1, 1.15)))

(bottom <- plot_grid(gg_unreg_dani_curve, gg_step_speed_cor,
                     rel_widths = c(1, 1),
                     nrow = 1,
                     labels = c("E", "F")))

png("figures/figure-2-new.png", width = 4.5, height = 6.25, units = "in", res = 300)
plot_grid(top, middle, bottom, nrow = 3, rel_heights = c(1.15,1, 1.1))
dev.off()

pdf("figures/figure-2-new.pdf", width = 4.5, height = 6.25)
plot_grid(top, middle, bottom, nrow = 3, rel_heights = c(1.15,1, 1.1))
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
