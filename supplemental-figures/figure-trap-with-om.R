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
            "#e7298a",
            "#d95f02")


spasm_files <- list("0 &micro;M" = "optical-trap/spasm-analysis/control/spasm-included.xlsx",
                    "1 &micro;M" = "optical-trap/spasm-analysis/1uM-dani/spasm_1uM-dani_10uM-atp.xlsx",
                    "10 &micro;M" = "optical-trap/spasm-analysis/10uM-dani/spasm-included-manual.xlsx",
                    "100 &micro;M" = "optical-trap/spasm-analysis/100uM-dani/spasm_100uM-dani_10uM-atp.xlsx",
                    "10 &micro;M OM" = "optical-trap/spasm-analysis/10uM-om/spasm-included.xlsx"
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

om <- spasm_ensemble_average(spasm_file = "optical-trap/spasm-analysis/10uM-om/combinedEnsembleAxesData.xlsx",
                              forward_time_filter = 0.3,
                              reverse_time_filter = 0.5,
                              step_estimation_method = "fit",
                              color = colorz[3],
                              textsize = 8,
                              ylim = c(-2, 9),
                              x_shift = 0.9,
                              title = "10 &micro;M OM"
                              )

con_f_final_nm <- mean(tail(con$ea_f$nanometers), 200)
con_f_final_nm
mean(tail(dani_1uM$ea_f$nanometers), 200)/con_f_final_nm
mean(tail(dani$ea_f$nanometers), 200) / con_f_final_nm
mean(tail(dani_100uM$ea_f$nanometers), 200) /con_f_final_nm


dat_labels <- data.frame(y = c(mean(tail(con$ea_f$nanometers, 200)),
                                   mean(tail(dani_1uM$ea_f$nanometers), 200) ,
                                   mean(tail(dani$ea_f$nanometers), 200)+0.06,
                                   mean(tail(dani_100uM$ea_f$nanometers), 200)-0.085,
                                   mean(tail(om$ea_f$nanometers), 200)),
                         x = 0.355,
                         id = c("0 &micro;M Dani", "1 &micro;M Dani", "10 &micro;M Dani", "100 &micro;M Dani", "10 &micro;M OM"))


con_ea_f <- con$ea_f
con_ea_f$id <- "0 &micro;M Dani"

dani_ea_f <- dani$ea_f
dani_ea_f$id <- "10 &micro;M Dani"


dani_1uM_ea_f <- dani_1uM$ea_f
dani_1uM_ea_f$id <- "1 &micro;M Dani"


dani_100uM_ea_f <- dani_100uM$ea_f
dani_100uM_ea_f$id <- "100 &micro;M Dani"


om_ea_f <- om$ea_f
om_ea_f$id <- "10 &micro;M OM"

ea_dat <- rbind(con_ea_f, dani_ea_f, dani_1uM_ea_f, dani_100uM_ea_f, om_ea_f)

ea_dat$id <- factor(ea_dat$id, levels = c("0 &micro;M Dani", "1 &micro;M Dani", "10 &micro;M Dani", "100 &micro;M Dani", "10 &micro;M OM"))

con_f_pred <- con$f_pred
dani_f_pred <- dani$f_pred
dani_1uM_f_pred <- dani_1uM$f_pred
dani_100uM_f_pred <- dani_100uM$f_pred
om_f_pred <- om$f_pred

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
  scale_color_manual(values = colorz)+
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




ton <- spasm_fit_ton(spasm_data = spasm_data, tmin = 0.01, colorz = colorz)

label_ton <- paste0("<span style='color:#666666'>0 &micro;M Dani = ", ton$boot_data$html_label[[1]],"  </span><br>",
                    "<span style='color:#f7b8d8'> 1 &micro;M Dani = ", ton$boot_data$html_label[[2]]," </span><br>",
                    "<span style='color:#ef72b2'> 10 &micro;M Dani = ", sub("24", "24.0", ton$boot_data$html_label[[3]]), " </span><br>",
                    "<span style='color:#e7298a'> 100 &micro;M Dani = ", ton$boot_data$html_label[[4]],"  </span><br>",
                    "<span style='color:#d95f02'> 10 &micro;M OM = ", ton$boot_data$html_label[[5]],"  </span><br>" )


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





png("supplemental-figures/trap-with-om.png", width = 6.5, height = 3, units = "in", res = 300)
plot_grid(ggea, ggton, labels = LETTERS)
dev.off()


cairo_pdf("supplemental-figures/trap-with-om.pdf", width = 6.5, height = 3)
plot_grid(ggea, ggton, labels = LETTERS)
dev.off()
