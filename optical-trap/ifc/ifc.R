library(data.table)
library(ggplot2)
library(grid)
library(cowplot)
library(ggtext)

ifc <- fread("optical-trap/lasertrapr/project_dani-ifc/summary/2024-06-07_project_dani-ifc_all-measured-events.csv")


gg <- readRDS("optical-trap/lasertrapr/project_dani-ifc/summary/figures/2024-06-07_ifc-bell-plots.rds")

gg$layers[[1]]$aes_params$size <- 1
gg$layers[[1]]$aes_params
gg$layers[[2]]$aes_params$size <-3

## png("figures/ifc.png", width = 6.5, height = 3, units = "in", res = 300)
## gg+
##   scale_color_manual(values = c("#666666", "#e7298a"))+
##   theme_cowplot(11)+
##   theme(strip.background = element_blank(),
##         strip.text = element_blank(),
##         legend.position = "none")
## dev.off()
gg_b <- ggplot_build(gg)

gg_b$data[[2]]$label <- sub("10uM-dani", "10 &micro;M danicamtiv", gg_b$data[[2]]$label)
gg_b$data[[2]]$label <- sub("control", "0 &micro;M danicamtiv", gg_b$data[[2]]$label)

gg_c <- ggplot_gtable(gg_b)

(gg <- gg+
  ## scale_color_manual(values = c("#666666", "#e7298a"))+
   facet_wrap(~conditions, scales ="free_y")+
  theme_cowplot(11)+
   scale_x_continuous(breaks = seq(0, 10, 2))+
  coord_cartesian(ylim = c(0.0029, 10), xlim = c(0, 10))+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")
  )


gg1 <-
ggplot(data = ifc[force_pn >= 0 & conditions == "control"], aes(x = force_pn, y = time_on_s ))+
  geom_point(color = "#666666", alpha = 0.5, size = 1)+
  annotate("richtext",
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           size = 3,
           label = gg_b$data[[2]]$label[[2]],
           label.color = "white",
           color = "#666666"
           )+
  xlab("Force (pN)")+
  ylab("Attachment Time (s)")+
  scale_y_log10(expand = expansion(c(0.17, 0.1), c(0.65, 0.1)),
                breaks = c(0.01, 0.1, 1, 10))+
   scale_x_continuous(breaks = seq(0, 10, 2))+
  coord_cartesian(ylim = c(0.01, 10), xlim = c(0, 10))+
  theme_cowplot(11)


gg2 <-
  ggplot(data = ifc[force_pn >= 0 & conditions == "10uM-dani"], aes(x = force_pn, y = time_on_s ))+
  geom_point(color = "#e7298a", alpha = 0.5, size = 1)+
  annotate("richtext",
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           size = 3,
           label = gg_b$data[[2]]$label[[1]],
           label.color = "white",
           color = "#e7298a"
           )+
  xlab("Force (pN)")+
  ylab("Attachment Time (s)")+
  scale_y_log10(expand = expansion(c(0.17, 0.1), c(0.65, 0.1)),
                breaks = c(0.01, 0.1, 1, 10))+
   scale_x_continuous(breaks = seq(0, 10, 2))+
  coord_cartesian(ylim = c(0.01, 10), xlim = c(0, 10))+
  theme_cowplot(11)



bottom_new <- plot_grid(gg1, gg2, labels = c("B", ""))

title_panel <- ggdraw() +
  draw_label(
    "Load Dependent Detachment Kinetics",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5, size = 10
  )+
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    ## plot.margin = margin(0, 0, 0, 200)
  )




## gg

## gg_c <- gg_b$plot+facet_wrap(~conditions, scales = "free_y")

## predict_df <- fread("~/lasertrapr/project_dani-ifc/summary/2024-06-07_ifc-predict-line.csv")
## ribbon_df <- fread("~/lasertrapr/project_dani-ifc/summary/2024-06-07_ifc-ci-ribbon.csv")

## predict_df$conditions <- factor(predict_df$conditions, levels = c("control", "10uM-dani"))
## ribbon_df$conditions <- factor(ribbon_df$conditions, levels = c("control", "10uM-dani"))

## bell_fits <-
## ggplot()+
##   geom_line(data = predict_df[F<=8], aes(F, t, color = conditions), linewidth = 1)+
##   geom_ribbon(data = ribbon_df[F<=8], aes(x = F, ymin = tmin, ymax = tmax, fill = conditions), alpha = 0.4)+
##   scale_fill_manual(values = c("#666666", "#e7298a"))+
##   scale_color_manual(values = c("#666666", "#e7298a"))+
##   ylab("")+
##   xlab("Force (pN)")+
##   scale_y_log10()+
##   coord_cartesian(ylim = c(0.0029, 10))+
##   theme_cowplot(11)+
##   theme(
##     axis.title.y = element_blank(),
##     axis.line.y = element_blank(),
##     axis.ticks.y = element_blank(),
##     axis.text.y = element_blank(),
##     legend.position = "none"
##     )

## (
## bottom <- plot_grid(gg_c, bell_fits, rel_widths = c(3, 1), labels = c("C", "D"), align = "h", axis = "b")
## )

(con_trace <- readRDS("optical-trap/lasertrapr/project_dani-ifc/summary/figures/control_2023-03-14_obs-02_9.2443-10.4391.rds"))

con_options <- fread("optical-trap/lasertrapr/project_dani-ifc/control/2023-03-14/obs-02/options.csv")
con_nm2pn <- con_options$nm2pn2

con_trace <-
con_trace+
  draw_line(c(0, 0.2), y = c(-150, -150))+
  annotate("text", x = 0.1, y = -180, label = "0.2 s", size = 3)+
  draw_line(-0.01, c(200, 300))+
  annotate("text", x = -0.04, y = 250, label = paste0(round(100*con_nm2pn, 1), " pN"), size = 3, angle = 90)+
  annotate("text", x = 1.225, y = 200, label = "M", size = 3)+
  annotate("text", x = 1.225, y = -50, label = "T", size = 3)+
  theme_void()
con_trace

dani_trace <- readRDS("optical-trap/lasertrapr/project_dani-ifc/summary/figures/10uM-dani_2023-11-20_obs-02_151.8853-152.926.rds")

dani_options <- fread("optical-trap/lasertrapr/project_dani-ifc/10uM-dani/2023-11-20/obs-02/options.csv")
dani_nm2pn <- dani_options$nm2pn

dani_trace <-
dani_trace+
  draw_line(c(0, 0.2), y = c(-190, -190))+
  annotate("text", x = 0.1, y = -215, label = "0.2 s", size = 3)+
  draw_line(-0.01, c(100, 183))+
  annotate("text", x = -0.04, y = 142.5, label = paste0(round(83*dani_nm2pn, 1), " pN"), size = 3, angle = 90)+
  annotate("text", x = 1.07, y = 100, label = "M", size = 3)+
  annotate("text", x = 1.07, y = -90, label = "T", size = 3)+
  theme_void()

gg_traces <- plot_grid(con_trace, dani_trace, labels = c("A", ""), nrow = 1)


png("figures/figure-3-new.png", width = 6.5, height = 4.2, units = "in", res = 300)
plot_grid(title_panel, gg_traces, bottom_new, nrow = 3, rel_heights = c(0.25, 1.1, 2))
dev.off()


pdf("figures/for-the-boss2/figure-4-new.pdf", width = 6.5, height = 4)
plot_grid(gg_traces, gg_c, nrow = 2, rel_heights = c(1, 2), labels = c("A", "B"))
dev.off()


pdf("figures/for-the-boss/ifc.pdf", width = 6.5, height = 4)
plot_grid(gg_traces, gg_c, nrow = 2, labels = c("", "C"), rel_heights = c(1, 2))
dev.off()


svg("figures/for-the-boss/ifc.svg", width = 6.5, height = 4)
plot_grid(gg_traces, gg_c, nrow = 2, labels = c("", "C"), rel_heights = c(1, 2))
dev.off()


##export for memlet

dat <- fread("~/lasertrapr/project_dani-ifc/summary/2024-06-07_ifc-bootstrap-bell-parameters.csv")

con_dat <- dat[conditions == "control", -1]
dani_dat <- dat[conditions == "10uM-dani", -1]

fwrite(con_dat, "optical-trap/ifc/control-ifc-bootstrap.csv", col.names = FALSE)
fwrite(dani_dat, "optical-trap/ifc/dani-ifc-bootstrap.csv", col.names = FALSE)


opt <- list.files("optical-trap/lasertrapr/project_dani-ifc",
           pattern = "options.csv",
           recursive = TRUE,
           full.names = TRUE)

opt_df <- rbindlist(lapply(opt, fread), fill = TRUE)

opt_df <- opt_df[include == TRUE & review == TRUE & status == "analyzed"]

con_opt <- opt_df[conditions == "control"]

unique(con_opt$date)


dani_opt <- opt_df[conditions == "10uM-dani"]


unique(dani_opt$date)
