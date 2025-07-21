library(data.table)
library(ggplot2)
library(grid)
library(cowplot)
library(ggtext)

ifc <- fread("optical-trap/ifc/lasertrapr/project_dani-ifc/summary/2024-06-07_project_dani-ifc_all-measured-events.csv")


gg <- readRDS("optical-trap/ifc/lasertrapr/project_dani-ifc/summary/figures/2024-06-07_ifc-bell-plots.rds")

gg$layers[[1]]$aes_params$size <- 1
gg$layers[[1]]$aes_params
gg$layers[[2]]$aes_params$size <-3

gg_b <- ggplot_build(gg)

gg_b$data[[2]]$label <- sub("10uM-dani", "10 &micro;M danicamtiv", gg_b$data[[2]]$label)
gg_b$data[[2]]$label <- sub("control", "0 &micro;M danicamtiv", gg_b$data[[2]]$label)

gg_c <- ggplot_gtable(gg_b)

(gg <- gg+
   facet_wrap(~conditions, scales ="free_y")+
   theme_cowplot(10)+
   scale_x_continuous(breaks = seq(0, 10, 2))+
   coord_cartesian(ylim = c(0.0029, 10), xlim = c(0, 10))+
   theme(strip.background = element_blank(),
         strip.text = element_blank(),
         legend.position = "none")
)

gg1 <-
ggplot(data = ifc[force_pn >= 0 & conditions == "control"], aes(x = force_pn, y = time_on_s ))+
  annotate("richtext",
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           size = 8/.pt,
           label = gg_b$data[[2]]$label[[2]],
           label.color = "transparent",
           color = "#666666"
           )+
  geom_point(color = "#666666", alpha = 0.5, size = 1)+
  xlab("Force (pN)")+
  ylab("Attachment Time (s)")+
  scale_y_log10(expand = expansion(c(0.17, 0.1), c(0.65, 0.1)),
                breaks = c(0.01, 0.1, 1, 10))+
   scale_x_continuous(breaks = seq(0, 10, 2))+
  coord_cartesian(ylim = c(0.01, 10), xlim = c(0, 10))+
  theme_cowplot(10)


gg2 <-
  ggplot(data = ifc[force_pn >= 0 & conditions == "10uM-dani"], aes(x = force_pn, y = time_on_s ))+
  annotate("richtext",
           x = -Inf,
           y = Inf,
           hjust = 0,
           vjust = 1,
           size = 8/.pt,
           label = gg_b$data[[2]]$label[[1]],
           label.color = "transparent",
           color = "#e7298a"
           )+
  geom_point(color = "#e7298a", alpha = 0.5, size = 1)+
  xlab("Force (pN)")+
  ylab("Attachment Time (s)")+
  scale_y_log10(expand = expansion(c(0.17, 0.1), c(0.65, 0.1)),
                breaks = c(0.01, 0.1, 1, 10))+
   scale_x_continuous(breaks = seq(0, 10, 2))+
  coord_cartesian(ylim = c(0.01, 10), xlim = c(0, 10))+
  theme_cowplot(10)



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


(con_trace <- readRDS("optical-trap/ifc/lasertrapr/project_dani-ifc/summary/figures/control_2023-03-14_obs-02_9.2443-10.4391.rds"))

con_options <- fread("optical-trap/ifc/lasertrapr/project_dani-ifc/control/2023-03-14/obs-02/options.csv")
con_nm2pn <- con_options$nm2pn2

con_trace <-
con_trace+
  draw_line(c(0, 0.2), y = c(-150, -150))+
  annotate("text", x = 0.1, y = -180, label = "0.2 s", size = 3)+
  draw_line(-0.01, c(200, 300))+
  annotate("text", x = -0.065, y = 250, label = paste0(round(100*con_nm2pn, 1), " pN"), size = 3, angle = 90)+
  annotate("text", x = 1.235, y = 200, label = "M", size = 3)+
  annotate("text", x = 1.235, y = -50, label = "T", size = 3)+
  theme_void()
con_trace

dani_trace <- readRDS("optical-trap/ifc/lasertrapr/project_dani-ifc/summary/figures/10uM-dani_2023-11-20_obs-02_151.8853-152.926.rds")

dani_options <- fread("optical-trap/ifc/lasertrapr/project_dani-ifc/10uM-dani/2023-11-20/obs-02/options.csv")
dani_nm2pn <- dani_options$nm2pn

dani_trace <-
dani_trace+
  draw_line(c(0, 0.2), y = c(-190, -190))+
  annotate("text", x = 0.1, y = -215, label = "0.2 s", size = 3)+
  draw_line(-0.01, c(100, 183))+
  annotate("text", x = -0.06, y = 142.5, label = paste0(round(83*dani_nm2pn, 1), " pN"), size = 3, angle = 90)+
  annotate("text", x = 1.08, y = 100, label = "M", size = 3)+
  annotate("text", x = 1.08, y = -90, label = "T", size = 3)+
  theme_void()

gg_traces <- plot_grid(con_trace, dani_trace, labels = c("A", ""), nrow = 1)


png("figures/figure-3.png", width = 4.5, height = 4, units = "in", res = 300)
plot_grid(title_panel, gg_traces, bottom_new, nrow = 3, rel_heights = c(0.1, 0.75, 1))
dev.off()

pdf("figures/figure-3.pdf", width = 4.5, height = 4)
plot_grid(title_panel, gg_traces, bottom_new, nrow = 3, rel_heights = c(0.1, 0.75, 1))
dev.off()

