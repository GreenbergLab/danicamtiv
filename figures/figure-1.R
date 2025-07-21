library(ggplot2)
library(cowplot)
library(ggtext)
library(magick)
library(rsvg)
library(dplyr)

theme_set(
  theme_cowplot(10)+
  theme(plot.title = element_markdown(hjust=0.5, size = 9),
        legend.position = "none")
)

gg_atpase <-
  readRDS("figures/gg_atpase.rds")+
  ggtitle("ATPase")+
  ylab("Rate (head<sup>-1</sup>&middot;s<sup>-1</sup>)")+
  theme_cowplot(10)+
          theme(
            plot.title = element_markdown(hjust=0.5, size = 9),
                axis.title = element_markdown(),
                plot.margin = margin(r=20,l=20, t=5),
                legend.position = "none")

gg_mot <- readRDS("figures/gg_mot.rds")+
  ggtitle("Motility")+
  ylab("Speed (nm/s)")+
  theme_cowplot(10)+
  theme(
    plot.title = element_markdown(hjust=0.5, size = 9),
        axis.title = element_markdown(),
    axis.text.x = element_markdown(),
        plot.margin = margin(r=20,l=20, t=5),
        legend.position = "none")

gg_adp <- readRDS("figures/gg_adp.rds")+
  ylab("Normalized <br> Fluorescence")+
  xlab("Time (ms)")+
  ggtitle("ADP Release <br> ( <i>k<sub>+5</sub>'</i> )")+
scale_x_continuous(labels = \(x) x*1000, breaks = c(0, 0.1, 0.2))+
  theme_cowplot(10)+
  theme(plot.title = element_markdown(hjust=0.5, size=9),
        axis.title = element_markdown(size=10),
        legend.position = "none")

gg_adp_affinity <- readRDS("figures/gg_adp_affinity.rds")+
  scale_y_continuous(breaks = seq(0, 75,by =  25))+
  ggtitle("ADP Affinity <br> ( K<sub>5</sub>' )")+
  xlab("ADP (&micro;M)")+
  theme_cowplot(10)+
  theme(plot.title = element_markdown(hjust=0.5, size = 9),
        axis.title = element_markdown(),
        legend.position = "none")


(
gg_kfast <- readRDS("figures/gg_fkast.rds")+
  ggtitle("ATP Induced <br> Dissociation <br> ( <i>k<sub>+2</sub>'</i> )")+
  xlab("ATP (mM)")+
  scale_x_continuous(labels = \(x) x/1000)+
  scale_y_continuous(breaks = seq(0, 1200, by = 300))+
  theme_cowplot(10)+
  theme(plot.title = element_markdown(hjust=0.5, size = 9),
        axis.title = element_markdown(size = 10),
        legend.position = "none")
  )

## gg_hydro <- readRDS("figures/gg-hydro.rds")

transient_plots <- plot_grid(gg_adp, gg_adp_affinity, gg_kfast, labels = c("D", "E", "F"), nrow = 1)
transient_plots2 <- plot_grid(gg_kfast, gg_adp, gg_adp_affinity, labels = c("D", "E", "F"), nrow = 1)
## bottom <- plot_grid(gg_adp_affinity, gg_hydro, align = "h", labels = c("C", "D"))
## stead_state_plots <- plot_grid(gg, mod,gg_atpase, gg_mot, labels = c("B", "C"))



## mod <- image_read_svg("figures/atpase-cycle-cartoon.svg", width = 5000, height = 2000)
mod <- image_read_svg("figures/scheme.svg", width = 5000, height = 2000)
## mod <- image_read_svg("figures/scheme2.svg", width = 5000, height = 2000)

(gg_mod <-
  ggplot()+
  draw_image(mod, clip = F, scale = 2)+
  theme_nothing()+
  theme(
    ## plot.background = element_rect(size = 1, color = "black"),
     plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)
                                        )# Left margin
  )

## top <- plot_grid(gg_atpase, gg_mot, labels = c("A", "B"))

stead_state_plots <- plot_grid(gg_atpase, gg_mot, labels = c("A", "B"), nrow = 1,
                               rel_widths = c(1,1))


png("figures/figure-1-new.png", width = 4.5, height = 4.5, units = "in", res = 300)
plot_grid(stead_state_plots, gg_mod,transient_plots2, ncol = 1, labels = c("", "C", ""),
          rel_heights = c(1, 0.55, 1))
dev.off()


svg("figures/figure-1-new.svg", width = 4.5, height = 4.5)
plot_grid(stead_state_plots, gg_mod,transient_plots2, ncol = 1, labels = c("", "C", ""),
          rel_heights = c(1, 0.55, 1))
dev.off()



pdf("figures/figure-1-new.pdf", width = 4.5, height = 4.5)
plot_grid(stead_state_plots, gg_mod,transient_plots2, ncol = 1, labels = c("", "C", ""),
          rel_heights = c(1, 0.55, 1))
dev.off()



cairo_ps("figures/figure-1-new.eps", width = 4.5, height = 4.5, fallback_resolution = 1000)
plot_grid(stead_state_plots, gg_mod,transient_plots2, ncol = 1, labels = c("", "C", ""),
          rel_heights = c(1, 0.55, 1))
dev.off()



cairo_pdf("figures/figure-1-new.pdf", width = 4.5, height = 4.5)
plot_grid(stead_state_plots, gg_mod,transient_plots2, ncol = 1, labels = c("", "C", ""),
          rel_heights = c(1, 0.55, 1))
dev.off()



png("figures/figure-1-alt.png", width = 4.5, height = 4.5, units = "in", res = 300)
plot_grid(stead_state_plots, gg_mod,transient_plots, ncol = 1, labels = c("", "C", ""),
          rel_heights = c(1, 0.55, 1))
dev.off()



cairo_pdf("figures/figure-1-alt.pdf", width = 4.5, height = 4.5)
plot_grid(stead_state_plots, gg_mod,transient_plots, ncol = 1, labels = c("", "C", ""),
          rel_heights = c(1, 0.55, 1))
dev.off()

cairo_pdf("figures/for-the-boss2/figure-1.pdf", width = 6.5, height = 6)
plot_grid(top, gg_mod, nrow = 2, labels = c("", "C"), rel_heights = c(1, 1))
## plot_grid(gg_mod, bottom, nrow = 2, labels = c("A", ""), rel_heights = c(1, 0.75))
dev.off()

postscript("figures/for-the-boss2/figure-1.eps", width = 6.5, height = 6)
plot_grid(top, gg_mod, nrow = 2, labels = c("", "C"), rel_heights = c(1, 1))
## plot_grid(gg_mod, bottom, nrow = 2, labels = c("A", ""), rel_heights = c(1, 0.75))
dev.off()
