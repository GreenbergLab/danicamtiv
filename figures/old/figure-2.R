library(ggplot2)
library(cowplot)
library(ggtext)
library(patchwork)

theme_set(theme_cowplot(11))
gg_adp <- readRDS("figures/gg_adp.rds")+ylab("Normalized Fluorescence")
gg_adp_affinity <- readRDS("figures/gg_adp_affinity.rds")+scale_y_continuous(breaks = seq(0, 75,by =  25))
gg_kfast <- readRDS("figures/gg_fkast.rds")+scale_y_continuous(breaks = seq(0, 1200, by = 300))
gg_hydro <- readRDS("figures/gg-hydro.rds")

top <- plot_grid(gg_kfast, gg_adp, align = "h", labels = c("A", "B"))
bottom <- plot_grid(gg_adp_affinity, gg_hydro, align = "h", labels = c("C", "D"))

left <- plot_grid(gg_adp, gg_kfast, labels = c("A", "C"), nrow = 2, align = "v")
right <- plot_grid(gg_adp_affinity, gg_hydro,labels = c("B", "D"), nrow = 2, align = "v")

png("figures/for-the-boss2/figure-2.png", width = 6.5, height = 5.25, units = "in", res = 300)
plot_grid(top, bottom, nrow = 2)
dev.off()

pdf("figures/for-the-boss2/figure-2.pdf", width = 6.5, height = 5.25)
plot_grid(top, bottom, nrow = 2)
dev.off()

postscript("figures/for-the-boss2/figure-2.eps", width = 6.5, height = 5.25)
plot_grid(top, bottom, nrow = 2)
dev.off()

## png("figures/figure-2.png", width = 6.5, height = 5.5, units = "in", res = 300)
## plot_grid(left, right, nrow = 1, align = "hv")
## dev.off()



## png("figures/figure-2.png", width = 6.5, height = 5.25, units = "in", res = 300)
## (gg_adp+gg_adp_affinity)/(gg_kfast+gg_hydro)+plot_annotation(tag_levels = "A")
## dev.off()

gg_kfast+
  theme_cowplot(12, font_family = "Arial")+
  scale_y_continuous(breaks = seq(0, 1250, 250))+
  theme(
    axis.title.y = element_markdown(),
    legend.position = "none")

ggsave("~/Downloads/kfast-all-days.pdf", device = cairo_pdf)
