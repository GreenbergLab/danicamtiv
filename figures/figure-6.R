library(ggplot2)
library(cowplot)
library(ggtext)

theme_set(theme_cowplot(11))

colorz <- c("#666666",
            "#f7b8d8",
            "#ef72b2",
            "#e7298a")


## gg_reg_speed <- readRDS("figures/gg_reg_speed.rds")
gg_reg_speed <- readRDS("figures/gg_mot_regulation.rds")
## gg_reg_norm_speed <- readRDS("figures/gg_reg_norm_speed.rds")

(gg_pca625_dani_titration <- readRDS("figures/gg_pca625_dani_titration.rds")+annotate("text", x = 2.5, y = Inf,
                                                                                     label = "pCa 6.25",
                                                                                     hjust = 0.5,
                                                                                     vjust = 2)+
   scale_color_manual(values = colorz)
  )

gg_myo_conc_speed <- readRDS("figures/gg_myo_conc_speed.rds")
gg_myo_conc_norm_speed <- readRDS("figures/gg_myo_conc_speed_norm.rds")+scale_y_continuous(breaks = seq(0, 1, by=0.25))
## gg_unreg_dani_curve <- readRDS("figures/gg_mot_unreg_dani_titration.rds")+theme(axis.text.x = element_markdown(size = 8))

fibersim_force_pca <- readRDS("figures/fibersim-force-pca.rds")
## fibersim_twitch1 <- readRDS("figures/fibersim-ggtwitch1.rds")
## fibersim_twitch2 <- readRDS("figures/fibersim-ggtwitch2.rds")
fibersim_twitch_integral <- readRDS("figures/fibersim-twitch-integral.rds")
fibersim_twitch <- readRDS("figures/fibersim-ggtwitch.rds")

title <- ggdraw() +
  draw_label(
    "Unregulated Actin Filaments",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5
  )+
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    ## plot.margin = margin(0, 0, 0, 200)
  )

(top0 <- plot_grid(gg_myo_conc_speed, gg_myo_conc_norm_speed,
                  labels = c("A", "B"), nrow = 1, rel_widths = c(1, 1)))

(top <- plot_grid(title, top0, nrow = 2, rel_heights = c(0.1, 0.9)))


title2 <- ggdraw() +
  draw_label(
    "Experimental Regulated Motility",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5,
    size = 10,
  )+
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    ## plot.margin = margin(0, 0, 0, 200)
  )

(middle0 <- plot_grid(gg_reg_speed, gg_pca625_dani_titration,
                  labels = c("A", "B"), nrow = 1))


(middle <- plot_grid(title2, middle0, nrow = 2, rel_heights = c(0.1, 0.9)))

twitch2 <- fibersim_twitch+draw_plot(fibersim_twitch_integral, x = 0.55, 0.3, 0.3, 0.8)


title3 <- ggdraw() +
  draw_label(
    "Computer Simulations",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5,
    size = 10
  )+
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    ## plot.margin = margin(0, 0, 0, 200)
  )

(bottom0 <- plot_grid(fibersim_force_pca, fibersim_twitch, fibersim_twitch_integral,
                     labels = c("F", "G", "H"),
                     nrow = 1,
                     rel_widths = c(1, 1, 0.5)))

(bottom <- plot_grid(title3, bottom0, nrow = 2, rel_heights = c(0.1, 0.9)))
## bottom <- plot_grid(gg_reg_speed, gg_reg_norm_speed, labels = c("C", "D"))
## top <- plot_grid(gg_myo_conc_speed, gg_myo_conc_norm_speed, labels = c("A", "B"))


png("figures/figure-5-new.png", width = 6.5, height = 5, res = 300, units = "in")
plot_grid(top, middle, nrow = 2, rel_heights = c(1, 1, 1))
dev.off()


pdf("figures/for-the-boss2/figure-6.pdf", width = 6.5, height = 6.5)
plot_grid(top, middle, bottom, nrow = 3, rel_heights = c(1, 1, 1))
dev.off()
