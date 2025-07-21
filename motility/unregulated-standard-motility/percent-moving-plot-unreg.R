library(data.table)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(ungeviz)
library(ggtext)

dat <- fread("motility/unregulated-standard-motility/percent-moving-hand-tracked.csv")
dat <- na.omit(dat)

dat[, percent_moving := number_moving/total_filaments]

dat_sum <- dat[, .(avg_percent_moving = mean(percent_moving),
                   sd_percent_moving = sd(percent_moving)),
               by = .(conditions)]

ttest <- t.test(dat[conditions == "0 &micro;M"]$percent_moving,
                dat[conditions == "10 &micro;M"]$percent_moving
                )

gg <-
ggplot()+
  geom_hpline(data = dat_sum,
              aes(conditions,
                  avg_percent_moving,
                  color = conditions),
              width = 0.3,
              alpha = 0.5)+
  geom_quasirandom(data = dat,
                   aes(conditions,
                       percent_moving,
                       color = conditions),
                   size = 1,
                   width = 0.2)+
  draw_line(c(1, 2), c(1.15, 1.15))+
  annotate("text", x = 1.5, y = 1.2, label = paste0("p = ", round(ttest$p.value, 2)), size = 3)+
  ylab("Percent Moving (%)")+
  xlab("[Danicamtiv]")+
  scale_color_manual(values = c("#666666", "#e7298a"))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  coord_cartesian(ylim = c(0, 1.2))+
  theme_cowplot(11)+
  theme(
    legend.position = "none",
    axis.text.x = element_markdown()
  )


## ggsave("supplemental-figures/unregulated-percent-moving.png", bg = "white")

png("supplemental-figures/unreg-percent-moving-update.png", width = 3, height = 2.5, units = "in", res = 300)
gg
dev.off()


pdf("supplemental-figures/unreg-percent-moving-update.pdf", width = 3, height = 2.5)
gg
dev.off()
