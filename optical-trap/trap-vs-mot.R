library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
library(cocor)

dat <- fread("dani-dose-motility-speed.csv")

lmod <- lm(step_rel ~ rel_con, data = dat)
lmod2 <- lm(step_nm ~ speed_nm_s, data = dat)
lmod3 <- lm(ton ~ speed_nm_s, data = dat)

pline <- predict(lmod, newdata = data.frame(rel_con = seq(0, 1, by = 0.01)))
pline2 <- predict(lmod2, newdata = data.frame(speed_nm_s = seq(0, 300, by = 1)))
pline3 <- predict(lmod3, newdata = data.frame(speed_nm_s = seq(100, 300, by = 1)))

pdat <- data.frame(x = seq(0, 1, by = 0.01),
                   y = pline)


pdat2 <- data.frame(x = seq(0, 300, by = 1),
                   y = pline2)


pdat3 <- data.frame(x = seq(100, 300, by = 1),
                   y = pline3)

ggplot(dat, aes(x = rel_con, y = step_rel))+
  geom_line(data = pdat, aes(x, y))+
  geom_point(aes(color = id))+
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))


colorz <- c("#666666",
            "#f7b8d8",
            "#ef72b2",
            "#e7298a")


step_cor <- cor.test(dat$speed_nm_s, dat$step_nm)
ton_cor <- cor.test(dat$speed_nm_s, dat$ton)

(
gg1 <-
ggplot(dat, aes(x = speed_nm_s, y = step_nm))+
  ## geom_line(data = pdat2, aes(x, y))+
  geom_point(aes(color = id), size = 7)+
  geom_richtext(aes(label = num),
                size = 2,
                fill = NA,
                label.color = NA, # remove background and outline
    label.padding = grid::unit(rep(0, 4), "pt"))+ # remove padding )+
  annotate("richtext", x = Inf, y = -Inf, label = paste0("<br> r = ", round(step_cor$estimate, 2), "<br>",
                                                         "p = ", round(step_cor$p.value, 2)),
           hjust = 1,
           vjust = 0,
           label.color = NA,
           size = 3)+
  scale_y_continuous(breaks = 0:6)+
  coord_cartesian(xlim = c(0, 300), ylim = c(0, 5.5))+
  scale_color_manual(values = colorz)+
  ylab("Working Stroke (nm)")+
  xlab("Motility Speed (nm/s)")+
  theme_cowplot(11)+
  theme(
    axis.title.y = element_markdown(size = 9),
  legend.position = "none")
  )

(
gg2 <-
ggplot(dat, aes(x = speed_nm_s, y = ton))+
  geom_line(data = pdat3, aes(x, y))+
  geom_point(aes(color = id), size = 7)+
  geom_richtext(aes(label = num),
                size = 2,
                fill = NA,
                label.color = NA, # remove background and outline
    label.padding = grid::unit(rep(0, 4), "pt"))+ # remove padding )+
  annotate("richtext", x = Inf, y = -Inf, label = paste0("<br> r = ", round(ton_cor$estimate, 2), "<br>",
                                                         "p = ", round(ton_cor$p.value, 2)),
           hjust = 1,
           vjust = 0,
           label.color = NA,
           size = 3)+
  ## scale_y_continuous(breaks = 0:6)+
  coord_cartesian(xlim = c(0, 300), ylim = c(0, 25))+
  scale_color_manual(values = colorz)+
  ylab("Detachment Rate (s<sup>-1</sup>)")+
  xlab("Motility Speed (nm/s)")+
  theme_cowplot(11)+
  theme(
  legend.position = "none",
  axis.title.y = element_markdown())
  )

## ggsave("supplemental-figures/ton-vs-mot-speed.png", bg = "white")
png("supplemental-figures/ton-vs-mot-speed.png", height = 2.5, width = 3, units = "in", res=300)
gg2
dev.off()


pdf("supplemental-figures/ton-vs-mot-speed.pdf", height = 2.5, width = 3)
gg2
dev.off()


(gg_together <- plot_grid(gg1, gg2, labels = c("F", "G"), nrow = 1))

saveRDS(gg_together, "figures/trap-mot-cor.rds")


cocor(~speed_nm_s + step_nm | speed_nm_s + ton, data = list(dat, dat))

cor(dat$speed_nm_s, dat$step_nm)
cor(dat$speed_nm_s, dat$ton)


