library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
library(ggbeeswarm)
library(ungeviz)

basesize <- 11

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

  paths <- file.path(folder, paths)

  data <- rbindlist(lapply(paths, read_mot))
### split the file path into new columns to get out conditions metadata
  data[, (id_split) := tstrsplit(path, "/", fixed = TRUE)]
}

data <- read_mtrackj(folder = "motility/unregulated-standard-motility",
                     id_split = c("mot", "unreg", "date", "id", "video", "filename"),
                     type = "Tracks")

data$id <- sub("uM-dani", " &micro;M", data$id)
data$id <- factor(data$id, levels = c("0 &micro;M", "10 &micro;M"))

## data[, conc_num := as.numeric(sub("ug-mL", "", conc))]


data_by_video <- data[, .(speed_nm_s = mean(`Mean v [nm/sec]`)), by = .(id, date, video)]
data_by_video[, relative_speed_nm_s := speed_nm_s/max(speed_nm_s), by = .(date, id, video) ]


## data_by_day <- data_by_video[, .(speed_nm_s = mean(`Mean v [nm/sec]`)), by = .(id, date)]
## data_by_day[, relative_speed_nm_s := speed_nm_s/max(speed_nm_s), by = .(date, id) ]

data_by_id <- data_by_video[, .(
  speed_nm_s = mean(speed_nm_s),
  sd_speed = sd(speed_nm_s),
  se_speed = sd(speed_nm_s)/sqrt(.N),
  relative_speed_nm_s = mean(relative_speed_nm_s),
  sd_rel_speed = sd(relative_speed_nm_s),
  se_rel_speed = sd(relative_speed_nm_s)/sqrt(.N)),
  by = .(id)]


colorz <- c("#666666", "#e7298a")

percent_drop <- round((data_by_id[id == "10 &micro;M"]$speed_nm_s/data_by_id[id == "0 &micro;M"]$speed_nm_s)*100,
                      0)


shapiro.test(data_by_video[id == "10 &micro;M"]$speed_nm_s)
shapiro.test(data_by_video[id == "0 &micro;M"]$speed_nm_s)

mot_pval <-
  t.test(data_by_video[id == "10 &micro;M"]$speed_nm_s, data_by_video[id == "0 &micro;M"]$speed_nm_s)

(
gg_mot <-
ggplot()+
  geom_errorbar(data = data_by_id, aes(id,
                                       y = speed_nm_s,
                                       ymin = speed_nm_s - sd_speed,
                                       ymax = speed_nm_s + sd_speed,
                                       color = id),
                width = 0.4,
                size = 1
                )+
  ## geom_col(data = data_by_id, aes(id, speed_nm_s, fill = id))+
  geom_hpline(data = data_by_id, aes(id, speed_nm_s, color = id), width = 0.3)+
  geom_quasirandom(data = data_by_video, aes(id, `speed_nm_s`, color = id),size = 1, width = 0.25)+
  ## annotate("richtext",
  ##          x = 2,
  ##          y = 25,
  ##          label = paste0("&darr;&ensp;", as.character(100-percent_drop), "%"),
  ##          fill = NA,
  ##          label.color = NA )+
  draw_line(c(1, 2), c(370, 370))+
  annotate("text", x = 1.5, y = 400, label = paste0("p < 0.001 "), size = 3)+
  ## geom_jitter(data = data_by_video, aes(id, `Mean v [nm/sec]`), width = 0.1, size = 0.2)+
  coord_cartesian(ylim = c(0, NA))+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.1)))+
  ylab("Motility Speed (nm/s)")+
  xlab("[danicamtiv]")+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
   axis.text.x = element_markdown(),
   axis.title.y = element_markdown()
   )
  )

saveRDS(gg_mot, "figures/gg_mot.rds")

png("figures/motility-and-atpase.png", width = 6.5, height = 3.5, units = "in", res = 300)
plot_grid(ggplot(), gg_mot, labels = c("A", "B"), rel_widths = c(1.5, 1))
dev.off()


pdf("figures/for-the-boss/unregulated-motility.pdf", width = 2.5, height = 3.5)
gg_mot
dev.off()



svg("figures/for-the-boss/unregulated-motility.svg", width = 2.5, height = 3.5)
gg_mot
dev.off()



## gg_speed <-
## ggplot(data_by_id)+
##   geom_point(aes(conc_num, speed_nm_s, color = id2))+
##   geom_errorbar(aes(x = conc_num,
##                     ymin = speed_nm_s - sd_speed,
##                     ymax = speed_nm_s + sd_speed,
##                     color = id2),
##                 width = 1)+
##   coord_cartesian(ylim = c(0, NA), xlim = c(0, NA))+
##   xlab("[Myosin] (&micro;g/mL)")+
##   ylab("Speed (nm/s)")+
##   scale_color_manual(values = colorz, name = "")+
##   theme_cowplot()+
##   theme(
##        legend.text = element_markdown(),
##        axis.title.x = element_markdown(),
##     legend.position = "none")

## gg_rel_speed <-
## ggplot(data_by_id)+
##   geom_point(aes(conc_num, relative_speed_nm_s, color = id2))+
##   geom_errorbar(aes(x = conc_num,
##                     ymin = relative_speed_nm_s - sd_rel_speed,
##                     ymax = relative_speed_nm_s + sd_rel_speed,
##                     color = id2),
##                 width = 1)+
##   coord_cartesian(ylim = c(0, NA), xlim = c(0, NA))+
##   xlab("[Myosin] (&micro;g/mL)")+
##   ylab("Relative Speed")+
##   scale_color_manual(values = colorz, name = "")+
##   theme_cowplot()+
##   theme(
##        legend.text = element_markdown(),
##        axis.title.x = element_markdown(),
##     legend.position = c(0.5, 0.4)
##     )

## plot_grid(gg_speed, gg_rel_speed)

## ggsave("~/Downloads/myosin-concentration-curve-dani.pdf", bg = "white")


## try <- data_by_id[id == "dmso"]
## vmax <- max(try$speed_nm_s)

## mod <- nls(speed_nm_s ~ vmax * (1 - (1-f)^conc_num), data = try, start = list(vmax = 200, f = 0.05))

## nd <- data.frame(conc_num = seq(0, 100, 0.1))
## nd$y <- predict(mod, newdata = nd)





## ggplot(try)+
##   geom_point(aes(conc_num, speed_nm_s, color = id2))+
##   geom_errorbar(aes(x = conc_num,
##                     ymin = speed_nm_s - sd_speed,
##                     ymax = speed_nm_s + sd_speed,
##                     color = id2),
##                 width = 1)+
##   geom_line(data = nd, aes(x = conc_num, y = y))+
##   coord_cartesian(ylim = c(0, NA), xlim = c(0, NA))+
##   scale_color_manual(values = colorz, name = "")+
##   theme_cowplot()+
##   theme(
##        legend.text = element_markdown(),
##        axis.title.x = element_markdown(),
##     legend.position = "none")
