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

data <- read_mtrackj(folder = "motility/unreg-dani-drug-curve",
                     id_split = c("mot", "folder", "date","unreg", "id", "video", "filename"),
                     type = "Tracks")

data$id <- sub("uM-dani", " &micro;M", data$id)
data$id <- factor(data$id, levels = c("0 &micro;M", "1 &micro;M", "10 &micro;M", "100 &micro;M"))

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

data_by_id[, rel_con := speed_nm_s/max(speed_nm_s)]


colorz <- c("#666666", alpha("#e7298a", 0.33), alpha("#e7298a", 0.66), "#e7298a")

## percent_drop <- round((data_by_id[id == "10 &micro;M"]$speed_nm_s/data_by_id[id == "0 &micro;M"]$speed_nm_s)*100,
##                       0)

gg_unreg_dani_curve <-
ggplot()+
  ## geom_errorbar(data = data_by_id, aes(id,
  ##                                      y = speed_nm_s,
  ##                                      ymin = speed_nm_s - sd_speed,
  ##                                      ymax = speed_nm_s + sd_speed,
  ##                                      color = id),
  ##               width = 0.4,
  ##               size = 1
  ##               )+
  ## geom_col(data = data_by_id, aes(id, speed_nm_s, fill = id))+
  ## geom_hpline(data = data_by_id, aes(id, speed_nm_s, color = id), width = 0.3)+
  ## geom_quasirandom(data = data_by_video, aes(id, `speed_nm_s`, color = id),size = 1, width = 0.25)+
  geom_quasirandom(data = data, aes(id, `Mean v [nm/sec]`, color = id),size = 1, width = 0.25)+
  ## annotate("richtext",
  ##          x = 2,
  ##          y = 25,
  ##          label = paste0("&darr;&ensp;", as.character(100-percent_drop), "%"),
  ##          fill = NA,
  ##          label.color = NA )+
  ## draw_line(c(1, 2), c(370, 370))+
  ## annotate("text", x = 1.5, y = 400, label = paste0("p < 0.001 "), size = 3)+
  ## geom_jitter(data = data_by_video, aes(id, `Mean v [nm/sec]`), width = 0.1, size = 0.2)+
  coord_cartesian(ylim = c(-0.1, NA))+
  ## scale_y_continuous(expand = expansion(c(0.02, 0.1), c(0, 0.1)))+
  ylab("Speed (nm/s)")+
  xlab("[danicamtiv]")+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
   axis.text.x = element_markdown(),
   axis.title.y = element_markdown()
   )
  ## )


saveRDS(gg_unreg_dani_curve, "figures/gg_mot_unreg_dani_titration.rds")

ggsave("supplemental-figures/unregulated-dani-motility-titration.png", bg = "white")


data$speed <- data$`Mean v [nm/sec]`

aov_mod <- aov(speed~id, data=data)
TukeyHSD(aov_mod)

krusk_mod <- kruskal.test(speed~id, data=data)
dunn.test::dunn.test(data$speed, data$id, method = "bonferroni")

ggplot(data = data, aes(sample = speed, color = id))+
  stat_qq()+
  stat_qq_line()+
  theme( legend.text = element_markdown())



ggsave("supplemental-figures/unregulated-dani-motility-titration-qq.png", bg = "white")
