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

data <- read_mtrackj(folder = "motility/pCa-6.25-drug-titrations",
                     id_split = c("mot", "folder", "date","pca", "id", "video", "filename"),
                     type = "Tracks")

data$id <- sub("uM-dani", " &micro;M Dani", data$id)
data$id <- sub("uM-om", " &micro;M OM", data$id)
data$id <- factor(data$id, levels = c("0 &micro;M Dani", "1 &micro;M Dani", "10 &micro;M Dani", "100 &micro;M Dani",
                                      "0 &micro;M OM", "1 &micro;M OM", "10 &micro;M OM", "100 &micro;M OM"))


data[, drug := ifelse(grepl("Dani", id, fixed = TRUE), "Danicamtiv", "OM")]

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

data_by_id[, drug := ifelse(grepl("Dani", id, fixed = TRUE), "Danicamtiv", "OM")]


colorz <- c("#666666", alpha("#e7298a", 0.33), alpha("#e7298a", 0.66), "#e7298a",
            "#666666", alpha("#d95f02", 0.33), alpha("#d95f02", 0.66), "#d95f02")

## percent_drop <- round((data_by_id[id == "10 &micro;M"]$speed_nm_s/data_by_id[id == "0 &micro;M"]$speed_nm_s)*100,
##                       0)


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
  facet_wrap(~drug, scales = "free_x")+
  ## annotate("richtext",
  ##          x = 2,
  ##          y = 25,
  ##          label = paste0("&darr;&ensp;", as.character(100-percent_drop), "%"),
  ##          fill = NA,
  ##          label.color = NA )+
  ## draw_line(c(1, 2), c(370, 370))+
  ## annotate("text", x = 1.5, y = 400, label = paste0("p < 0.001 "), size = 3)+
  ## geom_jitter(data = data_by_video, aes(id, `Mean v [nm/sec]`), width = 0.1, size = 0.2)+
  coord_cartesian(ylim = c(-0.1, 160))+
  scale_y_continuous(expand = expansion(c(0.02, 0.1), c(0, 0.1)))+
  ylab("Motility Speed (nm/s)")+
  xlab("")+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
   axis.text.x = element_markdown(),
   axis.title.y = element_markdown()
   )
  ## )


ggsave("figures/pCa-6.25-dani-titration.png", bg = "white")
