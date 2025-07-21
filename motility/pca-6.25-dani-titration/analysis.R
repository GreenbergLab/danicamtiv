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

data <- read_mtrackj(folder = "motility/pca-6.25-dani-titration",
                     id_split = c("mot", "folder", "date","pca", "id", "video", "filename"),
                     type = "Tracks")

data$id <- sub("uM-dani", " &micro;M", data$id)
data$id <- factor(data$id, levels = c("0 &micro;M", "1 &micro;M", "10 &micro;M", "100 &micro;M"))


colorz <- c("#666666", alpha("#e7298a", 0.33), alpha("#e7298a", 0.66), "#e7298a")

(
gg_dani_titration_pca625 <-
ggplot()+
  geom_quasirandom(data = data, aes(id, `Mean v [nm/sec]`, color = id),size = 1, width = 0.25)+
  coord_cartesian(ylim = c(-0.1, 160))+
  scale_y_continuous(expand = expansion(c(0.02, 0.1), c(0, 0.1)))+
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
)


saveRDS(gg_dani_titration_pca625, "figures/gg_pca625_dani_titration.rds")

## ggsave("figures/pCa-6.25-dani-titration.png", bg = "white")

##STATS

ktest <- kruskal.test(x = `Mean v [nm/sec]`, g = id, data = data)

xx <- dunn.test::dunn.test(data$`Mean v [nm/sec]`, data$id, method = "none")
