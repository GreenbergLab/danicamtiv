library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
library(drc)

theme_set(theme_cowplot(11))

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

data <- read_mtrackj(folder = "motility/reg-1uM-drugs",
                     id_split = c("mot", "reg", "date", "id", "pCa", "video", "filename"),
                     type = "Tracks")

## data$id <- sub("uM-OM", " &micro;M OM", data$id)
## data$id <- sub("control", "0 &micro;M", data$id)
## data$id <- factor(data$id, levels = c("0 &micro;M OM", "10 &micro;M OM"))

data[, pca := as.numeric(sub("pCa-", "",pCa))]

## data <- data[date != "2023-09-22" | pca != 4]


fit_hill <- function(data, response){
  drm(paste0(response, "~ pca"),
      data = data,
      fct = LL.3(names = c("hillslope", "vmax", "ec50")),
      logDose = 10)
}

predict_hill <- function(x){
  ## browser()
  dummy_x <- expand.grid(pca=exp(seq(log(4), log(9), length=1000)))
  predict_y <- predict(x, newdata = dummy_x, interval = "confidence")
  df <- as.data.table(predict_y)
  df$pca <- dummy_x$pca
  ## df <- data.table(dummy_x,
  ## y = predict_y)
  return(df)
}

get_vmax <- function(df){
  max(df$Prediction)
}

get_vmin <- function(df){
  min(df$Prediction)
  ## return(0)
}

normalize_data <- function(data, vmax, vmin){
 data[, .(pca = pca,
          norm_speed = (speed_nm_s-vmin)/(vmax-vmin),
          sd_norm_speed = (sd_speed_nm_s)/(vmax-vmin))]
}

normalize_predict <- function(data, vmax, vmin){
 data[, .(pca = pca,
          Prediction = (Prediction-vmin)/(vmax-vmin))]
}

data_by_video <- data[, .(speed_nm_s = mean(`Mean v [nm/sec]`)),
                       ## sd_speed_nm_s = sd(`Mean v [nm/sec]`),
                   by = .(id, pca, date, video)]



## data_by_video[, speed_nm_s := (speed_nm_s - min(speed_nm_s)),
##                        ## sd_speed_nm_s = sd(`Mean v [nm/sec]`)),
##                    by = .(id, date, video)]

## data_by_video <- data_by_video[video != 3 & date != "2023-09-22"]

## data_by_date <- data[, .(speed_nm_s = mean(`Mean v [nm/sec]`),
##                                   sd_speed_nm_s = sd(`Mean v [nm/sec]`)),
##                               by = .(id, pca, date)]

## data_by_date[, speed_nm_s2 := speed_nm_s-min(speed_nm_s), by = .(id, date)]

## data_by_id <- data_by_video[, .(speed_nm_s = mean(speed_nm_s),
##                        sd_speed_nm_s = sd(speed_nm_s)),
##                    by = .(id, pca)]

## data_by_id <- data[, .(speed_nm_s = mean(`Mean v [nm/sec]`),
##                        sd_speed_nm_s = sd(`Mean v [nm/sec]`)),
##                    by = .(id, pca)]


data_by_date <- data[, .(speed_nm_s = mean(`Mean v [nm/sec]`),
                       sd_speed_nm_s = sd(`Mean v [nm/sec]`)),
                   by = .(id, pca, date)]


data_by_id <- data_by_date[, .(speed_nm_s = mean(speed_nm_s),
                               sd_speed_nm_s = sd(speed_nm_s),
                               N = .N),
                           by = .(id, pca)]

## data_by_id <- data_by_id[id == "10uM-OM"]

## data_by_id <- data_by_date[, .(speed_nm_s = mean(speed_nm_s2),
##                        sd_speed_nm_s = sd(speed_nm_s2)),
##                    by = .(id, pca)]

## data_by_id[, `:=`(norm_speed = speed_nm_s/max(speed_nm_s),
##                   sd_norm_speed = sd_speed_nm_s/max(speed_nm_s)), by = id]

pca_nest <- data_by_id[, .(data = list(.SD)), by = id]
pca_nest[, hill_fit := lapply(data, fit_hill, response = "speed_nm_s")]
pca_nest[, predict_fit := lapply(hill_fit, predict_hill)]
pca_nest[, mod_table := lapply(hill_fit, broom::tidy)]
pca_nest[, vmax_val := lapply(predict_fit, get_vmax)]
pca_nest[, vmin_val := lapply(predict_fit, get_vmin)]
pca_nest[, norm_data := mapply(normalize_data,
                               data = data,
                               vmin = vmin_val,
                               vmax = vmax_val,
                               SIMPLIFY = FALSE)]
pca_nest[, norm_predict := mapply(normalize_predict,
                                  data = predict_fit,
                                  vmin = vmin_val,
                                  vmax = vmax_val,
                                  SIMPLIFY = FALSE)]

#norm
## pca_nest[, hill_fit_norm := lapply(data, fit_hill, response = "norm_speed")]
## pca_nest[, predict_fit_norm := lapply(hill_fit_norm, predict_hill)]

pca_lines <- pca_nest[, predict_fit[[1]], by = id]
## pca_lines_norm <- pca_nest[, predict_fit_norm[[1]], by = id]
pca_lines_norm <- pca_nest[, norm_predict[[1]], by = id]
data_norm <- pca_nest[, norm_data[[1]], by = id]

## fwrite(data_by_id, "motility/reg-1uM-drugs/data-by-id-1uM.csv", sep = ",")
## fwrite(data_norm, "motility/reg-1uM-drugs/data-norm-1uM.csv", sep = ",")
## fwrite(pca_lines, "motility/reg-1uM-drugs/pca-lines-1uM.csv", sep = ",")
## fwrite(pca_lines_norm, "motility/reg-1uM-drugs/pca-lines-norm-1uM.csv", sep = ",")

(gg_speed <-
ggplot()+
  geom_errorbar(data = data_by_id,
                aes(pca,
                    ymin = speed_nm_s-sd_speed_nm_s,
                    ymax = speed_nm_s+sd_speed_nm_s,
                    color = id),
                width = 0.1)+
  geom_point(data = data_by_id,
             aes(pca, speed_nm_s, color = id),
             size = 1.5)+
  geom_line(data = pca_lines,
            aes(pca, Prediction, color = id),
            linewidth = 0.8,
            alpha = 1)+
  scale_x_reverse()+
  ylab("Speed (nm/s)")+
  xlab("pCa")+
  ## scale_color_manual(values = c("#d95f02"))+
  theme(
    ## axis.title.y = element_markdown(),
    legend.position = "none"
    )
  )




## ggplot()+
##   ## geom_errorbar(data = data_by_id[id == "0 &micro;M"],
##   ##               aes(pca,
##   ##                   ymin = speed_nm_s-sd_speed_nm_s,
##   ##                   ymax = speed_nm_s+sd_speed_nm_s,
##   ##                   color = id),
##   ##               width = 0.1)+
##   ## geom_point(data = data_by_id[id == "0 &micro;M"],
##   ##            aes(pca, speed_nm_s, color = id),
##   ##            size = 1.5)+
##   geom_line(data = pca_lines[id == "0 &micro;M"],
##             aes(pca, Prediction, color = id),
##             linewidth = 1,
##             alpha = 1)+
##   scale_x_reverse()+
##   ylab("Speed (nm/s)")+
##   xlab("pCa")+
##   scale_color_manual(values = c("#666666", "#e7298a"))+
##   theme(
##     ## axis.title.y = element_markdown(),
##     legend.position = "none"
##     )

## ggsave("~/Downloads/dani-motility-pca-curves.pdf", device = cairo_pdf)

(gg_norm_speed <-
ggplot()+
  geom_errorbar(data = data_norm,
                aes(pca,
                    ymin = norm_speed-sd_norm_speed,
                    ymax = norm_speed+sd_norm_speed,
                    color = id),
                width = 0.1)+
  geom_point(data = data_norm,
             aes(pca, norm_speed, color = id),
             size = 1.5)+
  geom_line(data = pca_lines_norm,
            aes(pca, Prediction, color = id),
            linewidth = 0.8,
            alpha = 1)+
  scale_x_reverse()+
  ## scale_y_continuous(breaks = seq(0, 1.25, by = 0.25))+
  ylab("Normalized Speed")+
  xlab("pCa")+
  scale_color_manual(values = c("#d95f02"))+
  theme(
    ## axis.title.y = element_markdown(),
    legend.position = "none"
    )
  )

## saveRDS(gg_speed, "figures/gg_reg_speed.rds")
## saveRDS(gg_norm_speed, "figures/gg_reg_norm_speed.rds")

## ggsave("~/Downloads/dani-motility-pca-curves-normalized-min-max.pdf", device = cairo_pdf)


mod_table <- pca_nest[, mod_table[[1]], by = id]


mod_table2 <- dcast(mod_table,  id ~ term, value.var = "estimate")

mod_table2[, pCa := log10(ec50)]


global_mod <- drm(speed_nm_s ~ pca,
                  data = data_by_id,
                  curveid = id,
                  ## logDose = 10
                  fct = LL.4(names = c("hillslope", "min", "vmax", "ec50")),
                  )

plot(global_mod)


compParm(global_mod, strVal = "hillslope")
compParm(global_mod, strVal = "min")
compParm(global_mod, strVal = "vmax")
compParm(global_mod, strVal = "ec50", "-")


con7 <- data[id == "0 &micro;M" & pca == 7]
dani7 <- data[id == "10 &micro;M" & pca == 7]

t.test(con7$`Mean v [nm/sec]`, dani7$`Mean v [nm/sec]`)

con7[, .(mean = mean(`Mean v [nm/sec]`),
     sd = sd(`Mean v [nm/sec]`))]


dani7[, .(mean = mean(`Mean v [nm/sec]`),
     sd = sd(`Mean v [nm/sec]`))]
