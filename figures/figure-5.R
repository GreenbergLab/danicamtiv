library(data.table)
library(ggplot2)
library(cowplot)
library(readxl)
library(drc)
library(ggtext)
library(ggbeeswarm)

theme_set(theme_cowplot(10))


#################
# MOTILITY
#################

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

mot_data <- read_mtrackj(folder = "motility/regulated",
                     id_split = c("mot", "reg", "date", "id", "pCa", "video", "filename"),
                     type = "Tracks")

mot_data$id <- sub("uM-dani", " &micro;M Dani", mot_data$id)
mot_data$id <- sub("uM-OM", " &micro;M OM", mot_data$id)
mot_data$id <- sub("control", "0 &micro;M", mot_data$id)
mot_data$id <- factor(mot_data$id, levels = c("0 &micro;M",
                                      "1 &micro;M Dani",
                                      "10 &micro;M Dani",
                                      "100 &micro;M Dani",
                                      "1 &micro;M OM",
                                       "10 &micro;M OM"))

mot_data[, pca := as.numeric(sub("pCa-", "",pCa))]

mot_data <- mot_data[date != "2023-09-22" | pca != 4]
## data <- data[date != "2023-09-22" | pca != 4]
mot_data <- mot_data[date != "2023-10-12" | pca != 9]


## fit_hill <- function(data, response){
##   drm(paste0(response, "~ pca"),
##       data = data,
##       fct = LL.4(names = c("hillslope", "min", "vmax", "ec50")),
##       logDose = 10)
## }

fit_hill3 <- function(data, response){
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

## normalize_data <- function(data, vmax, vmin){
##  data[, .(pca = pca,
##           norm_speed = (speed_nm_s-vmin)/(vmax-vmin),
##           sd_norm_speed = (sd_speed_nm_s)/(vmax-vmin))]
## }

## normalize_predict <- function(data, vmax, vmin){
##  data[, .(pca = pca,
##           Prediction = (Prediction-vmin)/(vmax-vmin))]
## }

data_by_video <- mot_data[, .(speed_nm_s = mean(`Mean v [nm/sec]`)),
                       ## sd_speed_nm_s = sd(`Mean v [nm/sec]`),
                   by = .(id, pca, date, video)]



data_by_date <- mot_data[, .(speed_nm_s = mean(`Mean v [nm/sec]`),
                             sd_speed_nm_s = sd(`Mean v [nm/sec]`)),
                         by = .(id, pca, date)]


data_by_id <- data_by_date[, .(speed_nm_s = mean(speed_nm_s),
                               sd_speed_nm_s = sd(speed_nm_s),
                               N = .N),
                           by = .(id, pca)]

## data_by_id <- data_by_date[, .(speed_nm_s = mean(speed_nm_s2),
##                        sd_speed_nm_s = sd(speed_nm_s2)),
##                    by = .(id, pca)]

## data_by_id[, `:=`(norm_speed = speed_nm_s/max(speed_nm_s),
##                   sd_norm_speed = sd_speed_nm_s/max(speed_nm_s)), by = id]

pca_nest <- data_by_id[, .(data = list(.SD)), by = id]
pca_nest[, hill_fit := lapply(data, fit_hill3, response = "speed_nm_s")]
pca_nest[, predict_fit := lapply(hill_fit, predict_hill)]
pca_nest[, mod_table := lapply(hill_fit, broom::tidy)]
pca_nest[, vmax_val := lapply(predict_fit, get_vmax)]
## pca_nest[, vmin_val := lapply(predict_fit, get_vmin)]
## pca_nest[, norm_data := mapply(normalize_data,
##                                data = data,
##                                vmin = vmin_val,
##                                vmax = vmax_val,
##                                SIMPLIFY = FALSE)]
## pca_nest[, norm_predict := mapply(normalize_predict,
##                                   data = predict_fit,
##                                   vmin = vmin_val,
##                                   vmax = vmax_val,
##                                   SIMPLIFY = FALSE)]

#norm
## pca_nest[, hill_fit_norm := lapply(data, fit_hill, response = "norm_speed")]
## pca_nest[, predict_fit_norm := lapply(hill_fit_norm, predict_hill)]

mot_pca_lines <- pca_nest[, predict_fit[[1]], by = id]


(gg_speed_all <-
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
  geom_line(data = mot_pca_lines,
            aes(pca, Prediction, color = id),
            linewidth = 0.8)+
  scale_x_reverse()+
  ## facet_wrap(~id)+
  ylab("Speed (nm/s)")+
  xlab("pCa")+
  scale_color_manual(values = c("#666666",
                                alpha("#e7298a", 0.33),
                                alpha("#e7298a", 0.66),
                                "#e7298a",
                                alpha("#d95f02", 0.5),
                                "#d95f02"))+
  theme_cowplot(11)+
  theme(
    legend.text = element_markdown()
    ## axis.title.y = element_markdown(),
    ## legend.position = "none"
    )
  )


(gg_speed <-
ggplot()+
  geom_errorbar(data = data_by_id[id %in% c("0 &micro;M", "10 &micro;M Dani")],
                aes(pca,
                    ymin = speed_nm_s-sd_speed_nm_s,
                    ymax = speed_nm_s+sd_speed_nm_s,
                    color = id),
                width = 0.1)+
  geom_point(data = data_by_id[id %in% c("0 &micro;M", "10 &micro;M Dani")],
             aes(pca, speed_nm_s, color = id),
             size = 1.5)+
  geom_line(data = mot_pca_lines[id %in% c("0 &micro;M", "10 &micro;M Dani")],
            aes(pca, Prediction, color = id),
            linewidth = 0.8)+
  scale_x_reverse()+
  ## facet_wrap(~id)+
  ylab("Speed (nm/s)")+
  xlab("pCa")+
  scale_color_manual(values = c("#666666",
                                ## alpha("#e7298a", 0.33),
                                ## alpha("#e7298a", 0.66),
                                "#e7298a",
                                alpha("#d95f02", 0.5),
                                "#d95f02"))+
  theme(
    legend.text = element_markdown(),
    ## legend.position = "none"
    ## axis.title.y = element_markdown(),
    legend.position = "none"
    )
  )

################3
## global mod for compare
########################

mod_table <- pca_nest[, mod_table[[1]], by = id]


mod_table2 <- dcast(mod_table,  id ~ term, value.var = "estimate")

mod_table2[, pCa := log10(ec50)]


global_mod <- drm(speed_nm_s ~ pca,
                  data = data_by_id[id %in% c("0 &micro;M", "10 &micro;M Dani")],
                  curveid = id,
                  ## logDose = 10,
                  fct = LL.3(names = c("hillslope", "vmax", "ec50")),
                  )

plot(global_mod)


compParm(global_mod, strVal = "hillslope")
## compParm(global_mod, strVal = "min")
compParm(global_mod, strVal = "vmax")
compParm(global_mod, strVal = "ec50", "/")

#######################
## MOTILITY TITRATION #
#######################

titration_data <- read_mtrackj(folder = "motility/pca-6.25-dani-titration",
                     id_split = c("mot", "folder", "date","pca", "id", "video", "filename"),
                     type = "Tracks")

titration_data$id <- sub("uM-dani", " <br> &micro;M", titration_data$id)
titration_data$id <- factor(titration_data$id, levels = c("0 <br> &micro;M", "1 <br> &micro;M", "10 <br> &micro;M", "100 <br> &micro;M"))

colorz_titration <- c("#666666",
            "#f7b8d8",
            "#ef72b2",
            "#e7298a")


stat_df <- data.frame(x = c(3, 4), y = 140, label = c("*", "**"))

xx <- dunn.test::dunn.test(titration_data$`Mean v [nm/sec]`, titration_data$id, method = "none")

pstats <- xx$P[c(1, 2, 3)]

stats_df2 <- data.frame(x = c(1.5, 2, 2.5),
                        x1 = c(1, 1, 1),
                        x2 = c(2, 3, 4),
                        y = c(100, 125, 150),
                        label = c("p = 0.18", rep("p < 0.001", 2)))

(
gg_dani_titration_pca625 <-
ggplot()+
  geom_quasirandom(data = titration_data, aes(id, `Mean v [nm/sec]`, color = id),size = 1, width = 0.25)+
  geom_text(data = stats_df2, aes(x = x, y = y, label = label), size = 7/.pt)+
  geom_segment(data = stats_df2, aes(x = x1, xend = x2, y = y-12, yend = y-12))+
  annotate("text", x = 2.5, y = Inf, label = "pCa 6.25", size = 9/.pt, hjust= 0.5, vjust = 1)+
  coord_cartesian(ylim = c(-0.1, 165))+
  scale_y_continuous(expand = expansion(c(0.02, 0.1), c(0, 0.1)))+
  ylab("Speed (nm/s)")+
  xlab("[danicamtiv]")+
  scale_color_manual(values = colorz_titration)+
  scale_fill_manual(values = colorz_titration)+
  theme(
    legend.position = "none",
   axis.text.x = element_markdown(),
   axis.title.y = element_markdown()
   )
)

############################
## plot motility together ##
############################

top0 <- plot_grid(gg_speed, gg_dani_titration_pca625, labels = LETTERS)



title1 <- ggdraw() +
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

(top <- plot_grid(title1, top0, rel_heights = c(0.1, 1), ncol = 1))

######################
# FIBERSIM FORCE PCA #
######################

fit_hill <- function(data){
  drm(hs_force ~ hs_pCa,
      data = data,
      fct = LL.3(names = c("hillslope", "fmax", "ec50")),
      logDose = 10)
}

predict_hill <- function(x){
  ## browser()
  dummy_x <- expand.grid(hs_pCa=exp(seq(log(4), log(9), length=1000)))
  predict_y <- predict(x, newdata = dummy_x, interval = "confidence")
  df <- as.data.table(predict_y)
  df$hs_pCa <- dummy_x$hs_pCa
  ## df <- data.table(dummy_x,
  ## y = predict_y)
  return(df)
}

get_fmax <- function(df){
  max(df$Prediction)
}


normalize_data <- function(data, fmax){
 data[, .(hs_pCa = hs_pCa,
          norm_force = hs_force/fmax)]
}

normalize_predict <- function(data, fmax){
 data[, .(hs_pCa = hs_pCa,
          Prediction = Prediction/fmax)]
}

fs_force_data <- as.data.table(read_xlsx("fibersim/project_greenberg_fibersim/sim_data/force_pCa/sim_output/pCa_analysis.xlsx"))


fs_force_data <- fs_force_data[curve %in% c(1,6)]

fs_force_data$rep <- rep(1:5, length.out = nrow(fs_force_data))

fs_pca_nest <- fs_force_data[, .(data = list(.SD)), by = .(curve, rep)]


fs_pca_nest[, hill_fit := lapply(data, fit_hill)]
fs_pca_nest[, predict_fit := lapply(hill_fit, predict_hill)]
fs_pca_nest[, mod_table := lapply(hill_fit, broom::tidy)]
fs_pca_nest[, fmax_val := lapply(predict_fit, get_fmax)]
fs_pca_nest[, norm_data := mapply(normalize_data,
                               data = data,
                               fmax = fmax_val,
                               SIMPLIFY = FALSE)]
fs_pca_nest[, norm_predict := mapply(normalize_predict,
                                  data = predict_fit,
                                  fmax = fmax_val,
                                  SIMPLIFY = FALSE)]



fs_pca_lines <- fs_pca_nest[, predict_fit[[1]], by = curve]
## pca_lines_norm <- pca_nest[, predict_fit_norm[[1]], by = id]

fs_pca_lines_norm <- fs_pca_nest[, norm_predict[[1]], by = .(curve, rep)]
fs_data_norm <- fs_pca_nest[, norm_data[[1]], by = .(curve, rep)]


fs_data_fit_avg <- fs_pca_lines_norm[, .(norm_force = mean(Prediction)), by = .(curve, hs_pCa)]
fs_data_avg <- fs_pca_lines[, .(force = mean(Prediction)), by = .(curve, hs_pCa)]

## data_norm$id <- ifelse(data_norm$curve==1, "Base Model",
##                        "2X Attachment Rate")

## data_norm$id <- factor(data_norm$id, levels = c("Base Model", "2X Attachment Rate"))

## pca_lines_norm$id <- ifelse(pca_lines_norm$curve==1, "Base Model",
##                        "2X Attachment Rate")
## pca_lines_norm$id <- factor(pca_lines_norm$id, levels = c("Base Model", "2X Attachment Rate"))

colorz1 <- c("#666666",
             ## alpha("#e7298a", 0.2),
             ## alpha("#e7298a", 0.5),
             ## alpha("#e7298a", 0.6),
             ## alpha("#e7298a", 0.8),
             "#e7298a")


## pca50_values <- pca_nest[, mod_table[[1]], by = .(curve,rep)]
## pca50_values <- pca50_values[, -4 ]
## pca50_values <- pca50_values[term == "ec50"]
## pca50_values[, pca50 := log10(estimate)]

## pca50_values[, .(mean = mean(pca50),
##                  sd = sd(pca50)), by = curve]

## pca50_values1 <- pca50_values[curve == 1]
## pca50_values2 <- pca50_values[curve == 6 ]

## t.test(pca50_values1$pca50, pca50_values2$pca50)

make_force_pca_rel <- max(fs_data_avg[curve==1]$force)

(
force_pca <-
ggplot()+
  ## geom_point(data = data_norm,
  ##            aes(hs_pCa, norm_force,
  ##                color = as.factor(curve)),
  ##            size = 0.75,
  ##            shape = 16)+
  geom_line(data = fs_pca_lines,
            aes(hs_pCa, Prediction/make_force_pca_rel,
                color = as.factor(curve)),
            linewidth = 2,
            alpha = 0.5)+
  geom_line(data = fs_data_avg, aes(hs_pCa, force/make_force_pca_rel, color = as.factor(curve)))+
  scale_x_reverse()+
  ## scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  scale_color_manual(values = colorz1)+
  ylab("Relative Force")+
  xlab("pCa")+
  theme(
    axis.title.y = element_markdown(),
    legend.position = "none")
  )

##################################
#fibersim twitch#
##################################

twitch_files  <- list.files("fibersim/project_greenberg_fibersim/sim_data/twitch/sim_output",
                            pattern = "sim_prot",
                            recursive = TRUE,
                            full.names = TRUE)

fread_twitch <- function(x){
  dat <- fread(x)
  dat$id <- x
    return(dat)
}

twitch_data <- rbindlist(lapply(twitch_files, fread_twitch))
twitch_data[, c("fibersim", "p", "sim", "t", "so", "id2", "file") := tstrsplit(id,"/")]
twitch_data[, c("sim", "prot", "num", "num2", "rep") := tstrsplit(file,"_")]

max_control <- max(twitch_data[id2 == 1]$hs_1_force)

twitch_data <- twitch_data[id2 %in% c(1, 6)]

twitch_data$id3 <- ifelse(twitch_data$id2 == 1, "0 &micro;M", "10 &micro;M")

## twitch_data[, force := hs_1_force-min(hs_1_force)]

twitch_data_integral <- twitch_data[, .(integral = sum(hs_1_force*0.001)), by = .(id3, rep) ]
twitch_data_integral_avg <- twitch_data_integral[, .(integral = mean(integral),
                                           sd_integral = sd(integral)),
                                          by = .(id3) ]

twitch_average <- twitch_data[, .(force = mean(hs_1_force)), by = .(id2, time)]
twitch_average[, force0 := force-min(force), by = id2]

make_twitch_rel <- max(twitch_average[id2==1]$force0)

twitch_average[, force0_rel := force0/make_twitch_rel]

twitch_average[, norm_twitch := (force0-min(force0))/(max(force0)-min(force0)), by = id2]

twitch_stats <- t.test(twitch_data_integral[id3 == "0 &micro;M"]$integral,
                       twitch_data_integral[id3 == "10 &micro;M"]$integral)


twitch_data[, force0 := hs_1_force-min(hs_1_force), by = .(id2, rep)]
twitch_data[, force0_rel := force0/make_twitch_rel]



peak_time <- twitch_data[, .(peak_time_index = .SD$time[which.max(.SD$hs_1_force)]), by = .(id2, rep)]

peak_time_avg <- peak_time[ , .(mean = mean(peak_time_index),
               sd = sd(peak_time_index)), by = id2]

relax_vals <- vector("list")
for(i in seq_along(unique(twitch_average$id2))){
  dat <- copy(twitch_average[id2==unique(twitch_average$id2)[[i]]])
  dat <- dat[time > time[which(force0_rel==max(force0_rel))]]
  time_stamp <- dat$time[which(dat$force0_rel<=(0.25*max(dat$force0_rel)))][[1]]
  relax_vals[[i]] <- dat[time==time_stamp]
}


relax_vals_df <- rbindlist(relax_vals)

(
ggtwitch <-
ggplot(twitch_data)+
  geom_line(aes(x = time, force0_rel, color = id2), linewidth = 1, alpha = 0.5)+
  geom_line(data = twitch_average, aes(x = time, force0_rel, color = id2), alpha = 1)+
  geom_vline(data = peak_time_avg, aes(xintercept=mean, color = id2), linetype = "dashed")+
  geom_point(data = relax_vals_df, aes(x = time, y = force0_rel, color = id2), size = 2.5 )+
  scale_color_manual(values = colorz1)+
  xlab("Time (s)")+
  ylab("Relative Force")+
  scale_y_continuous(breaks = seq(0, 1.25, by = 0.25))+
  theme(legend.position = "none"
        ## axis.text.y.left = element_blank(),
        ## axis.ticks.y.left = element_blank(),
        )
  )



bottom0 <- plot_grid(force_pca, ggtwitch, labels = c("C", "D"))

title2 <- ggdraw() +
  draw_label(
    "Computer Simulations",
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

(bottom <- plot_grid(title2, bottom0, ncol = 1, rel_heights = c(0.1, 1)))

png("figures/figure-5-new.png", width = 4.5, height = 4, units = "in", res = 300)
plot_grid(top, bottom, ncol = 1)
dev.off()

pdf("figures/figure-5-new.pdf", width = 4.5, height = 4)
plot_grid(top, bottom, ncol = 1)
dev.off()


###SUPPP MOT CURVES
png("supplemental-figures/reg-motility-all-drug-supp.png", height = 2.5, width = 4, units = "in", res = 300)
gg_speed_all
dev.off()


pdf("supplemental-figures/reg-motility-all-drug-supp.pdf", height = 2.5, width = 4)
gg_speed_all
dev.off()
