library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
library(minpack.lm)

## basesize <- 12
## fz <- basesize/.pt
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

data <- read_mtrackj(folder = "motility/myosin-concentration-curves",
                     id_split = c("mot", "myo", "date", "id", "conc", "video", "filename"),
                     type = "Tracks")

data$id <- factor(data$id, levels = c("dmso", "10uM-dani", "control"))

data[, conc_num := as.numeric(sub("ug-mL", "", conc))]




data_by_day <- data[, .(speed_nm_s = mean(`Mean v [nm/sec]`)), by = .(conc_num, date, conc, id)]
data_by_day[, relative_speed_nm_s := speed_nm_s/max(speed_nm_s), by = .(date, id) ]

data_by_id <- data_by_day[, .(
  speed_nm_s = mean(speed_nm_s),
  sd_speed = sd(speed_nm_s),
  se_speed = sd(speed_nm_s)/sqrt(.N),
  relative_speed_nm_s = mean(relative_speed_nm_s),
  sd_rel_speed = sd(relative_speed_nm_s),
  se_rel_speed = sd(relative_speed_nm_s)/sqrt(.N)),
  by = .(conc, conc_num, id)]

data_by_id_norm <- data_by_id[, .(relative_speed_nm_s = speed_nm_s/max(speed_nm_s),
                                  sd_rel_speed = sd_speed/max(speed_nm_s),
                                  conc_num = conc_num,
                                  id2 = id2
                                  ),
                                  by = id]


colorz <- c("#666666", "#e7298a")

data_by_id <- data_by_id[id != "control"]
data_by_id$id2 <- ifelse(test = data_by_id$id == "dmso", yes = "0 &micro;M Dani", no = "10 &micro;M Dani")

ggplot(data_by_id)+
  geom_point(aes(conc_num, `speed_nm_s`, color = id))+
  coord_cartesian(ylim = c(0, NA))


(gg_speed <-
ggplot(data_by_id)+
  geom_errorbar(aes(x = conc_num,
                    ymin = speed_nm_s - sd_speed,
                    ymax = speed_nm_s + sd_speed,
                    color = id2),
                width = 1.5)+
  geom_point(aes(conc_num, speed_nm_s, color = id2), size = 1.5)+
  geom_line(aes(conc_num, speed_nm_s, color = id2), linewidth = 0.8)+
  coord_cartesian(ylim = c(0, NA), xlim = c(0, NA))+
  xlab("[Myosin] (&micro;g/mL)")+
  ylab("Speed (nm/s)")+
  scale_color_manual(values = colorz, name = "")+
  ## theme_cowplot(basesize)+
  theme(
       legend.text = element_markdown(),
       axis.title.x = element_markdown(),
       axis.title.y = element_markdown(),
    legend.position = "none")
  )

(
gg_rel_speed <-
ggplot(data_by_id_norm)+
  geom_point(aes(conc_num, relative_speed_nm_s, color = id2), size = 1.5)+
  geom_line(aes(conc_num, relative_speed_nm_s, color = id2), linewidth = 0.8)+
  geom_errorbar(aes(x = conc_num,
                    ymin = relative_speed_nm_s - sd_rel_speed,
                    ymax = relative_speed_nm_s + sd_rel_speed,
                    color = id2),
                width = 1.5)+
  coord_cartesian(ylim = c(0, NA), xlim = c(0, NA))+
  xlab("[Myosin] (&micro;g/mL)")+
  ylab("Normalized Speed")+
  scale_color_manual(values = colorz, name = "")+
  ## theme_cowplot(basesize)+
  theme(
       legend.text = element_markdown(),
       axis.title.x = element_markdown(),
    legend.position = "none"
    )
  )


## saveRDS(gg_speed, "figures/gg_myo_conc_speed.rds")
## saveRDS(gg_rel_speed, "figures/gg_myo_conc_speed_norm.rds")

## plot_grid(gg_speed, gg_rel_speed)

## ggsave("~/Downloads/myosin-concentration-curve-dani.pdf", bg = "white")

fit_duty_ratio <- function(data){
  vmax <- 1
  nlsLM(relative_speed_nm_s ~ vmax * (1 - (1-f)^conc_num), data = data, start = list(f = 0.05))
}

predict_duty <- function(mod){
nd <- data.frame(conc_num = 1:100)
p <- predict(mod, newdata = nd)
nd$y <- p
return(nd)
}

get_duty <- function(mod){
  coef(mod)[[1]]*100
}


duty_data <- data_by_id[, .(data_nest = list(.SD)), by = .(id, id2)]
duty_data[, mod := lapply(data_nest, fit_duty_ratio)]
duty_data[, predict_line := lapply(mod, predict_duty)]
duty_data[, duty_ratio_percent := lapply(mod, get_duty)]

duty_predict_line <- duty_data[, predict_line[[1]], by = .(id, id2)]


gg_rel_speed_fit <-
ggplot(data_by_id)+
  geom_point(aes(conc_num, relative_speed_nm_s, color = id2))+
  geom_errorbar(aes(x = conc_num,
                    ymin = relative_speed_nm_s - sd_rel_speed,
                    ymax = relative_speed_nm_s + sd_rel_speed,
                    color = id2),
                width = 1)+
  geom_line(data = duty_predict_line, aes(x = conc_num, y = y, color = id2))+
  coord_cartesian(ylim = c(0, NA), xlim = c(0, NA))+
  xlab("[Myosin] (&micro;g/mL)")+
  ylab("Normalized Speed")+
  scale_color_manual(values = colorz, name = "")+
  theme_cowplot(11)+
  theme(
       legend.text = element_markdown(),
       axis.title.x = element_markdown(),
    legend.position = "none"
    )

png("supplemental-figures/motility-duty-ratio.png", width = 3.42, height = 3, units = "in", res = 300)
gg_rel_speed_fit
dev.off()

pdf("supplemental-figures/motility-duty-ratio.pdf", width = 3.42, height = 3)
gg_rel_speed_fit
dev.off()

#############################################
############################################
#### DRUG AT PCA 6.25
############################################
############################################

data_reg_625 <- read_mtrackj(folder = "motility/drug-concentration-curves",
                     id_split = c("mot","fold", "pca", "date", "conc", "video", "filename"),
                     type = "Tracks")

data_reg_625 <- data_reg_625[conc %in% c("0uM-dani", "10uM-dani", "100uM-dani")]

data_reg_625$conc <- factor(data_reg_625$conc, levels = c("0uM-dani", "10uM-dani", "100uM-dani"))

data_reg_625[, conc_num := as.numeric(sub("uM-dani", "", conc))]

data_reg_625[, conc2 := sub("uM-dani", " &micro;M", conc)]

data_reg_625_by_day <- data_reg_625[, .(speed_nm_s = mean(`Mean v [nm/sec]`)), by = .(conc_num, date, conc, conc2)]
data_reg_625_by_video <- data_reg_625[, .(speed_nm_s = mean(`Mean v [nm/sec]`)), by = .(conc_num, date, conc, conc2, video)]
## data_reg_625_by_day[, relative_speed_nm_s := speed_nm_s/max(speed_nm_s), by = .(date, id) ]

data_reg_625_by_id <-
  data_reg_625_by_video[, .(
  speed_nm_s = mean(speed_nm_s),
  sd_speed = sd(speed_nm_s),
  se_speed = sd(speed_nm_s)/sqrt(.N)
  ## relative_speed_nm_s = mean(relative_speed_nm_s),
  ## sd_rel_speed = sd(relative_speed_nm_s),
  ## se_rel_speed = sd(relative_speed_nm_s)/sqrt(.N)
  ),
  by = .(conc, conc2, conc_num)]





## pca_625_pval <- t.test(
##                    data_reg_625_by_video[conc_num == 0]$speed_nm_s,
##                    data_reg_625_by_video[conc_num == 10]$speed_nm_s
##   )

pca_625_aov <- aov(speed_nm_s ~ conc, data = data_reg_625_by_video)
pca_625_speed_tukey <- TukeyHSD(pca_625_aov)


(gg_pca_625 <-
ggplot(data_reg_625_by_id[conc_num %in% c(0, 10, 100)])+
  geom_errorbar(aes(x = conc2,
                    ymin = speed_nm_s - sd_speed,
                    ymax = speed_nm_s + sd_speed,
                    color = conc),
                width = 0.2)+
  geom_col(aes(conc2, speed_nm_s, fill = conc))+
  geom_jitter(data = data_reg_625_by_video[conc_num %in% c(0, 10, 100)],
              aes(conc2, speed_nm_s),
              width = 0.2)+
  ## geom_jitter(data = data_reg_625[conc_num %in% c(0, 10, 100)], aes(conc2, `Mean v [nm/sec]`), width = 0.2)+
  draw_line(c(1, 2), c(110, 110))+
  annotate("text", x = 1.5, y = 115, label = paste0("p = ", round(pca_625_speed_tukey$conc[1, 4], 6)), size = fz)+
  draw_line(c(2, 3), c(120, 120))+
  annotate("text", x = 2.5, y = 125, label = paste0("p = ", round(pca_625_speed_tukey$conc[3, 4], 2)), size = fz)+
  draw_line(c(1, 3), c(130, 130))+
  annotate("text", x = 2, y = 135, label = paste0("p = ", round(pca_625_speed_tukey$conc[2, 4], 6)), size = fz)+
  coord_cartesian(ylim = c(0, NA))+
  xlab("[Dani]")+
  ylab("Speed (nm&middot;s<sup>-1</sup>)")+
  ggtitle('pCa 6.25')+
  scale_y_continuous(breaks = seq(0, 150, by = 25),expand = expansion(c(0, 0.1), c(0, 0.1)))+
  ## scale_x_log10()+
  scale_color_manual(values = c(colorz, colorz[2]), name = "")+
  scale_fill_manual(values = c(colorz, colorz[2]), name = "")+
  theme_cowplot(basesize)+
  theme(
       legend.text = element_markdown(),
       axis.title.x = element_markdown(),
    axis.text.x = element_markdown(),
       axis.title.y = element_markdown(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
    )
  )


plot_grid(gg_rel_speed, gg_pca_625)


###################
##################
##
####

read_wrmtrck <- function(folder, id_split){
### helper function to read in motility track data
  read_mot <- function(x){
    d <- fread(x)
    d$path <- x
    d
  }
### pCa curves
  wrmtrck_paths <- list.files(folder,
                              pattern = "_wrmtrck.txt",
                              recursive = TRUE,
                              full.names = FALSE)

  wrmtrck_paths <- file.path(folder, wrmtrck_paths)

  wrmtrck_data <- rbindlist(lapply(wrmtrck_paths, read_mot))
### split the file path into new columns to get out conditions metadata
  wrmtrck_data[, (id_split) := tstrsplit(path, "/", fixed = TRUE)]
}

analyze_wrmtrck <- function(data, moving_distance_threshold = 750, group, pct = 0.2){

  ## zero_motility <- function(data){
  ##   data[, AvgSpeed := ifelse(Distance < moving_distance_threshold, 0, AvgSpeed)]
  ## }

  ## zero_motility(data)
  group_vec <- c("date", group, "video")
  mot_filter <- data[Distance > moving_distance_threshold,
                     .(speed = mean(AvgSpeed),
                       n_moving = .N),
                     by = eval((group_vec))]
                                        #calcualte average speed of everything
  mot_no_filter <- data[,
                        .(speed = mean(AvgSpeed),
                          n_moving = 0),
                        by = eval((group_vec))]

  total_filaments <- data[, .(n_total = .N), by = eval((group_vec))]
                                        # anti join filtered vs not to get groups that are not in the filtered
  missing_rows <- mot_no_filter[!mot_filter, on = eval((group_vec))]

  mot <- rbind(mot_filter, missing_rows)

  mot_by_video <- mot[total_filaments, on = eval((group_vec))]

  mot_by_video[, percent_moving := round((n_moving/n_total)*100, 2)]

  group_vec2 <- c("date", group)
  mot_by_day <- mot_by_video[, .(speed = mean(speed),
                                 speed_sd = sd(speed),
                                 n_moving = sum(n_moving),
                                 n_total = sum(n_total),
                                 percent_moving = round(mean(percent_moving), 2),
                                 percent_moving_sd = sd(percent_moving)),
                             by = eval((group_vec2))]

  mot_all <- mot_by_day[, .(speed = mean(speed),
                            speed_sd = sd(speed),
                            n_moving = sum(n_moving),
                            n_total = sum(n_total),
                            percent_moving = round(mean(percent_moving), 2),
                            percent_moving_sd = sd(percent_moving),
                            n_days = .N),
                        by = eval((group))]


top_per_video <- data[order(AvgSpeed, decreasing = TRUE), head(.SD, round(pct*.N, 0)), by = eval((group_vec))]
names(top_per_video) <- make.unique(names(top_per_video))

top_per_day_avg <- top_per_video[, .(speed = mean(AvgSpeed),
                                 n_filaments = .N),
                                 by = eval((group_vec2))]


top_avg <- top_per_day_avg[, .(speed = mean(speed),
                               speed_sd = sd(speed),
                               speed_se = sd(speed)/sqrt(.N),
                               n_filaments = sum(n_filaments),
                               n_days = .N),
                           by = eval((group))]


  mot_all <- mot_by_day[, .(speed = mean(speed),
                            speed_sd = sd(speed),
                            n_moving = sum(n_moving),
                            n_total = sum(n_total),
                            percent_moving = round(mean(percent_moving), 2),
                            percent_moving_sd = sd(percent_moving),
                            percent_moving_se = sd(percent_moving)/sqrt(.N),
                            n_days = .N),
                        by = eval((group))]




  all_no_filter_day <- data[,
                     .(speed = mean(AvgSpeed),
                       speed_sd = sd(AvgSpeed)
                       ## n_filaments = .N
                       ),
                     by = eval((group_vec2))]


  all_no_filter <- all_no_filter_day[,
                     .(speed = mean(speed),
                       speed_sd = sd(speed),
                       speed_se = sd(speed)/sqrt(.N)),
                     by = eval((group))]


  return(
    list(
      all_no_filter = list(by_group = all_no_filter,
                           by_day = all_no_filter_day),
      filtered = list(by_group = mot_all,
                      by_day = mot_by_day,
                      by_video = mot_by_video),

      top  = list(filaments_per_video = top_per_video,
                  by_day = top_per_day_avg,
                  by_group = top_avg)
    )
  )
}


######################## pCa 6.25

reg_6.25_drug_curve <- read_wrmtrck(folder = "motility/drug-concentration-curves/pCa-6.25",
                               id_split = c("mot", "folder1", "pCa", "date", "id", "video", "filename"))


reg_6.25_drug_curve$id_num <- gsub(x = reg_6.25_drug_curve$id, pattern = "uM-dani", replacement = "")
reg_6.25_drug_curve$id_num <- as.numeric(reg_6.25_drug_curve$id_num)
reg_6.25_drug_curve$id_num2 <- gsub(x = reg_6.25_drug_curve$id, pattern = "uM-dani", replacement = " &micro;M")

mot_reg_6.25_drug_curve <- analyze_wrmtrck(data = reg_6.25_drug_curve,
                       moving_distance_threshold = 500,
                       group = c("id", "id_num", "id_num2"))

## mot_reg_625_curve_sum <-

(gg_6.25 <-
   ggplot()+
   geom_point(data = mot_reg_6.25_drug_curve$all_no_filter$by_group,
              aes(id_num2, speed, color = id_num2),
              alpha = 0.75,
              size = 3)+
   geom_errorbar(data = mot_reg_6.25_drug_curve$all_no_filter$by_group,
              aes(id_num2, ymin = speed-speed_sd, ymax = speed+speed_sd, color = as.factor(id_num2)),
              width = 0.1)+
   ylab("Speed (nm/s)")+
   xlab("[Dani] (&micro;M)")+
   ## ggtitle("All Filaments Speed")+
   coord_cartesian(ylim = c(0, NA))+
   ## scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100))+
   ## scale_y_continuous(expand = expansion(c(0, 1), c(0, 0.1)))+
   ## scale_x_reverse()+
scale_color_manual(values = c(colorz[1], rep(colorz[2], 7)), name = "")+
   ## scale_fill_manual(values = colorz, name = "")+
   theme_cowplot(12)+
   theme(
     legend.position = "none",
     axis.title.x = element_markdown(),
     plot.title = element_text(size = 12, face = "plain", hjust = 0.5)
)
  )

pm_aov <- aov(percent_moving ~ id, data = mot_reg_6.25_drug_curve$filtered$by_video[id_num%in%c(0, 10, 100)])

pm_tukey <- TukeyHSD(pm_aov)
## saveRDS(mot_speed, "motility/motility-speed-plot.rds")

(mot_pm_6.25_titration <-
ggplot()+
  geom_errorbar(data = mot_reg_6.25_drug_curve$filtered$by_group[id_num %in% c(0, 10, 100)],
                aes(id_num2,
                    ymin = percent_moving-percent_moving_sd,
                    ymax = percent_moving+percent_moving_sd,
                    color = id_num2,
                width = 0.2))+
  geom_col(data = mot_reg_6.25_drug_curve$filtered$by_group[id_num %in% c(0, 10, 100)],
             aes(id_num2,
                 percent_moving,
                fill = id_num2))+
  geom_jitter(data = mot_reg_6.25_drug_curve$filtered$by_video[id_num %in% c(0, 10, 100)],
              aes(x = id_num2, y = percent_moving), width = 0.2)+
  draw_line(c(1, 2), c(100, 100))+
  annotate("text", x = 1.5, y = 105, label = paste0("p = ", round(pm_tukey$id[2, 4], 4)), size = fz)+
  draw_line(c(2, 3), c(110, 110))+
  annotate("text", x = 2.5, y = 115, label = paste0("p = ", round(pm_tukey$id[3, 4], 4)), size = fz)+
  draw_line(c(1, 3), c(120, 120))+
  annotate("text", x = 2, y = 125, label = paste0("p = ", round(pm_tukey$id[1, 4], 6)), size = fz)+
   xlab("[Dani]")+
  ylab("Percent Moving (%)")+
  ggtitle("pCa 6.25")+
   coord_cartesian(ylim = c(0, NA))+
   ## ggtitle("Filaments moving >1000 um")+
   ## scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100))+
   scale_y_continuous(breaks = seq(0, 150, by = 25), expand = expansion(c(0, 0.1), c(0, 0.1)))+
   ## scale_x_reverse()+
   scale_color_manual(values = c(colorz[1], rep(colorz[2], 7)), name = "")+
   scale_fill_manual(values = c(colorz[1], rep(colorz[2], 7)), name = "")+
   ## scale_fill_manual(values = colorz, name = "")+
   theme_cowplot(basesize)+
   theme(
     legend.position = "none",
     axis.title.x = element_markdown(),
     axis.text.x = element_markdown(),
     plot.title = element_text(hjust = 0.5)
)
  )



top <- plot_grid(gg_speed, gg_rel_speed_fit, nrow = 1, labels = c("A", "B"))
bottom <- plot_grid(gg_pca_625, mot_pm_6.25_titration, labels = c("C", "D"))


png('figures/attachment-limited-motility.png', width = 6.5, height = 5.5, res = 300, units = "in")
plot_grid(top, bottom, nrow = 2, rel_heights = c(1, 1.25))
dev.off()

pdf('figures/for-the-boss/attachment-limited-motility.pdf', width = 6.5, height = 5.5)
plot_grid(top, bottom, nrow = 2, rel_heights = c(1, 1.25))
dev.off()



svg('figures/for-the-boss/attachment-limited-motility.svg', width = 6.5, height = 5.5)
plot_grid(top, bottom, nrow = 2, rel_heights = c(1, 1.25))
dev.off()

gg_6.25_all <- plot_grid(gg_6.25, mot_pm_6.25_titration, ncol = 2)
reg6.25_title <- ggdraw() +
  draw_label(
    "pCa 6.25 Dani Titration",
    fontface = 'bold',
    size = 12,
    x = 0.5,
    hjust = 0.5
  )


gg6.25 <- plot_grid(reg6.25_title, gg_6.25_all, ncol = 1, rel_heights = c(0.1, 1))
gg6.25
