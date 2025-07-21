########
# tension pca
#########

library(data.table)
library(ggplot2)
library(cowplot)
library(readxl)
library(drc)
library(ggtext)
library(ggbeeswarm)

theme_set(theme_cowplot(11))

## myosim_files <- list.files("~/MATMyoSim/danicamtiv/tension_pCa/test2",
##                            pattern = "pca-curve.csv",
##                            recursive = TRUE,
##                            full.names = TRUE)


## fit_files <- list.files("~/MATMyoSim/danicamtiv/tension_pCa/test2",
##                            pattern = "curve-fit.csv",
##                            recursive = TRUE,
##                            full.names = TRUE)

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

data <- as.data.table(read_xlsx("fibersim/project_greenberg_fibersim/sim_data/force_pCa/sim_output/pCa_analysis.xlsx"))

data <- data[curve %in% c(1,6)]

data$rep <- rep(1:5, length.out = nrow(data))

pca_nest <- data[, .(data = list(.SD)), by = .(curve, rep)]


pca_nest[, hill_fit := lapply(data, fit_hill)]
pca_nest[, predict_fit := lapply(hill_fit, predict_hill)]
pca_nest[, mod_table := lapply(hill_fit, broom::tidy)]
pca_nest[, fmax_val := lapply(predict_fit, get_fmax)]
pca_nest[, norm_data := mapply(normalize_data,
                               data = data,
                               fmax = fmax_val,
                               SIMPLIFY = FALSE)]
pca_nest[, norm_predict := mapply(normalize_predict,
                                  data = predict_fit,
                                  fmax = fmax_val,
                                  SIMPLIFY = FALSE)]



pca_lines <- pca_nest[, predict_fit[[1]], by = curve]
## pca_lines_norm <- pca_nest[, predict_fit_norm[[1]], by = id]

pca_lines_norm <- pca_nest[, norm_predict[[1]], by = .(curve, rep)]
data_norm <- pca_nest[, norm_data[[1]], by = .(curve, rep)]


data_fit_avg <- pca_lines_norm[, .(norm_force = mean(Prediction)), by = .(curve, hs_pCa)]
data_avg <- pca_lines[, .(force = mean(Prediction)), by = .(curve, hs_pCa)]

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


pca50_values <- pca_nest[, mod_table[[1]], by = .(curve,rep)]
pca50_values <- pca50_values[, -4 ]
pca50_values <- pca50_values[term == "ec50"]
pca50_values[, pca50 := log10(estimate)]

pca50_values[, .(mean = mean(pca50),
                 sd = sd(pca50)), by = curve]

pca50_values1 <- pca50_values[curve == 1]
pca50_values2 <- pca50_values[curve == 6 ]

t.test(pca50_values1$pca50, pca50_values2$pca50)


(
force_pca <-
ggplot()+
  ## geom_point(data = data_norm,
  ##            aes(hs_pCa, norm_force,
  ##                color = as.factor(curve)),
  ##            size = 0.75,
  ##            shape = 16)+
  geom_line(data = pca_lines,
            aes(hs_pCa, Prediction,
                color = as.factor(curve)),
            linewidth = 2,
            alpha = 0.5)+
  geom_line(data = data_avg, aes(hs_pCa, force, color = as.factor(curve)))+
  scale_x_reverse()+
  ## scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  scale_color_manual(values = colorz1)+
  ylab("Sim. Force (N&middot;m<sup>-2</sup>)")+
  xlab("pCa")+
  theme(
    axis.title.y = element_markdown(),
    legend.position = "none")
  )

## saveRDS(force_pca, "figures/gg_fibersim_force_pca.rds")

(
force_pca_norm <-
ggplot()+
  ## geom_point(data = data_norm,
  ##            aes(hs_pCa, norm_force,
  ##                color = as.factor(curve)),
  ##            size = 0.75,
  ##            shape = 16)+
  geom_line(data = pca_lines_norm,
            aes(hs_pCa, Prediction,
                color = as.factor(curve)),
            linewidth = 2,
            alpha = 0.5)+
  geom_line(data = data_fit_avg, aes(hs_pCa, norm_force, color = as.factor(curve)))+
  scale_x_reverse()+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  scale_color_manual(values = colorz1)+
  ylab("Normalized Force")+
  xlab("pCa")+
  theme(
    legend.position = "none")
  )


fig7_top <- plot_grid(force_pca, force_pca_norm, labels = LETTERS, nrow = 1)


## saveRDS(force_pca, "figures/fibersim-force-pca.rds")

################################3
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


max(twitch_average[id2==1]$force)*0.25

## relax_dt <- twitch_average[time > time[which(force==max(force))],
##                         .(relax = time[which(force <= (0.25*max(force)))][[1]],
##                         force = force[which(force<=0.25*max(force))][[1]]),
##                         by = id2]



relax_vals <- vector("list")
for(i in seq_along(unique(twitch_average$id2))){
  dat <- copy(twitch_average[id2==unique(twitch_average$id2)[[i]]])
  dat <- dat[time > time[which(force0_rel==max(force0_rel))]]
  time_stamp <- dat$time[which(dat$force0_rel<=(0.25*max(dat$force0_rel)))][[1]]
  relax_vals[[i]] <- dat[time==time_stamp]
}


relax_vals_df <- rbindlist(relax_vals)






(
twitch_integral <-
ggplot()+
  ## geom_errorbar(data = twitch_data_integral_avg,
  ##               aes(id3,
  ##                   ymin = (integral-sd_integral)/7144.99,
  ##                   ymax = (integral+sd_integral)/7144.99,
  ##                   color = id3),
  ##                   width = 0.25)+
  geom_quasirandom(data = twitch_data_integral, aes(id3, integral/7144.99, color = id3),
                   width = 0.2,
                   shape = 16,
                   size = 1.5)+
  draw_line(x = c(1, 2), y = c(3, 3))+
  annotate("text", x = 1.5, y = 3.2, label = "p < 0.001", size = 3)+
  scale_color_manual(values = colorz1)+
  coord_cartesian(ylim = c(0, 3.4))+
  ylab("Relative Integral")+
  xlab("[danicamtiv]")+
  theme(axis.text.x = element_markdown(),
        legend.position = "none")
  )
#39481.58

png("supplemental-figures/fibersim-twitch-integrals.png", height = 2.5, width = 3, units = "in", res = 300)
twitch_integral
dev.off()

plot_grid(force_pca, ggtwitch)

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

(
ggtwitch_norm <-
ggplot(twitch_data)+
  ## geom_line(aes(x = time, norm_twitch, color = id2), linewidth = 1, alpha = 0.5)+
  geom_line(data = twitch_average, aes(x = time, norm_twitch, color = id2), alpha = 1, linewidth = 0.9)+
  scale_color_manual(values = colorz1)+
  xlab("Time (s)")+
  ylab("Normalized Force")+
  coord_cartesian(xlim = c(0.1, 0.3))+
  scale_y_continuous(breaks = seq(0, 1.25, by = 0.25))+
  theme(legend.position = "none"
        ## axis.text.y.left = element_blank(),
        ## axis.ticks.y.left = element_blank(),
        )
)

fig7_bottom <-
plot_grid(ggtwitch, ggtwitch_norm, twitch_integral, nrow = 1, rel_widths = c(1, 1, 0.5), labels = c("C", "D", "E"))


## png("figures/figure-7.png", height = 4.5, width = 6.5, units = "in", res = 300)
## plot_grid(fig7_top, fig7_bottom, nrow = 2)
## dev.off()


ggtwitch1 <-
ggtwitch+
  annotate("rect",
           ymin = 0,
           ymax = 1.25,
           xmin = 0.1,
           xmax = 0.2,
           color = "black",
           fill = NA,
           linetype = "dashed")

ggtwitch2 <-
ggtwitch+coord_cartesian(xlim = c(0.1, 0.2), ylim = c(0, 1.25))+
scale_x_continuous(breaks = c(0.1, 0.15, 0.2))+
  theme(
        axis.line = element_blank(),
        panel.border = element_rect(linetype = "dashed", color = "black")
  )




## plot_grid(ggtwitch1, ggtwitch2)

## saveRDS(ggtwitch, "figures/fibersim-ggtwitch.rds")
## saveRDS(ggtwitch1, "figures/fibersim-ggtwitch1.rds")
## saveRDS(ggtwitch2, "figures/fibersim-ggtwitch2.rds")
## saveRDS(twitch_integral, "figures/fibersim-twitch-integral.rds")


##################
#############################3
###########################
####
#### SUPPLEMENTAL
#####################3
########################

data2 <- as.data.table(read_xlsx("fibersim/project_greenberg_fibersim/sim_data/force_pCa/sim_output/pCa_analysis.xlsx"))


data2$rep <- rep(1:5, length.out = nrow(data2))
pca_nest2 <- data2[, .(data = list(.SD)), by = .(curve,rep)]


pca_nest2[, hill_fit := lapply(data, fit_hill)]
pca_nest2[, predict_fit := lapply(hill_fit, predict_hill)]
pca_nest2[, mod_table := lapply(hill_fit, broom::tidy)]
pca_nest2[, fmax_val := lapply(predict_fit, get_fmax)]
pca_nest2[, norm_data := mapply(normalize_data,
                               data = data,
                               fmax = fmax_val,
                               SIMPLIFY = FALSE)]
pca_nest2[, norm_predict := mapply(normalize_predict,
                                  data = predict_fit,
                                  fmax = fmax_val,
                                  SIMPLIFY = FALSE)]



pca_lines2 <- pca_nest2[, predict_fit[[1]], by = .(curve, rep)]
## pca_lines_norm <- pca_nest[, predict_fit_norm[[1]], by = id]
pca_lines_norm2 <- pca_nest2[, norm_predict[[1]], by = .(curve, rep)]
data_norm2 <- pca_nest2[, norm_data[[1]], by = .(curve, rep)]


data_fit_avg2 <- pca_lines2[, .(force = mean(Prediction)), by = .(curve, hs_pCa)]


dark2 <- RColorBrewer::brewer.pal(8, "Dark2")
## colorz2 <- c("#666666", dark2[1:3], dark2[5], dark2[4])
colorz2 <- rainbow(6)

make_force_pca_rel <- max(data_fit_avg2[curve == 1]$force)

(
gg_force_pca_params <-
ggplot()+
  ## geom_point(data = data2,
  ##            aes(hs_pCa, hs_force,
  ##                color = as.factor(curve)),
  ##            size = 1.5,
  ##            shape = 16)+
  geom_line(data = pca_lines2,
            aes(hs_pCa,
                Prediction/make_force_pca_rel,
                color = as.factor(curve)),
            linewidth = 2,
            alpha = 0.5)+
  geom_line(data = data_fit_avg2, aes(hs_pCa,
                                      force/make_force_pca_rel,
                                      color = as.factor(curve)))+
  scale_x_reverse()+
  scale_y_continuous(breaks = seq(0, 1.25, by = 0.25))+
  coord_cartesian(ylim = c(0, 1.25))+
  scale_color_manual(name = "Simulation", values = colorz2)+
  ylab("Relative Force")+
  xlab("pCa")+
  theme(legend.position = c(0.1, 0.5)
        )
)

## saveRDS(force_pca, "figures/fibersim-force-pca.rds")

pca50_values2 <- pca_nest2[, mod_table[[1]], by = .(curve,rep)]
pca50_values2 <- pca50_values2[, -4 ]
pca50_values2 <- pca50_values2[term == "ec50"]
pca50_values2[, pca50 := log10(estimate)]

pca50_values2_sum <- pca50_values2[, .(mean = mean(pca50),
                 sd = sd(pca50)), by = curve]

gg_pca50_all_sims <-
ggplot()+
  geom_quasirandom(data = pca50_values2, aes(x = as.factor(curve), y = pca50, color = as.factor(curve)),size = 1, width = 0.25)+
  geom_text(data = pca50_values2_sum,
            aes(x = as.factor(curve),
                y = mean+0.065,
                label = round(mean, 2)),
            size = 3)+
  ## coord_cartesian(ylim = c(-0.1, NA))+
  scale_y_continuous(expand = expansion(c(0.1, 0.1), c(0, 0.1)))+
  ylab("pCa50")+
  xlab("Simulation")+
  scale_color_manual(values = colorz2)+
  ## scale_fill_manual(values = colorz)+
  theme_cowplot(11)+
  theme(
    legend.position = "none",
   axis.text.x = element_markdown(),
   axis.title.y = element_markdown()
   )

sim_force_pca <- plot_grid(gg_force_pca_params ,gg_pca50_all_sims, labels = c("A", "B"))

## ggsave("supplemental-figures/all-simulation-pca-curve-pca50-values.png", bg = "white")

################################3
##################################

twitch_files2  <- list.files("fibersim/project_greenberg_fibersim/sim_data/twitch/sim_output",
                            pattern = "sim_prot",
                            recursive = TRUE,
                            full.names = TRUE)

fread_twitch <- function(x){
  dat <- fread(x)
  dat$id <- x
    return(dat)
}

twitch_data2 <- rbindlist(lapply(twitch_files2, fread_twitch))
twitch_data2[, c("fibersim", "p", "sim", "t", "so", "id2", "file") := tstrsplit(id,"/")]

twitch_data2$rep <- rep(1:5, length.out = nrow(twitch_data2))

twitch_average2 <- twitch_data2[, .(hs_1_force = mean(hs_1_force)), by = .(id2, time)]

twitch_average2[, norm_twitch := (hs_1_force-min(hs_1_force))/(max(hs_1_force)-min(hs_1_force)), by = id2]
## max_control <- max(twitch_data[id2 == 1]$hs_1_force)

## twitch_data <- twitch_data[id2 %in% c(1, 6)]

twitch_data_integral2 <- twitch_data2[, .(integral = sum(hs_1_force*0.001)), by = .(id2, rep) ]

con_to_remove <- twitch_data_integral2[id2 == 1]$integral

peak_to_remove <- max(twitch_average2[id2==1]$hs_1_force)


(
ggintegral2 <-
ggplot(twitch_data_integral2)+
  geom_quasirandom(aes(id2, integral/con_to_remove, color = id2))+
  scale_color_manual(values = colorz2)+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.01)))+
  coord_cartesian(ylim = c(0, 5.25))+
  ylab("Normalized Integral")+
  xlab("Simulation")+
  theme(legend.position = "none")
  )

(
gg_twitch_params <-
ggplot(twitch_data2)+
  geom_line(aes(x = time, hs_1_force/peak_to_remove, color = id2), linewidth = 2, alpha = 0.5 )+
  geom_line(data = twitch_average2, aes(x = time, hs_1_force/peak_to_remove, color = id2), alpha = 1)+
  ## scale_color_manual(values = colorz1)+
  xlab("Time (s)")+
  ylab("Relative Force")+
  scale_color_manual(name = "Simulation", values = colorz2)+
  theme(legend.position = "none")
)



(
gg_twitch_norm <-
ggplot(twitch_data2)+
  ## geom_line(aes(x = time, hs_1_force, color = id2), linewidth = 2, alpha = 0.5 )+
  geom_line(data = twitch_average2,
            aes(x = time,
                norm_twitch,
                color = id2),
            alpha = 1,
            linewidth = 1.25)+
  ## scale_color_manual(values = colorz1)+
  xlab("Time (s)")+
  ylab("Normalized Force")+
  scale_color_manual(name = "Simulation", values = colorz2)+
  theme(
    legend.position = "none"
  )
)

gg_twitch_norm2 <-
  gg_twitch_norm+coord_cartesian(xlim = c(0.1, 0.2))+theme(legend.position = "right")






## png(filename = "supplemental-figures/supplemental-fibersim-5.png", width = 6.5, height = 2.5, res = 500, units = "in")
plot_grid(sim_force_pca,
          plot_grid(gg_twitch_params, ggintegral2, nrow = 1, labels = c("C", "D")),
          ## rel_widths = c(1, 1, 0.5),
          nrow = 2)
## dev.off()


## png(filename = "supplemental-figures/supplemental-fibersim-6.png", width = 6.5, height = 2.5, res = 500, units = "in")
## plot_grid(gg_twitch_norm, gg_twitch_norm2, nrow = 1, labels = LETTERS)
## dev.off()

peak_time2 <- twitch_data2[, .(peak_time_index = .SD$time[which.max(.SD$hs_1_force)]), by = .(id2, rep)]

peak_time2_sum <- peak_time2[ , .(mean = mean(peak_time_index),
               sd = sd(peak_time_index)), by = id2]


(
ggpeaktime <-
ggplot()+
  geom_quasirandom(data = peak_time2,
                 aes(id2,
                     peak_time_index,
                     color = id2))+
  scale_color_manual(values = colorz2)+
  ## scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.01)))+
  ## coord_cartesian(ylim = c(0, 5.25))+
  ylab("Time to peak (s)")+
  xlab("Simulation")+
  theme(legend.position = "none")
  )


gg_all_twitch <- plot_grid(gg_twitch_params, ggpeaktime, nrow = 1, labels = c("C", "D"))

png("supplemental-figures/fibersim-simulations-single-parameters.png", width = 6.5, height = 4.5, units = "in", res = 300)
plot_grid(sim_force_pca, gg_all_twitch, nrow = 2 )
dev.off()


pdf("supplemental-figures/fibersim-simulations-single-parameters.pdf", width = 6.5, height = 4.5)
plot_grid(sim_force_pca, gg_all_twitch, nrow = 2 )
dev.off()


## ggsave("supplemental-figures/all-simulations-twitches-integral-ttp.png", bg = "white")


##############################
###### KOOIKER/REGNIER SIMS###
##############################



data3 <- as.data.table(read_xlsx("fibersim/project_greenberg_fibersim_kooiker/sim_data/force_pCa/sim_output/pCa_analysis.xlsx"))

data3 <- data3[curve %in% c(1,2)]

data3$rep <- rep(1:5, length.out = nrow(data3))

pca_nest3 <- data3[, .(data = list(.SD)), by = .(curve, rep)]


pca_nest3[, hill_fit := lapply(data, fit_hill)]
pca_nest3[, predict_fit := lapply(hill_fit, predict_hill)]
pca_nest3[, mod_table := lapply(hill_fit, broom::tidy)]
pca_nest3[, fmax_val := lapply(predict_fit, get_fmax)]
pca_nest3[, norm_data := mapply(normalize_data,
                               data = data,
                               fmax = fmax_val,
                               SIMPLIFY = FALSE)]
pca_nest3[, norm_predict := mapply(normalize_predict,
                                  data = predict_fit,
                                  fmax = fmax_val,
                                  SIMPLIFY = FALSE)]



pca_lines3 <- pca_nest3[, predict_fit[[1]], by = .(curve, rep)]
## pca_lines_norm <- pca_nest[, predict_fit_norm[[1]], by = id]
pca_lines_norm3 <- pca_nest3[, norm_predict[[1]], by = .(curve, rep)]
data_norm3 <- pca_nest3[, norm_data[[1]], by = .(curve, rep)]


data_fit_avg3 <- pca_lines_norm3[, .(norm_force = mean(Prediction)), by = .(curve, hs_pCa)]
data_predict_avg3 <- pca_lines3[, .(force = mean(Prediction)), by = .(curve, hs_pCa)]
## data_norm$id <- ifelse(data_norm$curve==1, "Base Model",
##                        "2X Attachment Rate")

## data_norm$id <- factor(data_norm$id, levels = c("Base Model", "2X Attachment Rate"))

## pca_lines_norm$id <- ifelse(pca_lines_norm$curve==1, "Base Model",
##                        "2X Attachment Rate")
## pca_lines_norm$id <- factor(pca_lines_norm$id, levels = c("Base Model", "2X Attachment Rate"))

colorz3 <- c("#666666",
             ## alpha("#e7298a", 0.2),
             ## alpha("#e7298a", 0.5),
             ## alpha("#e7298a", 0.6),
             ## alpha("#e7298a", 0.8),
             "purple")

con_force <- max(data_predict_avg3[curve == 1]$force)

gg_kooiker1 <-
ggplot()+
  ## geom_point(data = data_norm,
  ##            aes(hs_pCa, norm_force,
  ##                color = as.factor(curve)),
  ##            size = 0.75,
  ##            shape = 16)+
  geom_line(data = pca_lines3,
            aes(hs_pCa, Prediction/con_force,
                color = as.factor(curve)),
            linewidth = 2,
            alpha = 0.5)+
  geom_line(data = data_predict_avg3, aes(hs_pCa, force/con_force, color = as.factor(curve)))+
  scale_x_reverse()+
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.25))+
  scale_color_manual(values = colorz3)+
  ylab("Relative Force")+
  xlab("pCa")+
  theme(
    legend.position = "none")



gg_kooiker2 <-
ggplot()+
  ## geom_point(data = data_norm,
  ##            aes(hs_pCa, norm_force,
  ##                color = as.factor(curve)),
  ##            size = 0.75,
  ##            shape = 16)+
  geom_line(data = pca_lines_norm3,
            aes(hs_pCa, Prediction,
                color = as.factor(curve)),
            linewidth = 2,
            alpha = 0.5)+
  geom_line(data = data_fit_avg3, aes(hs_pCa, norm_force, color = as.factor(curve)))+
  scale_x_reverse()+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  scale_color_manual(values = colorz3)+
  ylab("Normalized Force")+
  xlab("pCa")+
  theme(
    legend.position = "none")


## png(filename = "supplemental-figures/supplemental-fibersim-slowed-detachment.png", width = 6.5, height = 2.5, res = 500, units = "in")
## plot_grid(gg_kooiker1, gg_kooiker2, labels = LETTERS)
## dev.off()

## saveRDS(force_pca3, "figures/fibersim-force-pca-kooiker.rds")

#############
## kooiker twitch
###############333

twitch_files  <- list.files("fibersim/project_greenberg_fibersim_kooiker/sim_data/twitch/sim_output",
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

twitch_data <- twitch_data[id2 %in% c(1, 2)]

twitch_data$id3 <- ifelse(twitch_data$id2 == 1, "CON", "Slow Detach")

twitch_data_integral <- twitch_data[, .(integral = sum(hs_1_force*0.001)), by = .(id3, rep) ]
twitch_data_integral_avg <- twitch_data_integral[, .(integral = mean(integral),
                                           sd_integral = sd(integral)),
                                          by = .(id3) ]

twitch_average <- twitch_data[, .(hs_1_force = mean(hs_1_force)), by = .(id2, time)]


twitch_stats <- t.test(twitch_data_integral[id3 == "CON"]$integral,
                       twitch_data_integral[id3 == "Slow Detach"]$integral)

## (
## twitch_integral <-
## ggplot()+
##   ## geom_errorbar(data = twitch_data_integral_avg,
##   ##               aes(id3,
##   ##                   ymin = (integral-sd_integral)/7144.99,
##   ##                   ymax = (integral+sd_integral)/7144.99,
##   ##                   color = id3),
##   ##                   width = 0.25)+
##   geom_quasirandom(data = twitch_data_integral, aes(id3, integral/7144.99, color = id3),
##                    width = 0.2,
##                    shape = 16,
##                    size = 1)+
##   draw_line(x = c(1, 2), y = c(3, 3))+
##   annotate("text", x = 1.5, y = 3.2, label = "p < 0.001", size = 3)+
##   scale_color_manual(values = colorz3)+
##   ## coord_cartesian(ylim = c(0, 3.4))+
##   ylab("Normalized Integral")+
##   xlab("[danicamtiv]")+
##   theme(axis.text.x = element_markdown(),
##         legend.position = "none")
##   )

(
ggtwitch_kooiker <-
ggplot(twitch_data)+
  geom_line(aes(x = time, hs_1_force/39481.58, color = id2), linewidth = 1, alpha = 0.5)+
  geom_line(data = twitch_average, aes(x = time, hs_1_force/39481.58, color = id2), alpha = 1)+
  scale_color_manual(values = colorz3)+
  xlab("Time (s)")+
  ylab("Relative Force")+
  scale_y_continuous(breaks = seq(0, 5, by = 0.5))+
  theme(legend.position = "none"
        ## axis.text.y.left = element_blank(),
        ## axis.ticks.y.left = element_blank(),
        )
)



## ggtwitch1 <-
## ggtwitch+
##   annotate("rect",
##            ymin = 0,
##            ymax = 1.25,
##            xmin = 0.1,
##            xmax = 0.2,
##            color = "black",
##            fill = NA,
##            linetype = "dashed")

## ggtwitch2 <-
## ggtwitch+coord_cartesian(xlim = c(0.1, 0.2), ylim = c(0, 1.25))+
## scale_x_continuous(breaks = c(0.1, 0.15, 0.2))+
##   theme(
##         axis.line = element_blank(),
##         panel.border = element_rect(linetype = "dashed", color = "black")
##   )

## plot_grid(ggtwitch1, ggtwitch2)




png(filename = "supplemental-figures/supplemental-fibersim-slowed-detachment.png", width = 6.5, height = 2.25, res = 500, units = "in")
plot_grid(gg_kooiker1, ggtwitch_kooiker, labels = LETTERS, nrow = 1)
dev.off()


pdf("supplemental-figures/supplemental-fibersim-slowed-detachment.pdf", width = 6.5, height = 2.25)
plot_grid(gg_kooiker1, ggtwitch_kooiker, labels = LETTERS, nrow = 1)
dev.off()
