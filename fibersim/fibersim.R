########
# tension pca
#########

library(data.table)
library(ggplot2)
library(cowplot)
library(readxl)
library(drc)
library(ggtext)

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

data <- as.data.table(read_xlsx("fibersim/project_greenberg_fibersim_dose/sim_data/sim_output/pCa_analysis.xlsx"))

data <- data[curve %in% c(1,6)]

pca_nest <- data[, .(data = list(.SD)), by = curve]


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



## pca_lines <- pca_nest[, norm_predict[[1]], by = curve]
## pca_lines_norm <- pca_nest[, predict_fit_norm[[1]], by = id]
pca_lines_norm <- pca_nest[, norm_predict[[1]], by = curve]
data_norm <- pca_nest[, norm_data[[1]], by = curve]

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

force_pca <-
ggplot()+
  geom_point(data = data_norm,
             aes(hs_pCa, norm_force,
                 color = as.factor(curve)),
             size = 1.5,
             shape = 16)+
  geom_line(data = pca_lines_norm,
            aes(hs_pCa, Prediction,
                color = as.factor(curve)),
            linewidth = 0.8)+
  scale_x_reverse()+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  scale_color_manual(values = colorz1)+
  ylab("Normalized Force")+
  xlab("pCa")+
  theme(
    legend.position = "none")

saveRDS(force_pca, "figures/fibersim-force-pca.rds")

################################3
##################################

twitch_files  <- list.files("fibersim/project_greenberg_fibersim_dose_twitch/sim_data/twitch/sim_output",
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

max_control <- max(twitch_data[id2 == 1]$hs_1_force)

twitch_data <- twitch_data[id2 %in% c(1, 6)]

twitch_data$id3 <- ifelse(twitch_data$id2 == 1, "0 &micro;M", "10 &micro;M")

twitch_data_integral <- twitch_data[, .(integral = sum(hs_1_force*0.001)), by = id3 ]

(
twitch_integral <-
ggplot(twitch_data_integral)+
  geom_col(aes(id3, integral/min(integral), fill = id3))+
  scale_fill_manual(values = colorz1)+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.01)))+
  ylab("Normalized Integral")+
  xlab("[danicamtiv]")+
  theme(axis.text.x = element_markdown(),
        legend.position = "none")
  )

ggtwitch <-
ggplot(twitch_data)+
  geom_line(aes(x = time, hs_1_force/44409.2, color = id2), size = 0.8 )+
  scale_color_manual(values = colorz1)+
  xlab("Time (s)")+
  ylab("Normalized Force")+
  scale_y_continuous(breaks = seq(0, 1.25, by = 0.25))+
  theme(legend.position = "none"
        ## axis.text.y.left = element_blank(),
        ## axis.ticks.y.left = element_blank(),
        )

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

plot_grid(ggtwitch1, ggtwitch2)

saveRDS(ggtwitch, "figures/fibersim-ggtwitch.rds")
saveRDS(ggtwitch1, "figures/fibersim-ggtwitch1.rds")
saveRDS(ggtwitch2, "figures/fibersim-ggtwitch2.rds")
saveRDS(twitch_integral, "figures/fibersim-twitch-integral.rds")





####################



















ggsave("fibersim-force-pca.png")









myo_nest <- myosim_nest[fit_nest, on = "id"]

data_nest[, max_fit := lapply(fit, max_fit)]
data_nest[, norm_data := Map(normalize_data, x = data, max = max_fit)]
data_nest[, norm_fit := Map(normalize_fit, x = fit, max = max_fit)]


myo_dat <- myo_nest[, norm_data[[1]], by = id]
fit_dat <- myo_nest[, norm_fit[[1]], by = id]

max_fit <- function(x){
  max(x$y)
}

normalize_data <- function(x, max){
 x[, .(pCa,
       y = y/max)]
}


normalize_fit <- function(x, max){
 x[, .(x,
       y = y/max)]
}

myo_nest[, max_fit := lapply(fit, max_fit)]
myo_nest[, norm_data := Map(normalize_data, x = data, max = max_fit)]
myo_nest[, norm_fit := Map(normalize_fit, x = fit, max = max_fit)]


myo_dat <- myo_nest[, norm_data[[1]], by = id]
fit_dat <- myo_nest[, norm_fit[[1]], by = id]
## myosim$y1 <- myosim$y
## fit$y2 <- fit$y
## fit$pCa <- fit$x

## myosim <- myosim[fit, on = c("id", "pCa")]

## myosim[, norm_y := y/max(y, na.rm = TRUE), by = id]
## myosim[, norm_y2 := y2/max(y2), by = id]

## myosim <- myosim[id %in% c("5-XPS_200-k3_0.001-KCB_1-k1", "5-XPS_100-k3_0.001-KCB_1-k1")]
## fit <- fit[id %in% c("5-XPS_200-k3_0.001-KCB_1-k1", "5-XPS_100-k3_0.001-KCB_1-k1")]

ggplot()+
  geom_point(data = myo_dat, aes(pCa, y, color = id))+
  geom_line(data = fit_dat, aes(x, y, color = id))+
  scale_x_reverse()+
  scale_color_manual(values = c("#666666", "#e7298a"))+
  ylab("Muscle Force")+
  theme(
  legend.position = "right",
  legend.text = element_text(size =  8))

ggsave("~/Downloads/myosim-tension-pca.pdf", bg = "white")

##################
#############################3
###########################
####
#### SUPPLEMENTAL
#####################3
########################


data2 <- as.data.table(read_xlsx("fibersim/project_greenberg_fibersim/sim_data/force_pCa/sim_output/pCa_analysis.xlsx"))


pca_nest2 <- data2[, .(data = list(.SD)), by = curve]


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



pca_lines <- pca_nest2[, predict_fit[[1]], by = curve]
## pca_lines_norm <- pca_nest[, predict_fit_norm[[1]], by = id]
pca_lines_norm <- pca_nest2[, norm_predict[[1]], by = curve]
data_norm <- pca_nest2[, norm_data[[1]], by = curve]

## data_norm$id <- ifelse(data_norm$curve==1, "Base Model",
##                        "2X Attachment Rate")

## data_norm$id <- factor(data_norm$id, levels = c("Base Model", "2X Attachment Rate"))

## pca_lines_norm$id <- ifelse(pca_lines_norm$curve==1, "Base Model",
##                        "2X Attachment Rate")
## pca_lines_norm$id <- factor(pca_lines_norm$id, levels = c("Base Model", "2X Attachment Rate"))

## colorz2 <- c("#666666",
##              ## alpha("#e7298a", 0.2),
##              ## alpha("#e7298a", 0.5),
##              ## alpha("#e7298a", 0.6),
##              ## alpha("#e7298a", 0.8),
##              "#e7298a")

dark2 <- RColorBrewer::brewer.pal(8, "Dark2")
colorz2 <- c("#666666", dark2[1:3], dark2[5], dark2[4])

(
gg_force_pca_params <-
ggplot()+
  geom_point(data = data2,
             aes(hs_pCa, hs_force,
                 color = as.factor(curve)),
             size = 1.5,
             shape = 16)+
  geom_line(data = pca_lines,
            aes(hs_pCa, Prediction,
                color = as.factor(curve)),
            linewidth = 0.8)+
  scale_x_reverse()+
  ## scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  scale_color_manual(name = "Simulation", values = colorz2)+
  ylab("Simulated Force")+
  xlab("pCa")+
  theme(legend.position = c(0.1, 0.5)
        )
)

## saveRDS(force_pca, "figures/fibersim-force-pca.rds")

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

## max_control <- max(twitch_data[id2 == 1]$hs_1_force)

## twitch_data <- twitch_data[id2 %in% c(1, 6)]

twitch_data_integral2 <- twitch_data2[, .(integral = sum(hs_1_force*0.001)), by = id2 ]

con_to_remove <- twitch_data_integral2[id2 == 1]$integral

ggintegral2 <-
ggplot(twitch_data_integral2)+
  geom_col(aes(id2, integral/con_to_remove, fill = id2))+
  scale_fill_manual(values = colorz2)+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.01)))+
  ylab("Normalized Integral")+
  xlab("Simulation")+
  theme(legend.position = "none")

gg_twitch_params <-
ggplot(twitch_data2)+
  geom_line(aes(x = time, hs_1_force, color = id2), size = 0.8 )+
  ## scale_color_manual(values = colorz1)+
  xlab("Time (s)")+
  ylab("Simulated Force")+
  scale_color_manual(name = "Simulation", values = colorz2)+
  theme(legend.position = "none")




png(filename = "supplemental-figures/supplemental-fibersim-5.png", width = 6.5, height = 2.5, res = 500, units = "in")
plot_grid(gg_force_pca_params, gg_twitch_params, ggintegral2,
          rel_widths = c(1, 1, 0.5),
          nrow = 1,
          labels = c("A", "B", "C"))
dev.off()








  ## scale_y_continuous(breaks = seq(0, 1.25, by = 0.25))+
  theme(
        ## axis.text.y.left = element_blank(),
        ## axis.ticks.y.left = element_blank(),
        )

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

plot_grid(ggtwitch1, ggtwitch2)

saveRDS(ggtwitch, "figures/fibersim-ggtwitch.rds")
saveRDS(ggtwitch1, "figures/fibersim-ggtwitch1.rds")
saveRDS(ggtwitch2, "figures/fibersim-ggtwitch2.rds")
saveRDS(twitch_integral, "figures/fibersim-twitch-integral.rds")
