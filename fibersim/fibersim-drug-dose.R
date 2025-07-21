library(data.table)
library(ggplot2)
library(cowplot)
library(readxl)
library(drc)
library(ggtext)
library(ggbeeswarm)

theme_set(theme_cowplot(11))

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

## normalize_data <- function(data, fmax){
##  data[, .(hs_pCa = hs_pCa,
##           norm_force = hs_force/fmax)]
## }

## normalize_predict <- function(data, fmax){
##  data[, .(hs_pCa = hs_pCa,
##           Prediction = Prediction/fmax)]
## }

dani_data <- as.data.table(read_xlsx("fibersim/dose/project_greenberg_fibersim_dose/sim_data/sim_output/pCa_analysis.xlsx"))
dani_data[, id := "Danicamtiv"]
om_data <- as.data.table(read_xlsx("fibersim/dose/project_om_dose/sim_data/sim_output/pCa_analysis.xlsx"))
om_data[, id := "Omecamtiv Mecarbil"]

data <- rbind(dani_data, om_data)

data[, curve_id := paste0(id, curve)]


pca_nest <- data[, .(data = list(.SD)), by = .(curve, id, curve_id)]


pca_nest[, hill_fit := lapply(data, fit_hill)]
pca_nest[, predict_fit := lapply(hill_fit, predict_hill)]
pca_nest[, mod_table := lapply(hill_fit, broom::tidy)]
pca_nest[, fmax_val := lapply(predict_fit, get_fmax)]
## pca_nest[, norm_data := mapply(normalize_data,
##                                data = data,
##                                fmax = fmax_val,
##                                SIMPLIFY = FALSE)]
## pca_nest[, norm_predict := mapply(normalize_predict,
##                                   data = predict_fit,
##                                   fmax = fmax_val,
##                                   SIMPLIFY = FALSE)]


## pca_lines <- pca_nest[, norm_predict[[1]], by = curve]
## pca_lines_norm <- pca_nest[, predict_fit_norm[[1]], by = id]
pca_lines <- pca_nest[, predict_fit[[1]], by = .(curve, id, curve_id)]
## data_norm <- pca_nest[, norm_data[[1]], by = .(curve, id)]


## data_fit_avg <- pca_lines_norm[, .(norm_force = mean(Prediction)), by = .(curve, hs_pCa)]
## data_norm$id <- ifelse(data_norm$curve==1, "Base Model",
##                        "2X Attachment Rate")

## data_norm$id <- factor(data_norm$id, levels = c("Base Model", "2X Attachment Rate"))

## pca_lines_norm$id <- ifelse(pca_lines_norm$curve==1, "Base Model",
##                        "2X Attachment Rate")
## pca_lines_norm$id <- factor(pca_lines_norm$id, levels = c("Base Model", "2X Attachment Rate"))

colorz1 <- c("#666666",
             alpha("#e7298a", 0.25),
             alpha("#e7298a", 0.5),
             alpha("#e7298a", 0.75),
             "#e7298a",
             "#666666",
             alpha("#d95f02", 0.25),
             alpha("#d95f02", 0.5),
             alpha("#d95f02", 0.75),
             "#d95f02")

colorz2 <- colorz1[6:10]

(
force_pca <-
ggplot()+
  ## geom_point(data = data_norm,
  ##            aes(hs_pCa, norm_force,
  ##                color = as.factor(curve)),
  ##            size = 0.75,
  ##            shape = 16)+
  geom_line(data = pca_lines,
            aes(hs_pCa, Prediction/max(Prediction),
                color = as.factor(curve_id)),
            linewidth = 1)+
  ## geom_line(data = data_fit_avg, aes(hs_pCa, norm_force, color = as.factor(curve)))+
  facet_wrap(~id, scales = "free_y")+
  scale_x_reverse()+
  ## scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  scale_color_manual(values = colorz1)+
  ylab("Relative Force")+
  xlab("pCa")+
  theme(
    legend.position = "none",
    strip.background = element_blank()
    )
)


## png("~/Downloads/fibersim-dani-om-dose.png", bg = "white", width = 3, height = 3, units = "in")
## force_pca
## dev.off()

## cairo_pdf("~/Downloads/fibersim-om-only-dose.pdf", bg = "white", width = 3.25, height = 2.79)
## force_pca
## dev.off()

## ggsave("~/Downloads/fibersim-dani-om-dose.png", bg = "white")


###################################
dani_twitch_files  <- list.files("fibersim/dose/project_greenberg_fibersim_dose_twitch/sim_data/twitch/sim_output",
                            pattern = "sim_prot",
                            recursive = TRUE,
                            full.names = TRUE)


om_twitch_files  <- list.files("fibersim/dose/project_om_dose_twitch2/sim_data/twitch/sim_output",
                            pattern = "sim_prot",
                            recursive = TRUE,
                            full.names = TRUE)

twitch_files <- c(dani_twitch_files, om_twitch_files)

fread_twitch <- function(x){
  dat <- fread(x)
  dat$id <- x
    return(dat)
}

twitch_data <- rbindlist(lapply(twitch_files, fread_twitch))
twitch_data[, c("fibersim", "dose", "id2", "sd", "twitch", "so", "conc", "file") := tstrsplit(id,"/")]

twitch_data[, c("sim", "prot", "num", "num2", "rep") := tstrsplit(file,"_")]

twitch_data[, id3 := ifelse(id2 == "project_greenberg_fibersim_dose_twitch", "Danicamtiv", "Omecamtiv Mecarbil")]

twitch_data[, id_conc := paste0(id3, conc)]

## max_control <- max(twitch_data[id2 == 1]$hs_1_force)

## twitch_data <- twitch_data[id2 %in% c(1, 6)]

## twitch_data$id3 <- ifelse(twitch_data$id2 == 1, "0 &micro;M", "10 &micro;M")

## twitch_data_integral <- twitch_data[, .(integral = sum(hs_1_force*0.001)), by = .(id3, rep) ]
## twitch_data_integral_avg <- twitch_data_integral[, .(integral = mean(integral),
##                                            sd_integral = sd(integral)),
##                                           by = .(id3) ]

twitch_average <- twitch_data[, .(hs_1_force = mean(hs_1_force)), by = .(id3, id_conc, time, conc)]

twitch_relative <- max(twitch_average[id_conc == "Danicamtiv1"]$hs_1_force)

## twitch_stats <- t.test(twitch_data_integral[id3 == "0 &micro;M"]$integral,
##                        twitch_data_integral[id3 == "10 &micro;M"]$integral)

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
##   scale_color_manual(values = colorz1)+
##   coord_cartesian(ylim = c(0, 3.4))+
##   ylab("Normalized Integral")+
##   xlab("[danicamtiv]")+
##   theme(axis.text.x = element_markdown(),
##         legend.position = "none")
##   )

colorz_twitch <- c("#666666",
             alpha("#e7298a", 0.35),
             alpha("#e7298a", 0.45),
             alpha("#e7298a", 0.55),
             alpha("#e7298a", 0.65),
             alpha("#e7298a", 0.75),
             alpha("#e7298a", 0.85),
             "#e7298a",
             "#666666",
             alpha("#d95f02", 0.35),
             alpha("#d95f02", 0.45),
             alpha("#d95f02", 0.55),
             alpha("#d95f02", 0.65),
             alpha("#d95f02", 0.75),
             alpha("#d95f02", 0.85),
             "#d95f02")

twitch_colorz2 <- c(tail(RColorBrewer::brewer.pal(8, "Dark2"), 1), head(RColorBrewer::brewer.pal(8, "Dark2"), 7))


twitch_average[, conc2 := ifelse(conc == 1, "0%",
                                 ifelse(conc == 2, "5%",
                                        ifelse(conc == 3, "10%",
                                        ifelse(conc == 4, "20%",
                                               ifelse(conc == 5, "25%",
                                                      ifelse(conc == 6, "50%",
                                                             ifelse(conc == 7, "75%",
                                                                    ifelse(conc == 8, "100%", NA))))))))]

twitch_average$conc2 <- factor(twitch_average$conc2, levels = c("0%",
                                                                "5%",
                                                                "10%",
                                                                "20%",
                                                                 "25%",
                                                                "50%",
                                                                "75%",
                                                                "100%"))

ggtwitch <-
ggplot()+
  ## geom_line(data = twitch_data, aes(x = time, hs_1_force/twitch_relative, color = id_conc), linewidth = 1, alpha = 0.5)+
  geom_line(data = twitch_average, aes(x = time, hs_1_force/twitch_relative, color = conc2), linewidth = 1)+
  scale_color_manual(name = "% Myosin \n Heads", values = rev(rainbow(8)))+
  facet_wrap(~id3, scales = "free")+
  coord_cartesian(ylim = c(0, 1.5))+
  xlab("Time (s)")+
  ylab("Relative Force")+
  ## scale_y_continuous(breaks = seq(0, 1.25, by = 0.25))+
  theme(strip.background = element_blank(),
        ## legend.position = "none"
        )


png("supplemental-figures/fibersim-dose.png", width = 6, height = 2.5, units = "in", res = 300)
## plot_grid(force_pca, ggtwitch, nrow = 2, labels = LETTERS)
force_pca
dev.off()
