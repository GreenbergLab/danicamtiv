library(data.table)
library(ggplot2)
library(cowplot)
library(minpack.lm)
library(ggtext)
library(RcppRoll)
library(magick)
library(ggbeeswarm)
library(ungeviz)

theme_set(theme_cowplot(10))
#read data
folder <- file.path("stopped-flow", "detach-reattach")

files <- list.files(path = folder,
                    pattern = ".csv",
                    recursive = TRUE)

exclude <- c("2024-03-15/bio-rep-1/10uM-dani/10uM-actin/10uM-dani_5uM-myo_10uM-actin_3.csv",
             "2024-03-19/bio-rep-1/10uM-dani/20uM-actin/10uM-dani_5uM-myo_20uM-actin_4.csv",
             "2024-03-15/bio-rep-1/0uM-dani/15uM-actin/control_5uM-myo_15uM-actin_4.csv",
             "2024-03-18/bio-rep-1/0uM-dani/5uM-actin/dmso_5uM-myo_5uM-actin_0.csv")

files <- files[-c(which(files %in% exclude))]

files <- file.path(folder, files)


read_files <- function(x){
   d <- fread(x)
   d$path <- x
   return(d)
}

data <- rbindlist(lapply(files, read_files))

id_split <- c("sf", "dere", "date", "rep", "id", "actin", "filename")

data[, (id_split) := tstrsplit(path, "/", fixed = TRUE)]
data[, c("other", "filenum") := tstrsplit(filename, "actin_", fixed = TRUE)]
data[, filenum := as.factor(sub(".csv", "", filenum))]
data[, id2 := sub("uM-dani", " &micro;M", id)]

data[, `:=`(x_time_s = V1,
            y_signal = V2)]

data <- data[, .(id, id2, date, rep, actin, filenum, x_time_s, y_signal)]

#nest data
data_nest <- data[, list(spectra = list(.SD)), by = c("id", "id2", "date", "rep", "actin", "filenum")]

#prepare deta
#filter tmin
#normalize data
prep_data1 <- function(x){
  x <- x[1:10000]
  x <- x[x_time_s >= 0.020]
  x[, x_time_s := x_time_s - 0.020]
  x[, y_signal := y_signal-min(y_signal) ]
  x[, y_signal := y_signal/max(y_signal) ]
  return(x)
}

prep_data2 <- function(x){
  x <- x[1:10000]
  x1 <- x[1:5000]
  x2 <- x[5001:10000]
  x1_roll_mean <- na.omit(roll_meanl(x1$y_signal, n = 100, by = 100))
  x3 <- data.table(x_time_s = seq(0.02, 1, by = 0.02),
                   y_signal = x1_roll_mean
                   )
  x <- rbind(x3, x2)
  ## x <- x[c(seq(1, 5000, by = 100), 5001:10000)]
  x <- x[x_time_s >= 0.0205]
  x[, x_time_s := x_time_s - 0.0205]
  x[, y_signal := y_signal-min(y_signal) ]
  x[, y_signal := y_signal/max(y_signal) ]
  return(x)
}

data_nest[, spectra1 := lapply(spectra, prep_data1)]
data_nest[, spectra2 := lapply(spectra, prep_data2)]

fit_mod <- function(x){
  nlsLM(y_signal ~ (a1*exp(-k1*x_time_s)+(a2*exp(k2*x_time_s))+c),
       data = x,
       start = list(
         ## a = 1,
         a1 = 0.5,
         a2 = 0.5,
         c = 0.5,
         k1 = 5,
         k2 = -0.02),
       control = nls.lm.control(maxiter = 1024)
       )
}

predict_line <- function(mod){
  x <- seq(0, 101, by = 0.001)
  y <- predict(mod, newdata = data.frame(x_time_s = x))
  return(
    data.frame(x_time_s = x,
               y_signal = y,
               type = "prediction")
  )
}

get_k1 <- function(mod){
  coef(mod)[["k1"]]
}

get_k2<- function(mod){
  coef(mod)[["k2"]]
}


data_nest[, mod := lapply(spectra1, fit_mod)]
data_nest[, predict_df := lapply(mod, predict_line)]
data_nest[, coef_k1 := sapply(mod, get_k1)]
data_nest[, coef_k2 := sapply(mod, get_k2)]

data_predict <- data_nest[, predict_df[[1]], by = c("id", "id2", "date", "rep", "actin", "filenum")]

data_unnest <- data_nest[, spectra1[[1]], by = c("id","id2", "date", "rep", "actin", "filenum")]
data_unnest$type <- "real"

data_unnest2 <- data_nest[, spectra2[[1]], by = c("id", "id2", "date", "rep", "actin", "filenum")]

data_together <- rbind(data_unnest, data_predict)
data_together$type <- factor(data_together$type, levels = c("real", "prediction"))

## ggplot(data = data_unnest[id != "100uM-dani"],
##        aes(x = x_time_s,
##            y = y_signal,
##            color = filenum))+
##   geom_line(linewidth=0.2)+
##   facet_grid(date+rep+actin~id, scales = "free")+
##   scale_x_log10()+
##   theme_cowplot()+
##   theme(
##   strip.text = element_text(size = 6),
##   legend.position = "bottom")



## ggplot(data = data_together[id != "100uM-dani"],
##        aes(x = x_time_s,
##            y = y_signal,
##            color = filenum,
##            linetype = type,
##            alpha = type))+
##   geom_line(linewidth=0.2)+
##   facet_grid(date+rep+actin~id, scales = "free")+
##   scale_linetype_manual(values = c(:))+
##   scale_alpha_manual(values = c(1, 0.25))+
##   scale_x_log10()+
##   theme_cowplot()+
##   theme(
##     strip.text = element_text(size = 6),
##     legend.position = "bottom")


data_summary <- data_nest[, .(avg_k1 = mean(coef_k1),
                              sd_k1 = sd(coef_k1),
                              avg_k2 = mean(coef_k2),
                              sd_k2 = sd(coef_k2)),
                          by = c("id", "id2", "actin")]

setorder(data_summary, id, actin)



colors <- c("#666666", "#e7298a")

shapiro.test(data_nest[id == "0uM-dani"]$coef_k2)
shapiro.test(data_nest[id == "10uM-dani"]$coef_k2)

shapiro.test(data_nest[id == "0uM-dani"]$coef_k1)
shapiro.test(data_nest[id == "10uM-dani"]$coef_k1)

attach_pval <- kruskal.test(x = data_nest$coef_k2, data_nest$id)

detach_pval <- kruskal.test(x = data_nest$coef_k1, data_nest$id)

gg_attach <-
ggplot()+
  geom_col(data = data_summary[actin == "20uM-actin"],
           aes(x = id2,
               y = abs(avg_k2),
               fill = id2))+
  geom_errorbar(data = data_summary[actin == "20uM-actin"],
                aes(x = id2,
                    ymin = abs(avg_k2)-abs(sd_k2),
                    ymax = abs(avg_k2)+abs(sd_k2),
                    color = id2
                    ),
                width = 0.25)+
  geom_jitter(data = data_nest[actin == "20uM-actin"],
              aes(x = id2,
                  y = abs(coef_k2)),
              width = 0.1,
              size = 0.5)+
  ylab("<i>k<sub>att</sub></i> (s<sup>-1</sup>)")+
  xlab("[danicamtiv]")+
  draw_line(c(1, 2), y = 0.065)+
  annotate("text", x = 1.5, y = 0.07, label = paste0("p = ", round(attach_pval$p.value, 5)), size = 3)+
  ## ggtitle("Rate of Re-Attachment")+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0)))+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  theme_cowplot(11)+
  theme(
  axis.title.y = element_markdown(),
  axis.text.x = element_markdown(),
  legend.position = "none")
gg_attach



gg_detach <-
ggplot()+
  geom_col(data = data_summary[actin == "20uM-actin"],
           aes(x = id2,
               y = abs(avg_k1),
               fill = id2))+
  geom_errorbar(data = data_summary[actin == "20uM-actin"],
                aes(x = id2,
                    ymin = abs(avg_k1)-abs(sd_k1),
                    ymax = abs(avg_k1)+abs(sd_k1),
                    color = id2
                    ),
                width = 0.25)+
  geom_jitter(data = data_nest[actin == "20uM-actin"],
              aes(x = id2,
                  y = abs(coef_k1)),
              width = 0.1,
              size = 0.5)+
  ylab("<i>k<sub>det</sub></i> (s<sup>-1</sup>)")+
  xlab("[danicamtiv]")+
  draw_line(c(1, 2), y = 8.5)+
  annotate("text", x = 1.5, y = 9.25, label = paste0("p = ", round(detach_pval$p.value, 2)), size = 3)+
  ## ggtitle("Rate of Detachment")+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0)))+
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  theme_cowplot(11)+
  theme(
  axis.title.y = element_markdown(),
  axis.text.x = element_markdown(),
  legend.position = "none")
gg_detach



gg_transients <-
ggplot()+
  geom_line(data = data_unnest[id == "0uM-dani" &
                               actin == "20uM-actin" &
                               date == "2024-03-18" &
                               rep == "bio-rep-1" &
                               filenum == 1 #&
                               ## x_time_s <=20 &
                               ## y_signal >= 0
                               ],
            aes(x = x_time_s,
                y = y_signal),
            color = "#666666",
            alpha = 1,
            linewidth = 0.5)+
  geom_line(data = data_unnest[id == "10uM-dani" &
                               actin == "20uM-actin" &
                               date == "2024-03-18" &
                               rep == "bio-rep-1" &
                               filenum == 1 #&
                               ## x_time_s <= 20 &
                               ## y_signal >= 0
                               ],
            aes(x = x_time_s,
                y = y_signal),
            color = "#e7298a",
            alpha = 0.7,
            linewidth = 0.5)+
  ## geom_line(data = data_predict[id == "10uM-dani" &
  ##                              actin == "20uM-actin" &
  ##                              date == "2024-03-18" &
  ##                              rep == "bio-rep-1" &
  ##                              filenum == 1 #&
  ##                              ## x_time_s <= 20 &
  ##                              ## y_signal >= 0
  ##                              ],
  ##           aes(x = x_time_s,
  ##               y = y_signal),
  ##           color = "red",
  ##           alpha = 0.7,
  ##           linewidth = 0.5)+
  ylab("Normalized <br> Fluorescence")+
  xlab("Time (s)")+
  ## coord_cartesian(xlim = c(0, 20))+
  theme_cowplot(10)+
  theme(
    axis.title.y = element_markdown()
    ## axis.text.y = element_blank(),
    ## axis.ticks.y = element_blank()
  )


## gg_transients2 <-
## ggplot()+
##   geom_line(data = data_unnest[id == "0uM-dani" &
##                                actin == "20uM-actin" &
##                                date == "2024-03-18" &
##                                rep == "bio-rep-1" &
##                                filenum == 1 #&
##                                ## x_time_s <=20 &
##                                ## y_signal >= 0
##                                ],
##             aes(x = x_time_s,
##                 y = y_signal),
##             color = "#666666",
##             alpha = 0.5,
##             linewidth = 1)+
##   geom_line(data = data_predict[id == "0uM-dani" &
##                                actin == "20uM-actin" &
##                                date == "2024-03-18" &
##                                rep == "bio-rep-1" &
##                                filenum == 1 #&
##                                ## x_time_s <=20 &
##                                ## y_signal >= 0
##                                ],
##             aes(x = x_time_s,
##                 y = y_signal),
##             color = "black",
##             alpha = 1,
##             linewidth = 0.5,
##             linetype = "dashed")+
##   geom_line(data = data_unnest[id == "10uM-dani" &
##                                actin == "20uM-actin" &
##                                date == "2024-03-18" &
##                                rep == "bio-rep-1" &
##                                filenum == 2 #&
##                                ## x_time_s <= 20 &
##                                ## y_signal >= 0
##                                ],
##             aes(x = x_time_s,
##                 y = y_signal),
##             color = "#e7298a",
##             alpha = 0.5,
##             linewidth = 1)+
##   geom_line(data = data_predict[id == "10uM-dani" &
##                                actin == "20uM-actin" &
##                                date == "2024-03-18" &
##                                rep == "bio-rep-1" &
##                                filenum == 2 #&
##                                ## x_time_s <= 20 &
##                                ## y_signal >= 0
##                                ],
##             aes(x = x_time_s,
##                 y = y_signal),
##             color = "red",
##             alpha = 1,
##             linewidth = 0.5,
##             linetype = "dashed")+
##   ylab("Relative Fluorescence")+
##   xlab("log(Time (s))")+
##   coord_cartesian(
##     xlim = c(0.001, 100),
##                   ylim = c(0, NA))+
##   scale_x_log10()+
##   theme_cowplot(8)+
##   theme(
##     axis.title.y = element_blank(),
##     axis.ticks.y = element_blank(),
##     axis.line.y = element_blank(),
##     axis.text.y = element_blank()
##     ## axis.text.y = element_blank()
##     ## axis.text.y = element_blank(),
##     ## axis.ticks.y = element_blank()
##   )

## gg50 <-
## gg_transients <-
## ggplot(data = data_together[id == "10uM-dani" & actin == "20uM-actin" & date == "2024-03-19" & filenum==0],
##        aes(x = x_time_s,
##            y = y_signal,
##            color = filenum,
##            linetype = type,
##            alpha = type,
##            linewidth = type))+
##   geom_line()+
##   ylab("Relative Fluorescence")+
##   xlab("Time (s)")+
##   facet_grid(date+rep~id+filenum)+
##   scale_linetype_manual(values = c("solid", "solid"))+
##   scale_alpha_manual(values = c(0.5, 1))+
##   scale_linewidth_manual(values = c(1, 0.5))+
##   ## scale_x_log10()+
##   ## coord_cartesian(xlim = c(0, 10))+
##   theme_cowplot()+
##   theme(
##     strip.text = element_text(size = 12),
##     ## legend.position = "bottom",
##     strip.background = element_rect(fill = "white"),
##     legend.position = "none")


## plot_grid(gg100, gg50)



arrowup <- data_predict[id == "0uM-dani" &
                               actin == "20uM-actin" &
                               date == "2024-03-18" &
                               rep == "bio-rep-1" &
                               filenum == 1 #&
                               ## x_time_s <=20 &
                               ## y_signal >= 0
                        ]
arrowupy <- arrowup[15]$y_signal
arrowupx <- arrowup[15]$x_time_s
arrowupy2 <- arrowup[750]$y_signal
arrowupx2 <- arrowup[750]$x_time_s

gg_example <-
ggplot()+
  geom_line(
            aes(x = c(0, arrowup$x_time_s),
                y = c(0, arrowup$y_signal)),
            linewidth = 0.5)+
  annotate("richtext",
           x = 0.1,
           y = 0.45,
           label = "<i>k<sub>det</sub></i>",
           hjust = 0.5,
           vjust = 1,
           size = 10/.pt,
           label.padding = unit(c(0.1, 0.1, 0.1, 0.1), "lines"),
           label.margin = unit(c(0.05, 0.05, 0.05, 0.05), "lines"))+
  annotate("point",
           x = 0,
           y = 0.55,
           shape = 17,
           size = 1.5)+
  annotate("richtext",
           x = 30,
           y = 0.45,
           label = "<i>k<sub>att</sub></i>",
           hjust = 0.5,
           vjust = 1,
           label.padding = unit(c(0.1, 0.1, 0.1, 0.1), "lines"),
           label.margin = unit(c(0.05, 0.05, 0.05, 0.05), "lines"),
           size = 10/.pt)+
  annotate("point",
           x = 40,
           y = 0.28,
           shape = 17,
           size = 1.5)+
  theme_void()
gg_example

gg_example2 <-
  gg_example+
  coord_cartesian(xlim = c(-22, NA), ylim = c(-0.1, 1.1))+
  ## draw_line(x = c(-8, -8),
  ##           y = c(0.6, 0.8),
  ##           linetype = "solid",
  ##           linewidth = 0.8,
  ##           arrow = arrow(length = unit(0.3, "cm"), type = "closed"))+
  draw_image("detach-attach-1.1.png",
             x = -30,
             y = 0.1,
             width = 25,
             hjust = 0,
             vjust = 0.5)+
  draw_image("detach-attach-1.2.png",
             x = -30,
             y = 0.5,
             width = 25,
             hjust = 0,
             vjust = 0.5)+
  draw_image("detach-attach-1.3.png",
             x = 0,
             y = 0.49,
             width = 25,
             hjust = 0.5)+
  draw_image("detach-attach-1.2.png",
             x = 45,
             y = 0.5,
             width = 25,
             hjust = 0.5,
             vjust= 0.5)+
  draw_image("detach-attach-1.1.png",
             x = 75,
             y = 0.3,
             width = 25,
             hjust = 0.5,
             vjust = 0.5)
  ## coord_cartesian(xlim = c(-22, NA))
gg_example2

## gg_example2




gg_side <- plot_grid(gg_detach,
                     gg_attach,
                     ncol = 2,
                     nrow = 1,
                     rel_widths = c(1, 1),
                     labels = c("C", "D"))
gg_side

gg_top_right <-
ggdraw(gg_transients)+
  draw_plot(gg_transients2,
            x = 0.5,
            y = 0.5,
            width = 0.45,
            height = 0.45)
gg_top_right

gg_top_left <-
ggdraw(gg_example2)
  ## draw_plot(gg_side,
  ##           x = 0.4,
  ##           y = 0.5,
  ##           width = 0.6,
  ##           height = 0.5)
gg_top_left

gg_top <- plot_grid(gg_top_left, gg_top_right, nrow = 1, labels = c("A", "B"))
gg_top


## gg_hydro <- readRDS("figures/atp-hydrolysis.rds")


gg_bottom <- plot_grid(gg_side, gg_hydro, nrow = 1, rel_widths = c(1, 1.25))
gg_bottom

## plot_grid(gg_top, gg_bottom, nrow = 2, rel_heights = c(1.25, 1))

gg_transients+draw_plot(gg_example, height = 0.5, width = 0.5)

bottom <- plot_grid(gg_transients, gg_detach, gg_attach, nrow = 1, rel_widths = c(1, 0.5, 0.5), labels = c("B", "C", "D"))

plot_grid(gg_example2, bottom, nrow = 2, rel_heights = c(1, 1.5))
## plot_grid(gg_example2, gg_transients, gg_side, rel_widths = c(1.25, 1, 0.75), nrow = 1)
## pdf("figures/detach-attach.pdf", width = 6.5, height = 3)
## plot_grid(gg_left, gg_right, rel_widths = c(0.9, 1), nrow = 1, labels = c("A", "D"), label_size = 10)
## dev.off()
#############################33
###### SIMULATION FITS ##########
#############################

sim_dat <- readxl::read_excel("stopped-flow/detach-reattach/simulation-fit-values-detach-atttach.xlsx", sheet = 2)

## ggdetach_sim <-

sim_dat$id <- ifelse(sim_dat$Dani == 0, "0 <br> &micro;M", "10 <br> &micro;M")

sim_dat_sum <- sim_dat |>
  dplyr::group_by(Dani, id) |>
  dplyr::summarise(k1_avg = mean(k1),
                   k1_sd = sd(k1),
                   k2_avg = mean(k2),
                   k2_sd = sd(k2))



con_sim_data <- dplyr::filter(sim_dat, Dani == 0)
dani_sim_data <- dplyr::filter(sim_dat, Dani == 10)
detach_pval2 <- t.test(con_sim_data$k1, dani_sim_data$k1)
attach_pval2 <- t.test(con_sim_data$k2, dani_sim_data$k2)


ggdet2 <-
ggplot()+
  geom_quasirandom(data = sim_dat, aes(id, k1, color = id), size = 1, width = 0.25)+
  geom_errorbar(data = sim_dat_sum, aes(id, ymin = k1_avg-k1_sd, ymax = k1_avg+k1_sd, color = id),
                width = 0.4, size = 1)+
  geom_hpline(data  = sim_dat_sum, aes(id, k1_avg, color = id), width = 0.3)+
  draw_line(c(1, 2), y = 6.2)+
  annotate("text", x = 1.5, y = 6.5, label = paste0("ns"), size = 6/.pt)+
  coord_cartesian(ylim = c(0, 6.28))+
  scale_color_manual(values = colors)+
  theme_cowplot(8)+
  xlab("[danicamtiv]")+
  ylab("<i>k<sub>det</sub></i> (&micro;M<sup>-1</sup>&middot; s<sup>-1</sup>)")+
  theme(
  axis.text.x = element_markdown(),
  axis.title.y = element_markdown(),
  legend.position = "none"
  )



ggatt2 <-
  ggplot()+
  geom_quasirandom(data = sim_dat, aes(id, k2, color = id), size = 1, width = 0.25)+
  geom_errorbar(data = sim_dat_sum, aes(id, ymin = k2_avg-k2_sd, ymax = k2_avg+k2_sd, color = id), width = 0.4, size = 1)+
  geom_hpline(data  = sim_dat_sum, aes(id, k2_avg, color = id), width = 0.3)+
  draw_line(c(1, 2), y = 0.008)+
  annotate("text", x = 1.5, y = 0.0084, label = paste0(" p < 0.001 "), size = 6/.pt)+
  coord_cartesian(ylim = c(0, 0.0084))+
  scale_color_manual(values = colors)+
  theme_cowplot(8)+
  xlab("[danicamtiv]")+
  ylab("<i>k<sub>att</sub></i> (&micro;M<sup>-1</sup>&middot; s<sup>-1</sup>)")+
  theme(
  axis.text.x = element_markdown(),
  axis.title.y = element_markdown(),
  legend.position = "none"
  )


bottom2 <- plot_grid(gg_transients, ggdet2, ggatt2, nrow = 1, rel_widths = c(1, 0.5, 0.6), labels = c("B", "C", "D"))
bottom3 <- plot_grid(ggdet2, ggatt2, nrow = 1, rel_widths = c(1, 1), labels = c("C", "D"))

(bottom_plus <- plot_grid(gg_transients, bottom3, ncol = 1, rel_heights = c(1,1), labels = c("B", "")))


## ggsave("~/Downloads/detach-reattach-plot.pdf", bg = "white")


png("figures/figure-4-new.png", width = 3.42, height = 5, units = "in", res = 300)
ggdraw(plot_grid(gg_example2, gg_transients, ncol = 1, rel_heights = c(1, 1.5), labels = c("A", "B")))+
  draw_plot(bottom3, x = 0.3, y = 0.31, width = 0.7, height = 0.35)
dev.off()


pdf("figures/figure-4-new.pdf", width = 3.42, height = 5)
ggdraw(plot_grid(gg_example2, gg_transients, ncol = 1, rel_heights = c(1, 1.5), labels = c("A", "B")))+
  draw_plot(bottom3, x = 0.3, y = 0.31, width = 0.7, height = 0.35)
dev.off()




