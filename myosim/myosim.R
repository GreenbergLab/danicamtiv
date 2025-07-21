library(data.table)

## define a protocol
sarcomere_length_nm <- 2.1*1000

# sarcomere shorten distance
ssd <- sarcomere_length_nm*0.2

hz <- 2000

time_step <- 1/hz

#Mode -2 = length control
#mode -1 = length control with potentnial slack

protocol1 <-
data.table(dt = c(rep(time_step, 1*hz), time_step, rep(time_step, 40), time_step, rep(time_step, 1*hz)),
           dhsl = c(rep(0, 1*hz), -ssd, rep(0, 40), ssd, rep(0, 1*hz)),
           Mode = c(rep(-2, 1*hz), -1, rep(-1,40 ), -1, rep(-2, 1*hz)),
           pCa = 6.1)

protocol2 <-
data.table(dt = c(rep(time_step, 1*hz), time_step, rep(time_step, 40), time_step, rep(time_step, 1*hz)),
           dhsl = c(rep(0, 1*hz), -ssd, rep(0, 40), ssd, rep(0, 1*hz)),
           Mode = c(rep(-2, 1*hz), -1, rep(-1,40 ), -1, rep(-1, 1*hz)),
           pCa = 6.1)

## protocol <-
##   data.table(dt = c(8, time_step, 40/1000, time_step, 8),
##            dhsl = c(0, -ssd, 0, ssd, 0),
##            Mode = c(-2, -1, -1, -1, -2),
##            pCa = 6.1)


fwrite(protocol2, "myosim/protocol2.txt", sep = " ")



#####################
## read sim output ##
#####################
library(data.table)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

myosim_files <- list.files("~/MATMyoSim/danicamtiv/ktr/protocol-1/k_coop",
                           pattern = "myosim-output.csv",
                           recursive = TRUE,
                           full.names = TRUE)

fread2 <- function(x){
  y <- fread(x)
  y$filename <- x
  return(y)
}

myosim <- rbindlist(lapply(myosim_files, fread2))

myosim[, c("blank", "home", "brent", "mat", "dan", "ktr", "protocol", "coop", "id", "file") := tstrsplit(filename, "/", fixed = TRUE)]
## myosim[, c("xps", "cf", "kcb", "k3", "k1") := tstrsplit(id, "_", fixed = TRUE)]

## myosim <- myosim[protocol == "protocol-1"]

myosim[, norm_muscle_force := (muscle_force-min(muscle_force))/(max(muscle_force)-min(muscle_force)), by = .(id)]

## myosim_sub <- myosim[protocol == "protocol-1" & cf == "0.5-CF" & time_s >= 1.055]
## myosim_sub[, norm_muscle_force := (muscle_force-min(muscle_force))/(max(muscle_force)-min(muscle_force)), by = .(id, xps, cf)]

ggplot(myosim)+
  geom_line(aes(time_s,
                 norm_muscle_force,
                 color = id),
            linewidth = 1,
            alpha = 0.9)+
  ## scale_colour_manual(values = c("#666666", "#e7298a"))+
  ylab("Normalized Force")+
  xlab("Time (s)")+
  theme(
  ## legend.position  = c(0.5, 0.5)
  )




## ggplot(myosim_sub)+
##   geom_line(aes(time_s,
##                  muscle_force,
##                  color = id))


p1 <-
ggplot(myosim)+
  geom_line(aes(time_s, norm_muscle_force, color = id), linewidth = 1)+
  xlab("Time (s)")+
  ylab("Normalized Force")+
  ## scale_colour_manual(values = c("#666666", "#e7298a"))+
  ## facet_wrap(~protocol)+
  coord_cartesian(xlim = c(1.0, 1.5))


p2 <-
ggplot(myosim)+
  geom_line(aes(time_s, muscle_length, color = id), linewidth = 1)+
  xlab("Time (s)")+
  ylab("Muscle Length")+
  ## scale_colour_manual(values = c("#666666", "#e7298a"))+
  ## facet_wrap(~protocol)+
  coord_cartesian(xlim = c(1, 1.5))


plot_grid(p1, p2, nrow = 2)


ggplot(myosim)+
  geom_line(aes(time_s, muscle_force, color = id), linewidth = 1)+
  xlab("Time (s)")+
  ylab("Normalized Force")+
  ## scale_colour_manual(values = c("#666666", "#e7298a"))+
  ## facet_wrap(~protocol)+
  coord_cartesian(xlim = c(1.0, 1.5))

########
# tension pca
#########

library(data.table)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

myosim_files <- list.files("~/MATMyoSim/danicamtiv/tension_pCa/test2",
                           pattern = "pca-curve.csv",
                           recursive = TRUE,
                           full.names = TRUE)


fit_files <- list.files("~/MATMyoSim/danicamtiv/tension_pCa/test2",
                           pattern = "curve-fit.csv",
                           recursive = TRUE,
                           full.names = TRUE)


fread2 <- function(x){
  y <- fread(x)
  y$filename <- x
  return(y)
}

myosim <- rbindlist(lapply(myosim_files, fread2))
fit <- rbindlist(lapply(fit_files, fread2))

myosim[, c("blank", "home", "brent", "mat", "dan", "tension", "test", "id", "file") := tstrsplit(filename, "/", fixed = TRUE)]
fit[, c("blank", "home", "brent", "mat", "dan", "tension", "test", "id", "file") := tstrsplit(filename, "/", fixed = TRUE)]
## myosim[, c("xps", "cf") := tstrsplit(id, "_", fixed = TRUE)]

myosim_nest <- myosim[, .(data = list(.SD)) , by = id ]
fit_nest <- fit[, .(fit = list(.SD)) , by = id ]

myo_nest <- myosim_nest[fit_nest, on = "id"]

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
