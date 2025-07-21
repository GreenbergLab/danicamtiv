library(data.table)
library(minpack.lm)
library(ggplot2)
library(cowplot)
library(ggtext)


folder <- file.path("stopped-flow", "atp-hydrolysis")
files <- list.files(folder,
                    pattern = ".csv",
                    recursive = TRUE,
                    full.names = FALSE)


exclude <- c("2/10uM-dani/10uM-dani_2500uM-atp_0.csv")

files <- files[-c(which(files %in% exclude))]
files <- file.path(folder, files)

read_sf <- function(x){
  d <- fread(x)
  d$path <- x
  return(d)
}

data <- rbindlist(lapply(files, read_sf))

id_split <- c("sf", "hydro", "rep", "id", "filename")

data[, (id_split) := tstrsplit(path, "/", fixed = TRUE)]
data[, c("other", "filenum") := tstrsplit(filename, "atp_", fixed = TRUE)]
data[, filenum := as.factor(sub(".csv", "", filenum))]
data[, id2 := sub("uM-dani", " &micro;M", id)]

data[, `:=`(x_time_s = V1,
            y_signal = V2)]

data <- data[, .(id, id2, rep, filenum, x_time_s, y_signal)]

#nest data
data_nest <- data[, list(spectra = list(.SD)), by = c("id", "id2", "rep", "filenum")]

prep_data <- function(x){
  x <- x[1:5000]
  deadtime <- 0.002
  x <- x[x_time_s >= deadtime]
  ## x[, x_time_s := x_time_s - deadtime]
  x[, y_signal := y_signal-min(y_signal) ]
  x[, y_signal := y_signal/max(y_signal) ]
  return(x)
}

data_nest[, spectra_prepped := lapply(spectra, prep_data)]

fit_mod <- function(x){
  nlsLM(y_signal ~ (a*exp(-k*x_time_s)+c),
       data = x,
       start = list(
         a = -1,
         c = 5 ,
         k = 60),
       control = nls.lm.control(maxiter = 1024)
       )
}

predict_line <- function(mod){
  x <- seq(0, 0.25, by = 0.001)
  y <- predict(mod, newdata = data.frame(x_time_s = x))
  return(
    data.frame(x_time_s = x,
               y_signal = y,
               type = "prediction")
  )
}

get_k <- function(mod){
  coef(mod)[["k"]]
}


data_nest[, mod := lapply(spectra_prepped, fit_mod)]
data_nest[, predict_df := lapply(mod, predict_line)]
data_nest[, coef_k := sapply(mod, get_k)]

data_predict <- data_nest[, predict_df[[1]], by = c("id", "id2", "rep", "filenum")]

data_unnest <- data_nest[, spectra_prepped[[1]], by = c("id","id2", "rep", "filenum")]
data_unnest$type <- "real"


data_together <- rbind(data_unnest, data_predict)
data_together$type <- factor(data_together$type, levels = c("real", "prediction"))



data_summary <- data_nest[, .(avg_k = mean(coef_k),
                              sd_k = sd(coef_k)),
                          by = c("id", "id2")]

setorder(data_summary, id)

colorz <- c("#666666", "#e7298a")


shapiro.test(data_nest[id=="0uM-dani"]$coef_k)
shapiro.test(data_nest[id=="10uM-dani"]$coef_k)

hydro_pval <- t.test(data_nest[id == "0uM-dani"]$coef_k,
       data_nest[id == "10uM-dani"]$coef_k)


gg_hydro <-
ggplot()+
  geom_col(data = data_summary,
           aes(x = id2,
               y = avg_k,
               fill = id2))+
  geom_errorbar(data = data_summary,
                aes(x = id2,
                    ymin = abs(avg_k)-abs(sd_k),
                    ymax = abs(avg_k)+abs(sd_k),
                    color = id2
                    ),
                width = 0.25)+
  geom_jitter(data = data_nest,
              aes(x = id2,
                  y = abs(coef_k)),
              width = 0.1,
              size = 0.5)+
  ylab("<i>k<sub>hydrolysis</sub></i> (s<sup>-1</sup>)")+
  xlab("[Dani]")+
  draw_line(x = c(1, 2),
            y = rep(max(data_nest$coef_k)+10, 2))+
annotate("text", x = 1.5, y = max(data_nest$coef_k)+20, label = paste0("p = ", round(hydro_pval$p.value, 2)), size = 3)+
  ## ggtitle("Rate of Re-Attachment")+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0)))+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  theme_cowplot(11)+
  theme(
  axis.title.y = element_markdown(),
  axis.text.x = element_markdown(),
  legend.position = "none")


## ggplot(data_together[rep == 2 & id =="10uM-dani"])+
##   geom_line(aes(x_time_s,
##                 y_signal,
##                 color = filenum,
##                 linetype = type,
##                 alpha = type))+
##   facet_grid(id2~rep)+
##   ## scale_color_manual(values = colorz)+
##   scale_alpha_manual(values = c(0.7, 1))+
##   scale_linetype_manual(values = c("solid", "dashed"))+
##   theme_cowplot()


gg_transient <-
ggplot(data_together[rep == 1 & filenum == 1])+
  geom_line(aes(x_time_s*1000,
                y_signal,
                color = id2,
                linetype = type,
                alpha = type))+
  ## facet_grid(id2~rep)+
  xlab("Time (ms)")+
  ylab("Normalized Fluorescence")+
  scale_color_manual(values = colorz)+
  scale_alpha_manual(values = c(0.5, 1))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  theme_cowplot(11)+
  theme(
  axis.title.y = element_text(size = 10),
    legend.position = "none" )

## gg_atp_hydro <-
##   ggdraw(gg_transient)+
##   draw_plot(gg_hydro,
##             x = 0.5,
##             y = 0.1,
##             width = 0.5,
##             height = 0.5
##             )

## saveRDS(gg_transient, "figures/gg_hydrolysis.rds")

gg_atp_hydro <-
  plot_grid(gg_transient, gg_hydro, rel_widths = c(1.5, 1), nrow = 1, labels = c("E", "F"))
gg_atp_hydro

saveRDS(gg_atp_hydro, "figures/atp-hydrolysis.rds")

ggsave("~/Downloads/atp-hydrolysis-dani.png", bg = "white")

