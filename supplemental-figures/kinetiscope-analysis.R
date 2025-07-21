library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
library(RColorBrewer)

dat <- fread("supplemental-figures/kinetiscope-concentration-time.txt", skip = 8)

col_names <- c("dp", "time_s", "AM", "AM.Dani", "M", "M.Dani", "OFF", "OFF.Dani")

colnames(dat) <- col_names

dat_melt <- melt(data = dat, measure.vars = 3:8, variable.name = "population", value.name = "concentration")

dark2 <- brewer.pal(8, "Dark2")

colorz1 <- c(dark2[1:3], dark2[5], dark2[8], dark2[4])
colorz2 <- c(dark2[8], dark2[4])


gg1 <-
ggplot(dat_melt)+
  geom_line(aes(x = time_s,
                y = concentration,
                color = population))+
  xlab("Time (s)")+
  ylab("Conc. (mole/liter)")+
  scale_color_manual(values = colorz1)+
  theme_cowplot(12)

off_data <- dat_melt[grepl("OFF", population)]



off_nest <- off_data[, list(
  data_nest = list(
    data.table(
      time_s,
      concentration)
  )
),
by = .(population)]


fit_exp_fun <- function(x){
nls(concentration ~ a*exp(-k*time_s)+c,
    data = x,
    start = list(k = 10,
                 a = 1,
                 c = 1)
    )
}

predict_fun <- function(x){
  nd <- data.frame(time_s = seq(0, 1.20, by = 1/5000))
  data.table(time_s = nd$time_s,
            concentration = predict(x, newdata = nd))
}


get_coef <- function(x){
round(coef(x)[["k"]], 1)
}

off_nest[, exp_fit := lapply(data_nest, fit_exp_fun)]
off_nest[, predict_df := lapply(exp_fit, predict_fun)]
off_nest[, coef_k := sapply(exp_fit, get_coef)]

off_predict_df <- off_nest[, predict_df[[1]], by = population]
off_predict_df$id <- "fit"

off_data$id <- "real"

off_plot_data <- rbind(off_data, off_predict_df, fill= TRUE)

gg2 <-
ggplot(off_plot_data)+
  geom_line(aes(x = time_s,
                y = concentration,
                color = population,
                linetype = id,
                alpha = id))+
  annotate("richtext", x = 0.6 , y = 0.5, label = paste0("<br> ",
                                                   off_nest$population[[1]], " = ", off_nest$coef_k[[1]], " s<sup>-1</sup>",
                                                   "<br>",
                                                   off_nest$population[[2]], " = ", off_nest$coef_k[[2]], " s<sup>-1</sup>"
                                                   ),
           label.color = "transparent")+
  scale_linetype_manual(values = c("dashed", "solid"))+
  scale_alpha_manual(values = c(1, 0.5))+
  scale_color_manual(values = colorz2)+
  ylab("Conc. (mole/liter)")+
  xlab("Time (s)")+
  theme_cowplot(12)

scheme1 <- magick::image_read_pdf("supplemental-figures/kin-scheme1.pdf")

print(scheme1)


png("supplemental-figures/ggkinetiscope.png", width = 5, height = 5, units = "in", res = 500)
plot_grid(gg1, gg2, labels = c("A", "B"), nrow = 2)
dev.off()
