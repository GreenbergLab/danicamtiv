library(readxl)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
library(magick)
library(grid)
library(gridExtra)
library(minpack.lm)
library(RColorBrewer)

colorz <- brewer.pal(8, "Dark2")[c(8, 4)]
basesize <- 11

adp_data <- fread("stopped-flow/adp-affinity/adp-affinity-dani.csv")
adp_data[, adp_conc := adp_conc_uM/2]
adp_data <- adp_data[adp_conc <= 25, .(k = mean(k)), by = .(date, id, adp_conc)]


adp_data$id <- ifelse(adp_data$id == "control", "0 &micro;M", adp_data$id)
adp_data$id <- ifelse(adp_data$id == "10uM-dani", "10 &micro;M", adp_data$id)
## adp_data$id <- factor(adp_data$id, levels = c("control", "10uM-dani"))


fit_adp <- function(x){
  ## k0 <- mean(x[adp_conc == 0]$k)
  nlsLM(k ~ k0 / (1 + (adp_conc/K2dK1d)),
      data = x,
      start = list(k0 = 100, K2dK1d = 10))
}

predict_adp <- function(x){
  nd <- data.frame(adp_conc = seq(0, 25, by = 0.1))
  data.table(adp_conc = nd$adp,
             predict_y = predict(x, newdata = nd))
}


get_adp_coef <- function(x){
round(coef(x)[["K2dK1d"]], 1)
}


adp_data_nest <- adp_data[, .(data_nest = list(.SD)), by = id]
adp_data_nest[, adp_fit := lapply(data_nest, fit_adp)]
adp_data_nest[, predict_adp:= lapply(adp_fit, predict_adp)]
adp_data_nest[, coef_adp := sapply(adp_fit, get_adp_coef)]

adp_data_nest2 <- adp_data[, .(data_nest = list(.SD)), by = .(id, date)]
adp_data_nest2[, adp_fit := lapply(data_nest, fit_adp)]
adp_data_nest2[, predict_adp:= lapply(adp_fit, predict_adp)]
adp_data_nest2[, coef_adp := sapply(adp_fit, get_adp_coef)]

adp_sum <- adp_data_nest2[, .(adp_mean = mean(coef_adp),
                             adp_sd = sd(coef_adp)),
                         by = id]


adp_affinity_pval <- t.test(adp_data_nest2[id == "0 &micro;M"]$coef_adp,
                            adp_data_nest2[id == "10 &micro;M"]$coef_adp)


gg_adp_affinity1 <-
ggplot()+
  geom_col(data = adp_sum, aes(id, adp_mean, fill = id))+
  geom_errorbar(data = adp_sum, aes(id,
                                    ymax = adp_mean+adp_sd,
                                    ymin = adp_mean-adp_sd,
                                    color = id),
                width = 0.25)+
  geom_jitter(data = adp_data_nest2, aes(id, coef_adp), width = 0.25)+
  draw_line(c(1, 2), c(10, 10))+
  annotate("text", x = 1.5, y = 10.5, label = paste0("p = ", round(adp_affinity_pval$p.value, 2)), size = basesize/.pt)+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.1)))+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  xlab("Dani")+
ylab("<i>K</i><sub>5</sub>' (&micro;M)")+
ylab("ADP Affinity (&micro;M)")+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
    axis.title.y = element_markdown(),
    axis.text.x = element_markdown(size = 9),
    axis.title.x = element_text(color = "white")
  )
gg_adp_affinity1



adp_predict_df <- adp_data_nest[, predict_adp[[1]], by = .(id)]
adp_predict_df2 <- adp_data_nest2[, predict_adp[[1]], by = .(id, date)]
## adp_predict_df$id <- factor(adp_predict_df$id, levels = c("control", "10uM-dani"))

(adp <-
ggplot()+
  geom_point(data = adp_data[date == "12/8/2023"],
             aes(x = adp_conc,
                 y = k,
                 color = id,
                 shape = date),
             size = 3,
             ## shape = 16,
             alpha = 0.6
             )+
  geom_line(data = adp_predict_df2[date == "12/8/2023"],
            aes(x = adp_conc,
                y = predict_y,
                color = id),
            size = 0.9)+
  ## facet_wrap(~date)+
  ## ggtitle("K<sub>2D</sub><sup>'</sup>K<sub>1D</sub><sup>'</sup>")+
  ## ggtitle("ADP Affinity")+
  scale_y_continuous(breaks = seq(0, 150, 50))+
  coord_cartesian(ylim = c(0, NA))+
  xlab(~ADP~(mu*M))+
  ylab("<i>k<sub>obs</sub></i> (s<sup>-1</sup>)")+
  scale_color_manual(values = colorz)+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
  plot.title = element_markdown(),
  axis.title.y = element_markdown()
  )
  )

## adp <-
##   ggdraw(adp)+
##   draw_plot(gg_adp_affinity1, 0.5, 0.25, 0.45, 0.7)
## adp


## saveRDS(adp, "stopped-flow/k2dk1d.rds")
(gg_adp_affinity <- plot_grid(adp, gg_adp_affinity1, rel_widths = c(1, 0.5), align = "h"))
saveRDS(gg_adp_affinity, "stopped-flow/gg_adp_affinity.rds")

## ggsave("~/Downloads/k2dk1d.png")
  ## annotate("richtext",
  ##          x = 0,
  ##          y = 1200,
  ##          label = paste0("<span style = 'color: black;'> pH 7.0 = ", round(sf_data_nest$coef_k_2[[1]], 0), " s <sup>-1</sup></span><br>",
  ##                        "<span style = 'color: red;'> pH 6.5 = ", round(sf_data_nest$coef_k_2[[2]], 0), " s <sup>-1</sup></span><br>"),
  ##          color = "transparent",
  ##          fill = "transparent",
  ##          hjust = 0,
  ##          vjust = 1
  ##          )+
  ## annotate("text",
  ##          x = Inf,
  ##          y = -Inf,
  ##          label = list(bquote(atop(k[+2]==.(k2),
  ##                                   K[1]==.(K1)))),
  ##          parse = TRUE,
  ##          vjust = -0.5,
  ##          hjust = 1)+
  ## coord_cartesian(ylim = c(0, 650))+
  ## ggtitle("k<sub>+2</sub>'")+
  ## xlab(~ATP~(mu*M))+
  ## ylab(~k[fast]~(s^-1))+
  ## scale_color_manual(values = colorz)+
  ## theme_cowplot(basesize)+
  ## theme(
  ##   legend.position = "none",
  ## plot.title = element_markdown(),
  ## axis.text.x = element_text(size = 7)
  ## )
  ## )
