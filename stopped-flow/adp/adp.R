library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
library(minpack.lm)
library(extrafont)

loadfonts()

colorz <- c("#666666", "#e7298a")
basesize <- 11

theme_set(theme_cowplot(11))

adp_file_paths <- list(
    `control` = list.files("stopped-flow/adp/adp-release/adp-release-raw-data/control",
                             recursive = TRUE,
                             full.names = TRUE),
   `10uM-dani` = list.files("stopped-flow/adp/adp-release/adp-release-raw-data/10uM-dani",
                             recursive = TRUE,
                             full.names = TRUE)
  )

adp_data_list <- list()
num <- 1
for(myo in seq_along(adp_file_paths)){
  id <- names(adp_file_paths)[[myo]]
  for(csv in seq_along(adp_file_paths[[myo]])){
   adp_dat <- fread(adp_file_paths[[myo]][[csv]])#[1:5000,]
   adp_dat[, `:=`(time = V1,
                  dead_time = V1 - 0.003,
                  normalized_signal_0 = (V2-head(V2, 1)))]
   adp_dat[, normalized_signal_1 := normalized_signal_0/tail(normalized_signal_0, 100)]
   adp_dat$id <- id
   adp_dat$id_rep <- paste0(id, "-", csv)
   adp_dat$rep <- csv
   adp_dat$type <- "raw"
   adp_data_list[[num]] <- adp_dat
   num <- num + 1
   }
  }

adp_df <- rbindlist(adp_data_list)

## adp_df$id <- factor(adp_df$id, levels = c("control", "10uM-dani"))

adp_avg <- adp_df[, .(avg_signal = mean(V2),
                      type = "avg"), by = .(id, time, dead_time)]
adp_avg[, avg_signal_normalized_0 := (avg_signal-head(avg_signal, 1)), by = .(id)]
adp_avg[, avg_signal_normalized_1 := avg_signal_normalized_0/tail(avg_signal_normalized_0, 100), by = .(id)]


adp_avg_nest <- adp_avg[, list(
                           data_nest = list(
                                         data.table(
                                           dead_time,
                                           avg_signal_normalized_1)
                           )
            ),
  by = .(id)]


fit_exp_fun <- function(x){
nls(avg_signal_normalized_1 ~ a*exp(-k*dead_time)+c,
    data = x[dead_time>=0&dead_time<=0.25],
    start = list(k = 100,
                 c = 1,
                 a = -1)
    )
}

predict_fun <- function(x){
  nd <- data.frame(dead_time = seq(0, 0.20, by = 1/5000))
  data.table(dead_time = nd$dead_time,
            predict_y = predict(x, newdata = nd))
}


get_coef <- function(x){
round(coef(x)[["k"]], 0)
}

adp_avg_nest[, exp_fit := lapply(data_nest, fit_exp_fun)]
adp_avg_nest[, predict_df := lapply(exp_fit, predict_fun)]
adp_avg_nest[, coef_k := sapply(exp_fit, get_coef)]

adp_predict_df <- adp_avg_nest[, predict_df[[1]], by = id]

(gg_adp <-
  ggplot()+
  geom_line(data = adp_avg[dead_time>=0&dead_time<=0.20],
            aes(x = dead_time,
                y = avg_signal_normalized_1,
                color = id),
            alpha = 0.5,
            linewidth = 1)+
  geom_line(data = adp_predict_df,
            aes(dead_time, predict_y,
                color = id,
                linetype = id),
            alpha = 1)+
  ## annotate("richtext",
  ##         x = Inf, y = -Inf, color = "transparent", fill = "transparent",
  ##         label = paste0("<span style = 'color: black; float: left;'> pH 7.0 = ", adp_avg_nest$coef_k[[1]], " s<sup>-1</sup></span> <br>",
  ##                        "<span style = 'color: red; float: left;'> pH 6.5 = ", adp_avg_nest$coef_k[[2]], " s<sup>-1</sup></span> <br>  "),
  ##         hjust = 1,
  ##         vjust = 0.25
  ##         )+
  scale_color_manual(values = colorz, name = "")+
  xlab("Time (s)")+
  ylab("Normalized Flourescence")+
   scale_linetype_manual(values = c("solid", "dashed"))+
  ## ggtitle("ADP Release")+
  theme_cowplot(basesize)+
  theme(legend.position = "none",
        plot.title = element_markdown(),
        axis.title.y = element_text(size = 10)
        )
)

all_adp_rates <- fread("stopped-flow/adp/adp-release/dani-adp-release.csv")
all_adp_rates$id <- ifelse(all_adp_rates$id == "control", "0 &micro;M", all_adp_rates$id)
all_adp_rates$id <- ifelse(all_adp_rates$id == "10uM-dani", "10 &micro;M", all_adp_rates$id)
## all_adp_rates$id <- factor(all_adp_rates$id, levels = c("control", "10uM-dani"))
adp_summary <- all_adp_rates[, .(avg_k = mean(k),
                                 sd_k = sd(k)), by = id]

all_adp_con <- all_adp_rates[id == "0 &micro;M"]
all_adp_dani <- all_adp_rates[id == "10 &micro;M"]

adp_release_pval <- t.test(all_adp_con$k, all_adp_dani$k)

gg_adp2 <-
ggplot()+
  geom_errorbar(data = adp_summary, aes(id, ymin = avg_k-sd_k, ymax = avg_k+sd_k, color = id), width = 0.2)+
  geom_col(data = adp_summary, aes(id, avg_k, fill = id))+
  geom_jitter(data = all_adp_rates, aes(id, k, shape = date), width = 0.1, size = 1)+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  draw_line(c(1, 2), c(85, 85))+
  annotate("text", x = 1.5, y = 90, label = paste0("p = ", round(adp_release_pval$p.val, 2)), size = basesize/.pt)+
  xlab("[Dani]")+
  ylab("ADP Release Rate (s<sup>-1</sup>)")+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.1)))+
  theme_cowplot(basesize)+
  theme(legend.position = "none", plot.title = element_markdown(), axis.title.y = element_markdown(),
        axis.text.x = element_markdown(),
        axis.title.x = element_text())
gg_adp2

gg_adp_all <- plot_grid(gg_adp, gg_adp2, rel_widths = c(1, 0.5), align = "h")

## saveRDS(gg_adp_all, "stopped-flow/gg_adp.rds")





adp_data <- fread("stopped-flow/adp/adp-affinity/adp-affinity-dani.csv")
adp_data[, adp_conc := adp_conc_uM/2]
adp_data <- adp_data[
  adp_conc <= 80 & adp_conc !=0
,
  .(k = mean(k)), by = .(date, id, adp_conc)]


adp_data$id <- ifelse(adp_data$id == "control", "0 &micro;M", adp_data$id)
adp_data$id <- ifelse(adp_data$id == "10uM-dani", "10 &micro;M", adp_data$id)
## adp_data$id <- factor(adp_data$id, levels = c("control", "10uM-dani"))

adp_data_sum <- adp_data[, .(k = mean(k),
                             sd_k  = sd(k)),
                         by = .(id, adp_conc)]

adp_data_sum <- adp_data_sum[adp_conc != 0 ]

fit_adp <- function(x){
  ## k0 <- mean(x[adp_conc == 0]$k)
  nlsLM(k ~ k0 / (1 + (adp_conc/K2dK1d)),
      data = x,
      start = list(k0 = 100, K2dK1d = 10))
}

predict_adp <- function(x){
  nd <- data.frame(adp_conc = seq(0, 50, by = 0.1))
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
  geom_jitter(data = adp_data_nest2, aes(id, coef_adp), width = 0.25, size = 1)+
  draw_line(c(1, 2), c(25, 25))+
  annotate("text", x = 1.5, y = 27, label = paste0("p = ", round(adp_affinity_pval$p.value, 2)), size = basesize/.pt)+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.1)))+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  xlab("[Dani]")+
ylab("<i>K</i><sub>5</sub>' (&micro;M)")+
ylab("ADP Affinity (&micro;M)")+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
    axis.title.y = element_markdown(),
    axis.text.x = element_markdown(size = 9),
    axis.title.x = element_text()
  )
gg_adp_affinity1



adp_predict_df <- adp_data_nest[, predict_adp[[1]], by = .(id)]
adp_predict_df2 <- adp_data_nest2[, predict_adp[[1]], by = .(id, date)]
## adp_predict_df$id <- factor(adp_predict_df$id, levels = c("control", "10uM-dani"))

(gg_adp_affinity2 <-
ggplot()+
  geom_errorbar(data = adp_data_sum,
                aes(x = adp_conc,
                    ymin = k-sd_k,
                    ymax = k+sd_k,
                    color  = id
                    ),
                    width = 0.5)+
  geom_point(data = adp_data_sum,
             aes(x = adp_conc,
                 y = k,
                 color = id,
                 ## shape = date
                 ),
             size = 2,
             ## shape = 16,
             alpha = 0.75
             )+
  geom_line(data = adp_predict_df,
            aes(x = adp_conc,
                y = predict_y,
                color = id),
            linewidth = 0.9)+
  ## annotate("richtext", x = Inf, y = Inf, label = paste0("M&middot;D &#x21cc; M+D"),
  ##          hjust = 1, vjust     = 1,
  ##          label.color = NA)+
  ## facet_wrap(~date)+
  ## ggtitle("K<sub>2D</sub><sup>'</sup>K<sub>1D</sub><sup>'</sup>")+
  ## ggtitle("ADP Affinity")+
  ## scale_y_continuous(breaks = seq(0, 150, 50))+
  coord_cartesian(ylim = c(0, NA))+
  xlab(~ADP~(mu*M))+
  ylab("<i>k<sub>obs</sub></i> (s<sup>-1</sup>)")+
  scale_color_manual(values = colorz)+
  ## theme_cowplot(basesize)+
  theme(
    legend.position = "none",
  plot.title = element_markdown(),
  axis.title.y = element_markdown()
  )
  )


## ggsave("~/Downloads/adp-affinity.pdf", device = cairo_pdf)

## pdf("figures/adp-affinity-with-error.pdf", height = 4, width = 4)
## gg_adp_affinity2
## dev.off()

saveRDS(gg_adp_affinity2, "figures/gg_adp_affinity.rds")
saveRDS(gg_adp, "figures/gg_adp.rds")

## (gg_adp_affinity2 <-
## ggplot()+
##   geom_point(data = adp_data[date == "12/8/2023"],
##              aes(x = adp_conc,
##                  y = k,
##                  color = id,
##                  shape = date),
##              size = 3,
##              ## shape = 16,
##              alpha = 0.6
##              )+
##   geom_line(data = adp_predict_df2[date == "12/8/2023"],
##             aes(x = adp_conc,
##                 y = predict_y,
##                 color = id),
##             linewidth = 0.9)+
##   ## facet_wrap(~date)+
##   ## ggtitle("K<sub>2D</sub><sup>'</sup>K<sub>1D</sub><sup>'</sup>")+
##   ## ggtitle("ADP Affinity")+
##   scale_y_continuous(breaks = seq(0, 150, 50))+
##   coord_cartesian(ylim = c(0, NA))+
##   xlab(~ADP~(mu*M))+
##   ylab("<i>k<sub>obs</sub></i> (s<sup>-1</sup>)")+
##   scale_color_manual(values = colorz)+
##   theme_cowplot(basesize)+
##   theme(
##     legend.position = "none",
##   plot.title = element_markdown(),
##   axis.title.y = element_markdown()
##   )
##   )

## adp <-
##   ggdraw(adp)+
##   draw_plot(gg_adp_affinity1, 0.5, 0.25, 0.45, 0.7)
## adp


## saveRDS(adp, "stopped-flow/k2dk1d.rds")
## (gg_adp_affinity <- plot_grid(adp, gg_adp_affinity1, rel_widths = c(1, 0.5)))
## saveRDS(gg_adp_affinity, "stopped-flow/gg_adp_affinity.rds")

## saveRDS(gg_adp, "figures/gg_adp.rds")
## saveRDS(gg_adp_affinity2, "figures/gg_adp_affinity.rds")

gg_adp_right <- plot_grid(gg_adp2, gg_adp_affinity1, nrow = 2, align = "v", labels = c("B", "D"))
gg_adp_left <- plot_grid(gg_adp, gg_adp_affinity2, nrow = 2, align = "v", labels = c("A", "C"))


png("figures/adp.png", width = 6, height = 5.25, units = "in", res = 300)
plot_grid(gg_adp_left,
          gg_adp_right,
          nrow = 1,
          align = "h",
          rel_widths = c(2, 1))
dev.off()


pdf("figures/for-the-boss/adp.pdf", width = 6, height = 5.25)
plot_grid(gg_adp_left,
          gg_adp_right,
          nrow = 1,
          align = "h",
          rel_widths = c(2, 1))
dev.off()


svg("figures/for-the-boss/adp.svg", width = 6, height = 5.25)
plot_grid(gg_adp_left,
          gg_adp_right,
          nrow = 1,
          align = "h",
          rel_widths = c(2, 1))
dev.off()


###########
##error prop#
#       sd / mean        sd2  mean2

sqrt(((3.3/73.2)^2) + ((4.6/18.0)^2))
sqrt(((3.3/73.1)^2) + ((4.1/18.6)^2))

BSDA::tsum.test(mean.x = 4.1, s.x = 0.3, n.x = 3,
                mean.y = 3.9, s.y = 0.2, n.y = 3)
