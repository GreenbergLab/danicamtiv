library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)

colorz <- c("#666666", "#e7298a")
basesize <- 11

adp_file_paths <- list(
    `control` = list.files("stopped-flow/adp-release/adp-release-raw-data/control",
                             recursive = TRUE,
                             full.names = TRUE),
   `10uM-dani` = list.files("stopped-flow/adp-release/adp-release-raw-data/10uM-dani",
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
        plot.title = element_markdown()
        )
)

all_adp_rates <- fread("stopped-flow/adp-release/dani-adp-release.csv")
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
  geom_jitter(data = all_adp_rates, aes(id, k, shape = date), width = 0.1)+
  scale_color_manual(values = colorz)+
  scale_fill_manual(values = colorz)+
  draw_line(c(1, 2), c(85, 85))+
  annotate("text", x = 1.5, y = 90, label = paste0("p = ", round(adp_release_pval$p.val, 2)), size = basesize/.pt)+
  xlab("Dani")+
  ylab("ADP Release Rate (s<sup>-1</sup>)")+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.1)))+
  theme_cowplot(basesize)+
  theme(legend.position = "none", plot.title = element_markdown(), axis.title.y = element_markdown(),
        axis.text.x = element_markdown(),
        axis.title.x = element_text(color = "white"))
gg_adp2

gg_adp_all <- plot_grid(gg_adp, gg_adp2, rel_widths = c(1.5, 1), align = "h")

saveRDS(gg_adp_all, "stopped-flow/gg_adp.rds")
