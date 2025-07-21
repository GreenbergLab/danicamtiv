library(data.table)
library(ggplot2)
library(cowplot)
library(drc)


fit_hill <- function(data, response){
  if(data$avg_percent_moving[[1]] > 0.1){
  drm(paste0(response, "~ pca"),
      data = data,
      fct = LL.4(names = c("hillslope", "min", "vmax", "ec50")),
      logDose = 10)
  } else {
  drm(paste0(response, "~ pca"),
      data = data,
      fct = LL.3(names = c("hillslope", "vmax", "ec50")),
      logDose = 10)

  }
}

predict_hill <- function(x){
  ## browser()
  dummy_x <- expand.grid(pca=exp(seq(log(4), log(9), length=1000)))
  predict_y <- predict(x, newdata = dummy_x, interval = "confidence")
  df <- as.data.table(predict_y)
  df$pca <- dummy_x$pca
  ## df <- data.table(dummy_x,
  ## y = predict_y)
  return(df)
}









dat <- fread("motility/regulated/percent-moving-hand-tracked.csv")
dat <- na.omit(dat)


dat <- dat[date != "2023-10-12" | pca != 9]

dat[, percent_moving := number_moving/total_filaments]

dat_sum <- dat[, .(avg_percent_moving = mean(percent_moving),
                   sd_percent_moving = sd(percent_moving)),
               by = .(pca, conditions)]


pca_nest <- dat_sum[, .(data = list(.SD)), by = conditions]
pca_nest[, hill_fit := lapply(data, fit_hill, response = "avg_percent_moving")]
pca_nest[, predict_fit := lapply(hill_fit, predict_hill)]
pca_nest[, mod_table := lapply(hill_fit, broom::tidy)]

pca_lines <- pca_nest[, predict_fit[[1]], by = conditions]

gg <-
ggplot()+
  geom_point(data = dat_sum,
             aes(pca,
                 avg_percent_moving,
                 color = conditions),
             size  = 1.5)+
  ## geom_errorbar(aes(pca,
  ##                   ymin = avg_percent_moving-sd_percent_moving,
  ##                   ymax = avg_percent_moving+sd_percent_moving,
  ##               color = conditions))+
  geom_line(data = pca_lines,
            aes(pca, Prediction, color = conditions),
            linewidth = 0.8)+
  scale_x_reverse()+
  ylab("Percent Moving (%)")+
  xlab("pCa")+
  scale_color_manual(values = c("#666666", "#e7298a"))+
  theme_cowplot(11)+
  theme(
    ## axis.title.y = element_markdown(),
    legend.position = "none"
    )


## ggsave("supplemental-figures/regulated-percent-moving.png", bg = "white")

png("supplemental-figures/regulated-motility-percent-moving.png", bg = "white", width = 3, height = 2.5, units = "in", res = 300)
gg
dev.off()


pdf("supplemental-figures/regulated-motility-percent-moving.pdf", bg = "white", width = 3, height = 2.5)
gg
dev.off()
