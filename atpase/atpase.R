library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)

theme_set(theme_cowplot(11))
colorz <- c("#666666", "#e7298a")

dat <- fread("atpase/ATPase.csv")

dat$id <- factor(dat$id, levels = c("DMSO", "Dani"))

dat_sum <- dat[, .(rate = mean(rate),
                   sd_rate = sd(rate)), by = .(id, actin)]

dat_sum_nest <- dat_sum[, .(data = list(.SD)), by = id]

fit_mm <- function(data){
  nls(rate ~ (vmax*actin)/(km+actin),
      data = data,
      start = list(vmax = 5, km = 1))
}

draw_predict <- function(mod){
  nd <- data.frame(actin = seq(0, 80, by = 0.1))
  p <- predict(mod, newdata = nd)
  nd$rate <- p
  return(nd)
}


dat_sum_nest[, mod := lapply(data, fit_mm)]
dat_sum_nest[, predict := lapply(mod, draw_predict)]

predict_df <- dat_sum_nest[, predict[[1]], by = id]

gg_atpase <-
ggplot()+
  geom_errorbar(data = dat_sum,
               aes(actin,
                   ymin = rate-sd_rate,
                   ymax = rate+sd_rate,
                   color = id))+
  geom_point(data = dat_sum,
             aes(actin, rate, color = id),
             alpha = 1,
             size = 1.5)+
  geom_line(data = predict_df,
            aes(actin, rate, color = id),
            linewidth = 0.9)+
  ylab("ATPase Rate (head<sup>-1</sup>&middot;s<sup>-1</sup>)")+
  xlab("[Actin] (&micro;M)")+
  ## ylab(ATPase~Rate~(s^-1))+
  scale_color_manual(values = colorz)+
  theme(
  legend.position = "none",
   axis.title.y = element_markdown(),
   axis.title.x = element_markdown()
  )

saveRDS(gg_atpase, "figures/gg_atpase.rds")


dat_nest <- dat[, .(data = list(.SD)), by = .(id, trial)]

get_coef <- function(mod, n){
  coef(mod)[[n]]
}

dat_nest[, mod := lapply(data, fit_mm)]
dat_nest[, predict := lapply(mod, draw_predict)]
dat_nest[, coef1 := sapply(mod, get_coef, n = 1)]
dat_nest[, coef2 := sapply(mod, get_coef, n = 2)]


t.test(x = dat_nest[id == "DMSO"]$coef1,
 y = dat_nest[id == "Dani"]$coef1
       )


t.test(x = dat_nest[id == "DMSO"]$coef2,
 y = dat_nest[id == "Dani"]$coef2
       )


sd(dat_nest[id=="DMSO"]$coef2)
