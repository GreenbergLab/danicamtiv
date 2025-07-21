
library(data.table)
library(minpack.lm)
library(ggplot2)
library(cowplot)
library(ggtext)


fit_hydrolysis_nls <- function(x){
  nlsLM(rate ~ (`k_+H+k_-H`*atp_spec) / (K1+atp_spec),
      data = x,
      start = list(`k_+H+k_-H` = 70,
                   K1 = 1))
}

fit_hydrolysis_lm <- function(x){
lm(rate ~ 0+atp_spec, data = x[atp_spec < 150])
}

predict_nls_fun <- function(x){
  nd <- data.frame(atp_spec = seq(0, 3000, by = 1))
  data.table(atp_spec = nd$atp_spec,
             predict_y = predict(x, newdata = nd))
}

predict_hydrolysis_lm <- function(x){
data.frame(x = 1:150,
           y = predict(x, newdata = data.frame(atp_spec = 1:150)))
}


## get_kslow_coef <- function(x){
## round(coef(x)[[1]], 1)
## }


get_lm_coef <- function(x){
round(coef(x)[["atp_spec"]], 1)
}


get_hydro_coef <- function(x){
round(coef(x)[["k_+H+k_-H"]], 1)
}


## get_K1 <- function(x){
## round(coef(x)[["K1"]], 1)
## }



data <- fread("stopped-flow/atp-hydrolysis/dani-hydrolysis.csv")
data$id <- factor(data$id, levels = c("dmso", "dani"))


#############################
## FIT KFAST ##
#############################


data_nest <- data[, .(data_nest = list(.SD)), by = .(id)]

## data_nest[, lm_fit := lapply(data_nest, fit_hydrolysis_lm)]
## data_nest[, predict_lm := lapply(lm_fit, predict_hydrolysis_lm)]
## data_nest[, coef_lm_slope := sapply(lm_fit, get_lm_coef)]


## data_lm <- data_nest[, predict_lm[[1]], by = id]

## ggplot()+
##   geom_point(data = data[atp_spec < 200], aes(x = atp_spec, y = rate, color = id))+
##   geom_line(data = data_lm, aes(x = x, y = y, color = id))


## fit a hyperbolic (MM) to all the kfast rates (k1)
data_nest[, nls_fit := lapply(data_nest, fit_hydrolysis_nls)]
data_nest[, predict_nls := lapply(nls_fit, predict_nls_fun)]
data_nest[, coef_nls := sapply(nls_fit, get_hydro_coef)]


data_nls <- data_nest[, predict_nls[[1]], by = id]



colorz <- c("#666666", "#e7298a")

(
gg_hydro <-
ggplot()+
  geom_point(data = data, aes(x = atp_spec, y = rate, color = id),
             size = 2,
             shape = 16,
             alpha = 0.75)+
  geom_line(data = data_nls, aes(x = atp_spec, y = predict_y, color = id), linewidth = 0.9)+
  xlab(~ATP~(mu*M))+
  ylab("<i>k<sub>obs</sub></i> (s<sup>-1</sup>)")+
  scale_color_manual(values = colorz)+
  theme_cowplot(11)+
  theme(
    legend.position = "none",
  plot.title = element_markdown(),
  strip.background = element_rect(fill = "transparent"),
  axis.title.y = element_markdown()
  ## axis.text.x = element_text(size = 14)
  )
  )

png("supplemental-figures/atp-hydrolysis.png", width = 3.42, height = 3, units = "in", res = 300)
gg_hydro
dev.off()


cairo_pdf("supplemental-figures/atp-hydrolysis.pdf", width = 3.43, height = 3)
gg_hydro
dev.off()


## saveRDS(gg_hydro, "figures/gg-hydro.rds")

#################
## day FIT
#####################
coef_data_nest2 <- data[, .(data_nest = list(.SD)), by = .(id, date)]

coef_data_nest2[, nls_fit := lapply(data_nest, fit_hydrolysis_nls)]
coef_data_nest2[, predict_nls := lapply(nls_fit, predict_nls_fun)]
coef_data_nest2[, coef_nls := sapply(nls_fit, get_hydro_coef)]

## estimate second order rate of atp induced dissociation using linear regression
## on low concentration points

## fit a hyperbolic (MM) to all the kfast rates (k1)

coef_sum <- coef_data_nest2[, .(mean = mean(coef_nls),
                                sd = sd(coef_nls)),
                            by = id]

## second_order_ttest <-
  t.test(x = coef_data_nest2[id == "dmso"]$coef_nls,
         y = coef_data_nest2[id == "dani"]$coef_nls
         )


  t.test(x = coef_data_nest[id == "0uM-dani"]$coef_K1,
         y = coef_data_nest[id == "10uM-dani"]$coef_K1
         )

  t.test(x = coef_data_nest[id == "0uM-dani"]$coef_k_plus_2,
         y = coef_data_nest[id == "10uM-dani"]$coef_k_plus_2
         )

  ## shapiro.test(x = coef_data_nest[id == "0uM-dani"]$k_plus_2_over_K1
  ##        ## y = coef_data_nest[id == "10uM-dani"]$k_plus_2_over_K1
  ##        )

  ## shapiro.test(x = coef_data_nest[id == "10uM-dani"]$k_plus_2_over_K1
  ##        )

coef_data_summary <- coef_data_nest[, .(k_plus_2_over_K1 = mean(k_plus_2_over_K1),
                                        sd_k_plus_2_over_K1 = sd(k_plus_2_over_K1),
                                        coef_K1 = mean(coef_K1),
                                        sd_coef_K1 = sd(coef_K1),
                                        coef_k_plus_2 = mean(coef_k_plus_2),
                                        sd_coef_k_plus_2 = sd(coef_k_plus_2)),
                                    by = .(id, id2)]
###
