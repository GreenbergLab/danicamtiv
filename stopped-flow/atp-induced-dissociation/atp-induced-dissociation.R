library(data.table)
library(minpack.lm)
library(ggplot2)
library(cowplot)
library(ggtext)

fit_lmod <- function(x){
  lm(coef_kfast ~ atp_conc, data = x)
}

fit_min_mod <- function(x){
  ## browser()
min_mod <- nlsLM(y_signal ~ (a*exp(-k*x_time_s)+c),
                 data = x,
                 start = list(
                   a = -1,
                   c = 10 ,
                   k = 50),
                 control = nls.lm.control(maxiter = 1024)
                 )
  }

fit_2exp <- function(x, a1, c_est){
  ## browser()
  tryfit <- 1
  while(tryfit > 0){
  if(tryfit == 20){
    return(NA)
  }
    if(tryfit == 1){
      coefk1 <- 500
      coefk2 <- 50
      coefc <- c_est
      coefa2 <- -0.2
      coefa1 <- -2
    } else {
      coefk1 <- runif(1, 100, 1000)
      coefk2 <- runif(1, 1, 100)
      coefa1 <- runif(1, -10, -0.10)
      coefa2 <- runif(1, -1, -0.001)
      coefc <- runif(1, -10, 10)
    }
    ## while(try_num < try_times){
      mod <- try(nlsLM(y_signal ~ (a11*exp(-k1*x_time_s))+(a2*exp(-k2*x_time_s))+c,
                       data = x,
                       start = list(
                         a11 = coefa1,
                         a2 = coefa2,
                         c = coefc ,
                         k1 = coefk1,
                         k2 = coefk2),
                       ## control = nls.lm.control(maxiter = 1024)
                       control = nls.lm.control(maxiter = 1024)
                       ))

      if(class(mod) == "try-error"){
        tryfit <- tryfit + 1
      } else {
      coefk1 <- coef(mod)[["k1"]]
      coefk2 <- coef(mod)[["k2"]]
      coefc <- coef(mod)[["c"]]
      coefa2 <- coef(mod)[["a2"]]

      mod <- nlsLM(y_signal ~ (a1*exp(-k1*x_time_s))+(a2*exp(-k2*x_time_s))+c,
                       data = x,
                       start = list(
                         a2 = coefa2,
                         c = coefc ,
                         k1 = coefk1,
                         k2 = coefk2),
                       ## control = nls.lm.control(maxiter = 1024)
                       control = nls.lm.control(maxiter = 1024)
                       )
        tryfit <- -1
      }
    }
return(mod)
}

predict_line <- function(mod, deadtime){
  x <- seq(deadtime, 0.25, by = 0.001)
  y <- predict(mod, newdata = data.frame(x_time_s = x))
  return(
    data.table(x_time_s = x,
               y_signal = y,
               type = "prediction")
  )
}

get_coef <- function(mod, x){
  coef(mod)[[x]]
}

read_sf <- function(x){
  ## on this day, settings were set to collect 4000 points, other day 5000 points
  if(grepl("2023-09-25", x)){
    nrows = 4000
  } else {
    nrows = 5000
  }
  d <- fread(x, nrows = nrows)
  d$path <- x
  return(d)
}


prep_data <- function(x, deadtime){
  x <- x[x_time_s >= deadtime]
  ## x[, x_time_s := x_time_s - deadtime]
  ## x[, y_signal := (y_signal-min(y_signal))/(max(y_signal)-min(y_signal))]
  ## x[, y_signal := y_signal/max(y_signal) ]
  ## x[, y_signal := y_signal/max(y_signal) ]
  return(x)
}

fit_kfast_nls <- function(x){
  nlsLM(coef_kfast ~ (`k_+2`*atp_conc) / (K1+atp_conc),
      data = x,
      start = list(`k_+2` = 1000,
                   K1 = 1/300))
}

fit_kfast_lm <- function(x){
lm(coef_kfast ~ 0+atp_conc, data = x[atp_conc < 200])
}

predict_kfast_fun <- function(x){
  nd <- data.frame(atp_conc = seq(0, 6500, by = 1))
  data.table(atp_conc = nd$atp_conc,
             predict_y = predict(x, newdata = nd))
}

predict_kfast_lm <- function(x){
data.frame(x = 1:200,
           y = predict(x, newdata = data.frame(atp_conc = 1:200)))
}


fit_atp_induced_dissociation <- function(data_nest){
  ## prep data
  ## browser
  data <- copy(data_nest)
  data <- data[, spectra_prepped := lapply(spectra, prep_data, deadtime = 0.003)]
  setorder(data, atp)
  ## calculate a single exp for the minimum atp concentration
  data_min <- data[1]
  data_rest <- data[-1]
  print("prior min")
  ## fit single exponential to lowest ATP conc
  data_min[, mod := lapply(spectra_prepped, fit_min_mod)]
  data_min[, predict_df := lapply(mod, predict_line, deadtime = 0.003)]
  data_min[, coef_kfast := sapply(mod, get_coef, x = "k")]
  data_min[, coef_a1 := sapply(mod, get_coef, x = "a")]
  data_min[, coef_c := sapply(mod, get_coef, x = "c")]
  data_min[, exp_type := "1-exp"]

  print("after min")
  ##get a1 value to lock for double exp
  a1 <- data_min$coef_a
  c_est <- data_min$coef_c

  print("before rest")
  ## fit the other atp concentration locking to min atp a1
  data_rest[, mod := lapply(spectra_prepped, fit_2exp, a1 = a1, c_est = c_est)]
  data_rest[, predict_df := lapply(mod, predict_line, deadtime = 0.00)]
  data_rest[, coef_kfast := sapply(mod, get_coef, x = "k1")]
  data_rest[, coef_kslow := sapply(mod, get_coef, x = "k2")]
  data_rest[, coef_a2 := sapply(mod, get_coef, x = "a2")]
  data_rest[, coef_a1 := a1]
  data_rest[, coef_c := sapply(mod, get_coef, x = "c")]
  data_rest[, exp_type := "2-exp"]

  print("after rest")
  ## unravel fit lines for plot
  data_predict <- data_rest[, predict_df[[1]], by = "atp"]

  ## unravel data for plotting
  data_unnest <- data_rest[, spectra_prepped[[1]], by = "atp"]
  data_unnest[, type := "real"]
  ## data_unnest$type <- "real"


  data_together <- rbind(data_unnest, data_predict)
  data_together$type <- factor(data_together$type, levels = c("real", "prediction"))

  data_to_return <- rbind(data_min, data_rest, fill = TRUE)

  return(list(
    data_nested_mod = data_to_return,
    data_to_plot = data_together)
    )
}

fit_kslow <- function(x){
nlsLM(coef_kslow ~ (k_alpha * atp_conc) / (k_not_used + atp_conc),
                 data = x,
                 start = list(k_alpha = 100,
                              k_not_used = 100))
}

predict_kslow <- function(x){
  data.table(x = seq(1, 6500, by = 1),
             y = predict(x, newdata = data.frame(atp_conc = seq(1, 6500, by = 1))))
}


## get_kslow_coef <- function(x){
## round(coef(x)[[1]], 1)
## }


## get_kfast_lm_coef <- function(x){
## round(coef(x)[["conc"]], 1)
## }


## get_k_2 <- function(x){
## round(coef(x)[["k_+2"]], 1)
## }


## get_K1 <- function(x){
## round(coef(x)[["K1"]], 1)
## }

folder <- file.path("stopped-flow", "atp-induced-dissociation")
files <- list.files(folder,
                    pattern = "data.csv",
                    recursive = TRUE,
                    full.names = FALSE)



files <- file.path(folder, files)


data <- rbindlist(lapply(files, read_sf))
path_split <- c("sf", "aid", "parent1", "parent2", "filename")
filename_split <- c("date", "id", "atp", "rep", "data")

data[, (path_split) := tstrsplit(path, "/", fixed = TRUE)]
data[, (filename_split) := tstrsplit(filename, "_", fixed = TRUE)]
## data[, rep := as.factor(sub(".csv", "", rep))]
data[, id2 := sub("uM-dani", " &micro;M", id)]

data[, atp := as.numeric(sub("uM-atp", "", atp))]

data[, `:=`(x_time_s = V1,
            y_signal = V2)]


data <- data[, .(y_signal = mean(y_signal)), by = c("date", "id", "id2", "atp", "x_time_s")]



#nest data
data <- data[, list(spectra = list(.SD)), by = c("date", "id", "id2", "atp")]

data <- data[,list(data_nest = list(.SD)), by = .(date, id, id2)]

data[, analyzed := lapply(data_nest, fit_atp_induced_dissociation)]

data_mod <- data[, analyzed[[1]][[1]], by = .(date, id, id2)]

## get spec files
spec_files <- list.files(folder,
                    pattern = "spec.csv",
                    recursive = TRUE,
                    full.names = FALSE)

spec_files <- file.path(folder, spec_files)
spec_data <- rbindlist(lapply(spec_files, fread))
spec_data <- spec_data[, .(spec_reading = mean(spec_reading)), by = .(date, atp, dilution_multiplier)]
spec_data[, date := as.character(date)]
spec_data[, atp := as.numeric(atp)]
spec_data[, atp_conc := (spec_reading*65*dilution_multiplier)/2]




#############################

#############################
## PLOT TRANSIENTS ##
#############################
data2plot <- data[, analyzed[[1]][[2]], by = .(date, id, id2)]
data2plot <- data2plot[x_time_s >= 0.00]
data2plot[, y_signal_norm := y_signal-min(y_signal), by = .(date, id, id2, atp)]
data2plot[, y_signal_norm := y_signal_norm/max(y_signal_norm), by = .(date, id, id2, atp)]
data2plot[, y_signal2 := y_signal/max(y_signal), by = .(date, id, id2, atp)]
data2plot[, x_time_ms := x_time_s*1000]
data2plot <- data2plot[spec_data, on = c("date", "atp"), nomatch = NULL]
data2plot[, atp_conc_label := paste0(round(atp_conc, 0), " &micro;M ATP")]
data2plot$atp_conc_label <- factor(data2plot$atp_conc_label, levels = c("23 &micro;M ATP",
                                                                        "102 &micro;M ATP",
                                                                         "1028 &micro;M ATP"))
data2plot$id2 <- factor(data2plot$id2, levels = c("0 &micro;M", "10 &micro;M"))

colorz <- c("#666666", "#e7298a")
basesize <- 11

(gg_transients <-
   ggplot(data2plot[date == "2024-03-25"  & atp %in% c(20, 100, 1000)])+
   geom_line(aes(x_time_ms,
                 y_signal_norm,
                 color = id2,
                 linetype = type,
                 alpha = type,
                 linewidth = type))+
   facet_wrap(~atp_conc_label)+
   scale_color_manual(values = colorz)+
   coord_cartesian(c(0, 122))+
   xlab("Time (ms)")+
   ylab("Normalized Fluorescence")+
   ## ggtitle("ATP Induced Dissociation")+
   scale_alpha_manual(values = c(0.5, 1))+
   scale_linetype_manual(values = c("solid", "dashed"))+
   scale_linewidth_manual(values = c(1, 0.5))+
   theme_cowplot(basesize)+
   theme(
     legend.position = "none",
     strip.background = element_rect(fill = "white"),
     strip.text = element_markdown(),
     plot.title = element_text(hjust = 0.5)
   )
)


#############################
## FIT KFAST ##
#############################

coef_data <- data_mod[, .(date, id, id2, atp, coef_kfast, coef_kslow, coef_a1, coef_a2 )]

# merge with spec files
coef_data <-  merge(coef_data, spec_data, by = c("date", "atp"))
coef_data$id2 <- factor(coef_data$id2, levels = c("0 &micro;M", "10 &micro;M"))

coef_data_nest <- coef_data[, .(data_nest = list(.SD)), by = .(date, id, id2)]

## estimate second order rate of atp induced dissociation using linear regression
## on low concentration points
coef_data_nest[, kfast_lm_fit := lapply(data_nest, fit_kfast_lm)]
coef_data_nest[, predict_kfast_lm := lapply(kfast_lm_fit, predict_kfast_lm)]
coef_data_nest[, coef_kfast_slope := sapply(kfast_lm_fit, get_coef, x = "atp_conc")]

## fit a hyperbolic (MM) to all the kfast rates (k1)
coef_data_nest[, kfast_nls_fit := lapply(data_nest, fit_kfast_nls)]
coef_data_nest[, predict_kfast_nls := lapply(kfast_nls_fit, predict_kfast_fun)]
coef_data_nest[, coef_k_plus_2 := sapply(kfast_nls_fit, get_coef, x = "k_+2")]
coef_data_nest[, coef_K1 := sapply(kfast_nls_fit, get_coef, x = "K1")]
coef_data_nest[, k_plus_2_over_K1 := coef_k_plus_2/coef_K1]


#################
## GLOBAL FIT
#####################
coef_data_nest2 <- coef_data[, .(data_nest = list(.SD)), by = .(id, id2)]

## estimate second order rate of atp induced dissociation using linear regression
## on low concentration points
coef_data_nest2[, kfast_lm_fit := lapply(data_nest, fit_kfast_lm)]
coef_data_nest2[, predict_kfast_lm := lapply(kfast_lm_fit, predict_kfast_lm)]
coef_data_nest2[, coef_kfast_slope := sapply(kfast_lm_fit, get_coef, x = "atp_conc")]

## fit a hyperbolic (MM) to all the kfast rates (k1)
coef_data_nest2[, kfast_nls_fit := lapply(data_nest, fit_kfast_nls)]
coef_data_nest2[, predict_kfast_nls := lapply(kfast_nls_fit, predict_kfast_fun)]
coef_data_nest2[, coef_k_plus_2 := sapply(kfast_nls_fit, get_coef, x = "k_+2")]
coef_data_nest2[, coef_K1 := sapply(kfast_nls_fit, get_coef, x = "K1")]
coef_data_nest2[, k_plus_2_over_K1 := coef_k_plus_2/coef_K1]


#############pseudo r2#
############################

modelr::rsquare(model = coef_data_nest2$kfast_nls_fit[[1]], data = coef_data_nest2$data_nest[[1]] )
modelr::rsquare(model = coef_data_nest2$kfast_nls_fit[[2]], data = coef_data_nest2$data_nest[[2]] )


coef_data_nest2[, kfast_lmod := lapply(data_nest, fit_lmod)]

summary(coef_data_nest2$kfast_lmod[[1]])
summary(coef_data_nest2$kfast_lmod[[2]])


##################

second_order_ttest <-
  t.test(x = coef_data_nest[id == "0uM-dani"]$k_plus_2_over_K1,
         y = coef_data_nest[id == "10uM-dani"]$k_plus_2_over_K1
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


gg_second_order_stats <-
ggplot()+
  geom_errorbar(data = coef_data_summary,
                aes(x = id2,
                    ymin = k_plus_2_over_K1 - sd_k_plus_2_over_K1,
                    ymax = k_plus_2_over_K1 + sd_k_plus_2_over_K1,
                    color = id2),
                width = 0.2)+
                geom_col(data = coef_data_summary,
                         aes(x = id2,
                             y = k_plus_2_over_K1,
                             fill = id2))+
                geom_jitter(data = coef_data_nest,
                            aes(x = id2, y =
                                           k_plus_2_over_K1),
                            width = 0.2, sized = 1)+
  draw_line(c(1, 2), c(6.1, 6.1))+
  annotate("text", x = 1.5, y = 6.5, label = paste0("p = ", round(second_order_ttest$p.value, 2)), size = 10/.pt)+
  xlab("[Dani]")+
  ylab("2<sup>nd</sup> Order Rate (&micro;M<sup>-1</sup>&middot;s<sup>-1</sup>)")+
  scale_y_continuous(expand = expansion(c(0, 0.1), c(0, 0.1)))+
  scale_fill_manual(values = colorz)+
  scale_color_manual(values = colorz)+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
  axis.text.x = element_markdown(),
  axis.title.y = element_markdown(size = 9),
  axis.title.x = element_markdown()
  )




kfast_lm_predict_df <- coef_data_nest[, predict_kfast_lm[[1]], by = c("date", "id", "id2")]
kfast_nls_predict_df <- coef_data_nest[, predict_kfast_nls[[1]], by = c("date", "id", "id2")]

kfast_lm_predict_df2 <- coef_data_nest2[, predict_kfast_lm[[1]], by = c("id", "id2")]
kfast_nls_predict_df2 <- coef_data_nest2[, predict_kfast_nls[[1]], by = c("id", "id2")]




(gg_second_order_1 <-
ggplot()+
   ## annotate("rect",
   ##          xmin = 0,
   ##          xmax = 300,
   ##          ymin = 0 ,
   ##          ymax = 450,
   ##          fill = "grey",
   ##          alpha = 0.6)+
  geom_point(data = coef_data,#[date == "2023-09-25"],
             aes(x = atp_conc,
                 y = coef_kfast,
                 color = id2),
             size = 2,
             shape = 16,
             alpha = 0.75)+
  geom_line(data = kfast_nls_predict_df2,
            aes(x = atp_conc,
                y = predict_y,
                color = id2),
            linewidth = 0.9)+
  ## facet_grid(~date)+
  ## ggtitle("k<sub>+2</sub>'")+
  ## ggtitle("ATP Induced Dissociation")+
  xlab(~ATP~(mu*M))+
  ylab("<i>k<sub>fast</i></sub> (s<sup>-1</sup>)")+
  scale_color_manual(values = colorz)+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
  plot.title = element_markdown(),
  strip.background = element_rect(fill = "transparent"),
  axis.title.y = element_markdown()
  ## axis.text.x = element_text(size = 14)
  )
)



## saveRDS(gg_second_order_1, "figures/gg_fkast.rds")










## (gg_second_order_1 <-
## ggplot()+
##    ## annotate("rect",
##    ##          xmin = 0,
##    ##          xmax = 300,
##    ##          ymin = 0 ,
##    ##          ymax = 450,
##    ##          fill = "grey",
##    ##          alpha = 0.6)+
##   geom_point(data = coef_data,#[date == "2023-09-25"],
##              aes(x = atp_conc,
##                  y = coef_kfast,
##                  color = id2),
##              size = 2,
##              shape = 16,
##              alpha = 1)+
##   geom_line(data = kfast_nls_predict_df[date == "2023-09-25"],
##             aes(x = atp_conc,
##                 y = predict_y,
##                 color = id2),
##             linewidth = 0.9)+
##   ## facet_grid(~date)+
##   ## ggtitle("k<sub>+2</sub>'")+
##   ## ggtitle("ATP Induced Dissociation")+
##   xlab(~ATP~(mu*M))+
##   ylab("<i>k<sub>fast</i></sub> (s<sup>-1</sup>)")+
##   scale_color_manual(values = colorz)+
##   theme_cowplot(basesize)+
##   theme(
##     legend.position = "none",
##   plot.title = element_markdown(),
##   strip.background = element_rect(fill = "transparent"),
##   axis.title.y = element_markdown()
##   ## axis.text.x = element_text(size = 14)
##   )
## )


(gg_second_order_2 <- gg_second_order_1+
   coord_cartesian(xlim = c(0, 150), ylim = c(0, 400))+
   ggtitle("")+
   xlab("")+
   ylab("")+
   theme_cowplot(10)+
  theme(
    legend.position = "none",
  plot.title = element_markdown(),
  axis.text.x = element_text(size = 10)
  )
  )




gg_kfast <- ggdraw(gg_second_order_1)+draw_plot(gg_second_order_stats, 0.65, 0.125, 0.25, 0.6)
gg_kfast

gg_kfast_with_stats <- plot_grid(gg_second_order_1,
                                 gg_second_order_stats,
                                 rel_widths = c(1, 0.5),
                                 align = "h",
                                 labels = c("B", "C"))

plot_grid(gg_transients, gg_kfast_with_stats,
          nrow = 2, rel_heights = c(1, 1), labels = c("A", ""))


pdf("figures/atp-induced-dissociation.pdf", width = 6, height = 5)
plot_grid(gg_transients, gg_kfast_with_stats,
          nrow = 2, rel_heights = c(1, 1), labels = c("A", ""))
dev.off()


png("figures/atp-induced-dissociation.png", width = 6, height = 5, units = "in", res = 300)
plot_grid(gg_transients, gg_kfast_with_stats,
          nrow = 2, rel_heights = c(1, 1), labels = c("A", ""))
dev.off()



pdf("figures/for-the-boss/atp-induced-dissociation.pdf", width = 6, height = 5)
plot_grid(gg_transients, gg_kfast_with_stats,
          nrow = 2, rel_heights = c(1, 1), labels = c("A", ""))
dev.off()



svg("figures/for-the-boss/atp-induced-dissociation.svg", width = 6, height = 5)
plot_grid(gg_transients, gg_kfast_with_stats,
          nrow = 2, rel_heights = c(1, 1), labels = c("A", ""))
dev.off()



#####################
## FIT KSLOW ##
####################


sf_data_no_10 <- coef_data[!is.na(as.numeric(coef_kslow))]

sf_data_no_10_nest <- sf_data_no_10[, list(data_nest = list(.SD)), by = .(date, id, id2)]
sf_data_no_10_nest[, kslow_mod := lapply(data_nest, fit_kslow)]
sf_data_no_10_nest[, predict_kslow_df := lapply(kslow_mod, predict_kslow)]
sf_data_no_10_nest[, coef_k_alpha := sapply(kslow_mod, get_coef, x = "k_alpha")]

kslow_predict_df <- sf_data_no_10_nest[, predict_kslow_df[[1]], by = .(date, id, id2)]



sf_data_no_10_nest2 <- sf_data_no_10[, list(data_nest = list(.SD)), by = .(id, id2)]
sf_data_no_10_nest2[, kslow_mod := lapply(data_nest, fit_kslow)]
sf_data_no_10_nest2[, predict_kslow_df := lapply(kslow_mod, predict_kslow)]
sf_data_no_10_nest2[, coef_k_alpha := sapply(kslow_mod, get_coef, x = "k_alpha")]
kslow_predict_df2 <- sf_data_no_10_nest2[, predict_kslow_df[[1]], by = .(id, id2)]

(gg_slow <-
ggplot()+
  geom_point(data = sf_data_no_10,
             aes(x = atp,
                 y = coef_kslow,
                 color = id2),
             alpha = 0.5,
             size = 2,
             ## color = "blue",
             shape = 16)+
  geom_line(data = kslow_predict_df2,
            aes(x = x,
                y = y,
                color = id2))+
  ## annotate("richtext",
  ##          x = Inf,
  ##          y = Inf,
  ##          label = paste0("<span style = 'color: black;'> pH 7.0 = ", sf_data_no_10_nest$coef_kslow[[1]], " s <sup>-1</sup></span><br>",
  ##                        "<span style = 'color: red;'> pH 6.5 = ", sf_data_no_10_nest$coef_kslow[[2]], " s <sup>-1</sup></span><br>"),
  ##          color = "transparent",
  ##          fill = "transparent",
  ##          hjust = 1,
  ##          vjust = 1
  ##          )+
  ## annotate("text",
  ##         x = Inf,
  ##         y = -Inf,
  ##         vjust = -1,
  ##         hjust =1,
  ##         parse = TRUE,
  ##         label = list(bquote(k[alpha]==.(k_alpha))))+
  ## coord_cartesian(ylim = c(0, 350))+
  ggtitle("k<sub>+&alpha;</sub>'")+
  xlab(~ATP~(mu*M))+
  ylab(~k[slow]~(s^-1))+
  ## scale_color_manual(values = colorz)+
  theme_cowplot(basesize)+
  theme(
    legend.position = "none",
  ## axis.text.x = element_text(size = 7),
  plot.title = element_markdown())
  )



kslow_data_sum <- sf_data_no_10_nest[, .(coef_k_alpha = mean(coef_k_alpha),
                                         sd_coef_k_alpha = sd(coef_k_alpha)),
                                     by = .(id, id2)]

t.test(sf_data_no_10_nest[id == "0uM-dani"]$coef_k_alpha,
       sf_data_no_10_nest[id == "10uM-dani"]$coef_k_alpha)

ratio_df <- na.omit(coef_data[, .(x = atp,
                                a1 = as.numeric(coef_a1),
                                a2 = as.numeric(coef_a2),
                                id = id,
                                date = date)])

ratio_df[, ratio_a := a1/a2]

ratio_plateau <- ratio_df[x >= 500, .(avg = mean(ratio_a)), by = .(date, id)]


ratio_plateau_sum <- ratio_plateau[, .(ratio = mean(avg),
                                  sd_ratio = sd(avg)),
                                  by = id]




t.test(ratio_plateau[id == "0uM-dani"]$avg,
       ratio_plateau[id == "10uM-dani"]$avg)



kslow_df <- ratio_plateau[sf_data_no_10_nest[, .(id, date, coef_k_alpha)], on = c("id", "date")]

kslow_df[, k_minus_alpha := coef_k_alpha/avg]

kslow_df_sum <- kslow_df[, .(k_minus_alpha = mean(k_minus_alpha),
                                  sd_k_minus_alpha = sd(k_minus_alpha)),
                                  by = id]


t.test(kslow_df[id == "0uM-dani"]$k_minus_alpha,
       kslow_df[id == "10uM-dani"]$k_minus_alpha    )



(gg_ratio <-
ggplot(ratio_df, aes(x, ratio_a, color = id))+
  geom_point(shape = 16, alpha = 0.5, size = 2)+
  geom_segment(data = ratio_plateau,
           aes(x = 500,
           xend = 6500,
           y = avg,
           yend = avg,
           color = id),
           linetype = "dashed")+
  ## annotate("richtext",
  ##          x = Inf,
  ##          y = Inf,
  ##          label = paste0("<span style = 'color: black;'> pH 7.0 = ", round(ratio_plateau$avg[[1]], 1), " s <sup>-1</sup></span><br>",
  ##                        "<span style = 'color: red;'> pH 6.5 = ", round(ratio_plateau$avg[[2]], 1), " s <sup>-1</sup></span><br>"),
  ##          color = "transparent",
  ##          fill = "transparent",
  ##          hjust = 1,
  ##          vjust = 1
  ##          )+
  ## scale_color_manual(values = colorz)+
## coord_cartesian(ylim = c(0, 20))+
  ## annotate("text",
  ##          x = Inf,
  ##          y = Inf,
  ##          vjust = 1,
  ##          hjust = 1,
  ##          label = paste0("Amplitude Ratio = ", round(ratio_plateau$avg, 1)))+
  ggtitle("K<sub>&alpha;</sub>'")+
  xlab(~ATP~(mu*M))+
  ylab(~A[fast]:A[slow])+
  theme_cowplot(basesize)+
  theme(
    plot.title = element_markdown(),
  axis.text.x = element_text(size = 7),
   legend.position = "none")
  )
