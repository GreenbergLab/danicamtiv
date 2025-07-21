library(data.table)
library(ggplot2)
library(ggtext)
library(cowplot)
library(RColorBrewer)
library(fitdistrplus)
library(magick)


theme_set(theme_cowplot(11))
colors <- brewer.pal(8, "Dark2")
colorz <- c(colors[8], colors[1], colors[4], colors[3], colors[2])


con_trace <- readRDS("~/lasertrapr/project_dani-single-molecule/summary/figures/0uM-dani_10uM-atp_2024-04-05_obs-07_50.157-52.6062.rds")+
  coord_cartesian(xlim = c(-0.02, 2.25))+
  draw_line(x = c(0, 0.2), y = -70)+
  annotate("text", x = 0.1, y = -75, label = "0.2 s", size = 3)+
  draw_line(x = -0.05, y = c(-30, -10))+
  annotate("text", x = -0.1, y = -20, label = "20 nm", angle = 90, size = 3)+
  ggtitle("0 &micro;M Dani")+
  theme_void(11)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold"))





con_trace_data <- con_trace$data

con_trace <-
  ggplot(con_trace_data)+
  geom_point(aes(new_time_index/20000, y = processed_bead_1, color = population),
             size = 0.25,
             alpha = 0.5,
             shape = 16)+
  coord_cartesian(xlim = c(-0.02, 2.25))+
  draw_line(x = c(0, 0.2), y = -70)+
  annotate("text", x = 0.1, y = -75, label = "0.2 s", size = 3)+
  draw_line(x = -0.05, y = c(-30, -10))+
  annotate("text", x = -0.1, y = -20, label = "20 nm", angle = 90, size = 3)+
  scale_color_manual(values = c("black", rep("#666666", 100)))+
  ggtitle("0 &micro;M Dani")+
  theme_void(10)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold"))




dani_trace <- readRDS("~/lasertrapr/project_dani-single-molecule/summary/figures/10uM-dani_10uM-atp_2023-11-06_obs-14_28.2474-30.0897.rds")+
  coord_cartesian(xlim = c(-0.04, 1.8))+
  draw_line(x = c(0, 0.2), y = -80)+
  annotate("text", x = 0.1, y = -85, label = "0.2 s", size = 3)+
  draw_line(x = -0.05, y = c(-37, -17))+
  annotate("text", x = -0.1, y = -27, label = "20 nm", angle = 90, size = 3)+
  ggtitle("10 &micro;M Dani")+
  theme_void(11)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold"))
dani_trace

dani_trace_data <- dani_trace$data

dani_trace <-
  ggplot(dani_trace_data)+
  geom_point(aes(new_time_index/20000, y = processed_bead_1, color = population),
             size = 0.25,
             alpha = 0.5,
             shape = 16)+
  coord_cartesian(xlim = c(-0.04, 1.8))+
  draw_line(x = c(0, 0.2), y = -80)+
  annotate("text", x = 0.1, y = -85, label = "0.2 s", size = 3)+
  draw_line(x = -0.05, y = c(-37, -17))+
  annotate("text", x = -0.1, y = -27, label = "20 nm", angle = 90, size = 3)+
  ggtitle("10 &micro;M Dani")+
  scale_color_manual(values = c("black", rep("#e7298a", 100)))+
  theme_void(10)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold"))

fig_top <- plot_grid(con_trace, dani_trace, labels = c("A", "B"))


gg_ea <- readRDS("~/lasertrapr/project_dani-single-molecule_10uM-atp/summary/figures/2024-07-10_15-10-46.249368_ensemble-average-plot.rds")

# remove grobs


gg_ea <-
gg_ea+
  facet_wrap(~conditions, labeller = as_labeller(c("0uM-dani_10uM-atp" = "0 &micro;M Dani",
                                                  "10uM-dani_10uM-atp" = "10 &micro;M Dani"
                                                  ## "10uM-OM_10uM-atp" = "10 &micro;M OM"
                                                  )))+
  xlab("")+
  ylab("Displacement (nm)")+
  draw_line(x = c(0, 0.1), y = -1.1)+
  annotate("text", x = 0.05, y = -1.6, label = "100 ms", size = 9/.pt)+
  scale_y_continuous(breaks = seq(0, 7, 1))+
  theme_cowplot(11)+
  theme(strip.text = element_markdown(face = "bold"),
        strip.background = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")


substeps_df <- fread("~/lasertrapr/project_dani-single-molecule_10uM-atp/summary/2024-06-03_all-substeps-ensemble-average.csv")

substeps_df$conditions <- factor(substeps_df$conditions, levels = c("0uM-dani_10uM-atp", "10uM-dani_10uM-atp", "10uM-OM_10uM-atp"))
substeps_df$conditions2 <- sub("_10uM-atp","", substeps_df$conditions)
## substeps_df$conditions2 <- sub("10uM-atp_","", substeps_df$conditions2)

substeps_df$conditions2 <- factor(substeps_df$conditions2, levels = c("0uM-dani",
                                                     "10uM-dani",
                                                       "10uM-OM"))

substeps_df$conditions3 <- sub("uM-dani", " &micro;M <br> Dani", substeps_df$conditions2)
substeps_df$conditions3 <- sub("uM-OM", " &micro;M <br> OM", substeps_df$conditions3)

substeps_df$conditions3 <- factor(substeps_df$conditions3, levels = c("0 &micro;M <br> Dani",
                                                                      "10 &micro;M <br> Dani",
                                                                      "10 &micro;M <br> OM"))


substeps_df <- substeps_df[, .(prior_unbound_position_nm = mean(prior_unbound_position_nm),
                                 bead_position_substep_1_nm = mean(bead_position_substep_1_nm),
                                 substep_1_nm = mean(substep_1_nm),
                                 bead_position_substep_2_nm = mean(bead_position_substep_2_nm),
                                 substep_2_nm = mean(substep_2_nm),
                                 after_unbound_position_nm = mean(after_unbound_position_nm),
                                 substep_2_nm_alt = mean(substep_2_nm_alt),
                                 total_step_nm_alt = mean(total_step_nm_alt),
                                 total_step_nm = mean(total_step_nm)),
                             by = c("project", "conditions", "conditions2", "conditions3",
                                    "date", "obs", "event_id")
                             ]

dani_10uM_substeps_summary <- substeps_df[, .(
  total = mean(total_step_nm),
  total_sd = sd(total_step_nm),
  total_se = sd(total_step_nm)/sqrt(.N),
  total2 = mean(total_step_nm_alt, na.rm = TRUE),
  total_sd2 = sd(total_step_nm_alt, na.rm = TRUE),
  total_se2 = sd(total_step_nm_alt, na.rm = TRUE)/sqrt(.N),
  s1 = mean(substep_1_nm),
  s1_sd = sd(substep_1_nm),
  s1_se = sd(substep_1_nm)/sqrt(.N),
  s2 = mean(substep_2_nm),
  s2_sd = sd(substep_2_nm, na.rm = TRUE),
  s2_se = sd(substep_2_nm, na.rm = TRUE)/sqrt(.N),
  s22 = mean(substep_2_nm_alt, na.rm = TRUE),
  s2_sd2 = sd(substep_2_nm_alt, na.rm = TRUE),
  s2_se2 = sd(substep_2_nm_alt, na.rm = TRUE)/sqrt(.N),
  n = paste0("n = ", .N)),
  by = .(project, conditions, conditions2)
  ]


t.test(substeps_df[conditions2=="0uM-dani"]$substep_2_nm,
       substeps_df[conditions2=="10uM-dani"]$substep_2_nm)


## dani_10uM_substeps_summary$conditions <- factor(dani_10uM_substeps_summary$conditions, levels = c("0uM-dani_10uM-atp",
##                                                                                                   "10uM-dani_10uM-atp",
##                                                                                                   "100uM-dani_10uM-atp"))

## dani_10uM_substeps_summary$conditions3 <- gsub("_10uM-atp", "", dani_10uM_substeps_summary$conditions)
## dani_10uM_substeps_summary$conditions3 <- factor(dani_10uM_substeps_summary$conditions2, levels = c("0uM-dani",
##                                                                                                     "10uM-dani",
##                                                                                                     "100uM-dani"))
## dani_10uM_substeps$conditions3 <- gsub("_10uM-atp", "", dani_10uM_substeps$conditions)
## dani_10uM_substeps$conditions3 <- factor(dani_10uM_substeps$conditions2, levels = c("0uM-dani",
##                                                                                                     "10uM-dani",
##                                                                                                     "100uM-dani"))
## ## dani_10uM_substeps_summary[, total_step := s1+s2]
## dani_10uM_substeps_summary_to_display <- dani_10uM_substeps_summary[, .(conditions,
##                                                              total_step = round(total, 1),
##                                                              substep1 = round(s1, 1),
##                                                              substep2 = round(s2, 1))]

## setorder(dani_10uM_substeps_summary, conditions)


rnorm_list <- vector("list")
  for(i in 1:nrow(dani_10uM_substeps_summary)){
    rnorm_list[[i]] <-
      data.table(total = rnorm(10000,
                                       mean = dani_10uM_substeps_summary$total[[i]],
                                       sd = dani_10uM_substeps_summary$total_sd[[i]]),
                 s1 = rnorm(10000,
                                     mean = dani_10uM_substeps_summary$s1[[i]],
                                     sd = dani_10uM_substeps_summary$s1_sd[[i]]
                                     ),
                 s2 = rnorm(10000,
                                     mean = dani_10uM_substeps_summary$s2[[i]],
                                     sd = dani_10uM_substeps_summary$s2_sd[[i]]
                                     ),
                 conditions = dani_10uM_substeps_summary$conditions[[i]],
                 conditions2 = dani_10uM_substeps_summary$conditions2[[i]]
                 )
  }

rnorm_df_10uM_atp_substeps <- do.call("rbind", rnorm_list)



## ggplot(substeps_df, aes(x = substep_1_nm, y = substep_2_nm))+
##   ## geom_tile(aes(x = substep_1_nm, y = substep_2_nm, fill = after_stat(density)))+
##   stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
##   scale_fill_viridis(discrete  = FALSE, option = "turbo", direction = 1)+
##   facet_wrap(~conditions2, ncol = 1)+
##   xlab("Substep 1 (nm)")+
##   ylab("Substep 2 (nm)")+
##   coord_cartesian(xlim = c(-20, 20), y = c(-20, 20))+
##   theme_cowplot()+
##   theme(
##    strip.background = element_rect(fill = "white")
##   )

(gg_total_ecdf <-
   ggplot(data = substeps_df[conditions2 != "10uM-OM"])+
   stat_ecdf(aes(total_step_nm, color = conditions2), linewidth = 1, alpha = 0.5)+
   stat_ecdf(data = rnorm_df_10uM_atp_substeps[conditions2 != "10uM-OM"],
               aes(total, color = conditions2),
             linetype = "dashed")+
   scale_color_manual(name = "", values = c(colorz[1], colorz[3], colorz[5]))+
   xlab("Total Displacement (nm)")+
   ylab("Probability")+
   ## ggtitle("Total Powerstroke")+
   theme(legend.position = "none",
         plot.title = element_text(hjust= 0.5))+
   coord_cartesian(xlim = c(-20, 30))
)



   ## ggplot(data = substeps_df)+
   ## geom_histogram(aes(total_step_nm, fill = conditions2), linewidth = 0.5, alpha = 1, color = "black", binwidth = 1)+
   ##   facet_wrap(~conditions2)
## basesize <- 11

me <- fread("~/lasertrapr/project_dani-single-molecule_10uM-atp/summary/2024-06-03_project_dani-single-molecule_10uM-atp_all-measured-events.csv")

me$conditions <- factor(me$conditions, levels = c("0uM-dani_10uM-atp",
                                                  "10uM-dani_10uM-atp",
                                                  "10uM-OM_10uM-atp"))

me[, ton := (attachment_duration_bead_1_s+attachment_duration_bead_2_s)/2]


## wilcox.test(x = me[conditions=="0uM-dani_10uM-atp"]$ton,
##             y = me[conditions=="10uM-dani_10uM-atp"]$ton)




fit_ton <- function(data, tmin){
                                        #from Schulte/Scott 2023 Biophysical Journal
                                        #get attachement times
  ton <- data$ton
  fit_ton_pdf_mle <- function(ton, k){
                                        # pass the pars through to the negative log likelihood function
                                        # optim will optimize these
                                        # the variables ton will be inherited from parent function (no need to pass directly)
    nll_fun <- function(k){
                                        # PDF function from SPASM
      -sum(log( (k*exp(-k*ton) )/ ( exp(-k*tmin) ) ))
    }

    fit <- optimize(nll_fun, lower = 0, upper = 100)

    return(fit)
  } #close fit_ton_pdf_mle

                                        # find k
  mod <- fit_ton_pdf_mle(ton = ton, k = 5)

                                        # extract fitted value from the model
  k <- mod$minimum[1]

                                        # function to generate fitted values to the exponential cumulative distribution
  fit_cdf <- function(x, k){ 1-exp(-k*x) }

                                        #calculate number of missing events
  n_missing <- fit_cdf(min(ton), k)*length(ton)

                                        #layout the range of x values for the real data bound by upper/lower bounds of data
  x_range <- seq(min(ton), max(ton), by = 1/1000)

                                        # generate the cdf values for the data
  cdf <- sapply(x_range, \(x) (sum(ton<=x)+n_missing) / (length(ton)+n_missing) )

  real_cdf <- data.frame(x = x_range,
                         y = cdf)



                                        # predict the fitted values from optimized points
  predict_x_range <- seq(0, max(ton), by = 1/1000)
  predict_y <- fit_cdf(k = k, x = predict_x_range)

  predict_df <- data.frame(x = predict_x_range,
                           y = predict_y)

#### BOOTSTRAP ####
  boostrap_ton_ci <- function(ton){

    boot_ton_fit <- function(ton){
      s <- sample(1:length(ton), replace = TRUE)
      new_ton <- ton[s]
      mod <- fit_ton_pdf_mle(ton = new_ton,
                             k = 5)
    } #close boot_ton_fit

    boot <- replicate(1000, boot_ton_fit(ton), simplify = FALSE)

    boot_df <- data.frame(k = sapply(boot, \(x) x$minimum[1]))

    ks <- sort(boot_df$k)

    k_lower <- ks[25]
    k_upper <- ks[975]

    return(list(boot_df = boot_df,
                k_95 = c(k_lower, k_upper)))

  } #close bootstrap_ton_ci

  ci <- boostrap_ton_ci(ton = ton)

  k_low_diff <- round(mod$minimum[[1]] - ci$k_95[[1]], 2)
  k_high_diff <- round(ci$k_95[[2]] - mod$minimum[[1]], 2)

  html_label <- paste0(round(mod$minimum[1], 1),
                       " (-",
                       round(k_low_diff, 1),
                       "/+",
                       round(k_high_diff, 1),
                       ")")

  parse_label <- paste0(round(mod$minimum[1], 1),
                        "~(-",
                        round(k_low_diff, 1),
                        "/+",
                        round(k_high_diff, 1),
                        ")")

  list(
    data_df = real_cdf,
    mod = mod,
    rate = mod$minimum[1],
    predict_df = predict_df,
    boot_df = ci$boot_df,
    boot_ci = list(k1_low = k_low_diff,
                   k1_up = k_high_diff),
    html_label = html_label,
    parse_label = parse_label
  )
}




                                        # fit the data
ton_boot <- me[, .(dat = list(.SD)), by = conditions]

ton_boot[, fit := lapply(dat, fit_ton, tmin = 0.01)]

ton_boot[, ton_rate := lapply(fit, `[[`, "rate")]
ton_boot[, boot_ci := lapply(fit, `[[`, "boot_ci")]
ton_boot[, predict := lapply(fit, `[[`, "predict_df")]
ton_boot[, parse_label  := sapply(fit, `[[`, "parse_label")]
ton_boot[, html_label  := sapply(fit, `[[`, "html_label")]
ton_boot[, boot_df := lapply(fit, `[[`, "boot_df")]
ton_boot[, cdf_data := lapply(fit, `[[`, "data_df")]

  ## ton_rate <-
  ##   ton_boot |>
  ##   dplyr::select(conditions, html_label)|>
  ##   tidyr::unnest(cols = c( html_label))


  ## ton_ci <-
  ##   ton_boot |>
  ##   dplyr::select(conditions, boot_ci)|>
  ##   tidyr::unnest(cols = c(boot_ci))
                                        # unravel data from the nest back to long data for ggplot
                                        # gets the real emperical CDF
  ton_real_df <-
    ton_boot[, cdf_data[[1]], by = conditions]

                                        # unravel data from the nest back to long data for ggplot
                                        # gets the fit prediction line
  ton_predict_df <-
    ton_boot[, predict[[1]], by = conditions]

                                        # make the plots
  ton_cdf <-
    ggplot()+
    geom_step(data = ton_real_df[conditions != "10uM-OM_10uM-atp"],
              aes(x = x*1000,
                  y = y,
                  color = conditions),
              alpha = 0.5,
              linewidth= 1)+
      ## stat_ecdf(aes(x), color = "red", pad = FALSE)+
    geom_line(data = ton_predict_df[conditions != "10uM-OM_10uM-atp"],
              aes(x = x*1000,
                  y = y,
                  color = conditions),
              linetype = "dashed",
              linewidth = 0.4)+
      coord_cartesian(xlim = c(0, 1000 ))+
    scale_color_manual(values = c(colorz[[1]], colorz[[3]], colorz[[5]]))+
    ## ggtitle("Attachment Durations")+
    ylab("Probability")+
    xlab("Attachment Duration (ms)")+
    theme_cowplot(11)+
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none")



gg_middle <- plot_grid(gg_total_ecdf, ton_cdf, labels = c("C", "D"))


gg_trap <- image_read_svg("figures/trap-cartoon.svg", height = 10000, width = 10000)

gg_trap <- ggplot()+draw_image(gg_trap)+theme_nothing()

traces_right <- plot_grid(con_trace, dani_trace, nrow = 2)

top <- plot_grid(gg_trap, traces_right, rel_widths = c(1, 1.2), labels = c("A", "B"))

bottom <- plot_grid(gg_middle, gg_ea, rel_widths = c(1.25, 1), nrow = 1, labels = c("", "E"))

plot_grid(top, bottom, nrow = 2, rel_heights = c(1, 1), labels = c("", "", ""))

png("figures/optical-trap-standard2.png", width = 6.5, height = 4.5, units = "in", res = 1000)
plot_grid(top, bottom, nrow = 2, rel_heights = c(1.25, 1), labels = c("", "", ""))
dev.off()



plot_grid(top, gg_middle, gg_ea, nrow = 3, rel_heights = c(1, 1, 1), labels = c("", "", "E"))











boot_df <- ton_boot[, boot_df[[1]], by = conditions]

## fwrite(boot_df[conditions == "0uM-dani_10uM-atp", .(k)], "optical-trap/0-dani-ton-boot.csv", col.names = FALSE)
## fwrite(boot_df[conditions == "10uM-dani_10uM-atp", .(k)], "optical-trap/10-dani-ton-boot.csv", col.names = FALSE)

## x <- vector()
## for(i in seq_len(1000)){
## x[[i]] <- rexp(1, rate = 70)+rexp(1, rate = 40)
## }


#stats step

t.test(x = substeps_df[conditions=="0uM-dani_10uM-atp"]$total_step_nm,
       y = substeps_df[conditions=="10uM-dani_10uM-atp"]$total_step_nm)

####################
#################
## 1mM ATP ####
#################
#################

me1 <- fread("~/lasertrapr/project_dani-1mM-atp/summary/2024-06-03_project_dani-1mM-atp_all-measured-events.csv")

me1$conditions <- factor(me1$conditions, levels = c("control",
                                                  ## "10uM-dani_10uM-atp",
                                                  "10uM-dani_1mM-atp"))

me1[, ton := (attachment_duration_bead_1_s+attachment_duration_bead_2_s)/2]
me1[, step := (displacement_bead_1_nm+displacement_bead_2_nm)/2]



t.test(me1[conditions=="control"]$step,
       me1[conditions=="10uM-dani_1mM-atp"]$step)

ton_boot1 <- me1[, .(dat = list(.SD)), by = conditions]

ton_boot1[, fit := lapply(dat, fit_ton, tmin = 0.01)]

ton_boot1[, ton_rate := lapply(fit, `[[`, "rate")]
ton_boot1[, boot_ci := lapply(fit, `[[`, "boot_ci")]
ton_boot1[, predict := lapply(fit, `[[`, "predict_df")]
ton_boot1[, parse_label  := sapply(fit, `[[`, "parse_label")]
ton_boot1[, html_label  := sapply(fit, `[[`, "html_label")]
ton_boot1[, boot_df := lapply(fit, `[[`, "boot_df")]
ton_boot1[, cdf_data := lapply(fit, `[[`, "data_df")]

  ## ton_rate <-
  ##   ton_boot |>
  ##   dplyr::select(conditions, html_label)|>
  ##   tidyr::unnest(cols = c( html_label))


  ## ton_ci <-
  ##   ton_boot |>
  ##   dplyr::select(conditions, boot_ci)|>
  ##   tidyr::unnest(cols = c(boot_ci))
                                        # unravel data from the nest back to long data for ggplot
                                        # gets the real emperical CDF
  ton_real_df1 <-
    ton_boot1[, cdf_data[[1]], by = conditions]

                                        # unravel data from the nest back to long data for ggplot
                                        # gets the fit prediction line
  ton_predict_df1 <-
    ton_boot1[, predict[[1]], by = conditions]

                                        # make the plots
    ggplot()+
    geom_step(data = ton_real_df1,
              aes(x = x,
                  y = y,
                  color = conditions),
              alpha = 0.5,
              linewidth= 1)+
      ## stat_ecdf(aes(x), color = "red", pad = FALSE)+
    geom_line(data = ton_predict_df1,
              aes(x = x,
                  y = y,
                  color = conditions),
              linetype = "dashed",
              linewidth = 0.4)+
      coord_cartesian(xlim = c(0, 1 ))+
    scale_color_manual(values = c(colorz[[1]], colorz[[3]], colorz[[5]]))+
    ## ggtitle("Attachment Durations")+
    ylab("Probability")+
    xlab("Attachment Duration (ms)")+
    theme_cowplot(11)+
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none")





boot_df1 <- ton_boot1[, boot_df[[1]], by = conditions]

fwrite(boot_df1[conditions == "control", .(k)], "optical-trap/0-dani-ton-boot_1mM-atp.csv", col.names = FALSE)
fwrite(boot_df1[conditions == "10uM-dani_1mM-atp", .(k)], "optical-trap/10-dani-ton-boot_1mM-atp.csv", col.names = FALSE)
