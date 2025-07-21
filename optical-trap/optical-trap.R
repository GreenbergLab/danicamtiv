library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)
library(RColorBrewer)
library(viridis)

## basesize <- 11

theme_set(theme_cowplot(11))

## ea <- fread("~/lasertrapr/project_dani-single-molecule/summary/2024-05-06_ensemble-averages.csv")
ea <- fread("~/lasertrapr/project_dani-single-molecule/summary/ensemble-averages.csv")
ea_om <- fread("~/lasertrapr/project_10uM-OM_10uM-atp/summary/ensemble-averages.csv")

ea <- rbind(ea, ea_om)

me <- fread("~/lasertrapr/project_dani-single-molecule/summary/2024-05-31_project_dani-single-molecule_all-measured-events.csv")
me_sum <- me[, .(n = .N), by = conditions]


me_om <- fread("~/lasertrapr/project_10uM-OM_10uM-atp/summary/2024-05-02_project_10uM-OM_10uM-atp_all-measured-events.csv")
me_om_sum <- me_om[, .(n = .N), by = conditions]

me_sum <- rbind(me_sum, me_om_sum)
## ea <- ea[conditions != "1uM-dani_10uM-atp"]

ea_sum <- ea[direction == "forward" & ensemble_index < 6000, .(y = mean(tail(avg))), by = .(conditions)]
ea_sum[, labels := sub("uM-dani_10uM-atp", " &micro;M Dani", conditions)]
ea_sum[, labels := sub("uM-OM_10uM-atp", " &micro;M OM", labels)]

ea_sum$labels <- factor(ea_sum$labels, levels = c("0 &micro;M Dani",
                                                  "1 &micro;M Dani",
                                                  "10 &micro;M Dani",
                                                  "100 &micro;M Dani",
                                                  "10 &micro;M OM"))

ea_sum <- ea_sum[me_sum, on = "conditions"]

ea_sum[, labels2 := paste0(labels, "<br> (n = ", n, ")")]


ea$conditions <- factor(ea$conditions, levels = c("0uM-dani_10uM-atp",
                                                  "1uM-dani_10uM-atp",
                                                  "10uM-dani_10uM-atp",
                                                  "100uM-dani_10uM-atp",
                                                  "10uM-OM_10uM-atp"))

colors <- brewer.pal(8, "Dark2")
colorz <- c(colors[8], colors[1], colors[4], colors[3], colors[2])

ea <- ea[conditions %in% c("0uM-dani_10uM-atp", "10uM-dani_10uM-atp", "10uM-OM_10uM-atp")]
ea_sum <- ea_sum[conditions %in% c("0uM-dani_10uM-atp", "10uM-dani_10uM-atp", "10uM-OM_10uM-atp")]

(
gg_ea <-
ggplot(ea[
  ## direction == "forward" &
          ## ensemble_index > -4
])+
          ## ensemble_index < 3000])+
  geom_point(aes(forward_backward_index/20000,
                 avg,
                 color = conditions),
             alpha = 0.8,
             ## shape = 16,
            linewidth = 0.5)+
  ## geom_richtext(data = ea_sum, aes(x = 0.37, y = y-0.15,
  ##                                  label = labels2,
  ##                                  color = conditions),
  ##               label.color = NA,
  ##               ## label.fill = NA,
  ##               label.padding = grid::unit(rep(0, 4), "pt"),
  ##               size = 3
  ##               ## label.padding = rep(0, 4)
  ##               )+
  scale_color_manual(values = c(colorz[1], colorz[3], colorz[5]))+
  scale_y_continuous(breaks = 0:7, expand = expansion(c(0, 0.01), c(0, 0.01)))+
  ## coord_cartesian(xlim = c(-0.002, 0.45), ylim = c(-0.8, 7))+
  draw_line(c(0.01, 0.11), c(-0.5, -0.5))+
annotate("text", x = 0.055, y = -0.7, label = "100 ms", size = 3)+
  xlab("")+
  ylab("Displacement (nm)")+
  ## theme_cowplot(12, font_family = "Arial")+
  theme(
  axis.line.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  legend.position = "none")
)
## ggsave("~/Downloads/2024-05-17_ensemble-average_control-10dani-10OM.pdf", device = cairo_pdf)


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
  theme_void(11)+
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
  theme_void(11)+
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold"))

fig_top <- plot_grid(con_trace, dani_trace, labels = c("A", "B"))


#######
######
#####
###SUBSTEPS
#####
#############

## options_paths <- list.files("~/lasertrapr/project_dani-single-molecule",
##                             pattern = "options.csv",
##                             full.names = TRUE,
##                             recursive = TRUE)

## options_data <- rbindlist(lapply(options_paths, fread), fill = TRUE)
## options_data <- options_data[include == TRUE & review == TRUE & report == "success"]

## me <- vector("list")
## for(o in 1:nrow(options_data)){
##   print(o)
##   project <- options_data$project[[o]]
##   conditions <- options_data$conditions[[o]]
##   date <- options_data$date[[o]]
##   obs <- options_data$obs[[o]]
##   original_filename <- options_data$original_filename[[o]]
##   path_to_read <- file.path("~", "lasertrapr", project, conditions, date, obs, "substeps.csv")
##   data_in <- fread(path_to_read)
##   data_in[, original_filename := original_filename]
##   me[[o]] <- data_in
## }


## me <- rbindlist(me, fill = TRUE)




substeps_df <- fread("~/lasertrapr/project_dani-single-molecule/summary/2024-05-31_all-substeps-ensemble-average.csv")

substeps_df_om <- fread("~/lasertrapr/project_10uM-OM_10uM-atp/summary/2024-05-02_all-substeps-ensemble-average.csv")

substeps_df <- rbind(substeps_df, substeps_df_om, fill = TRUE )

substeps_df$conditions <- factor(substeps_df$conditions, levels = c("0uM-dani_10uM-atp", "1uM-dani_10uM-atp", "10uM-dani_10uM-atp", "100uM-dani_10uM-atp", "10uM-OM_10uM-atp"))
substeps_df$conditions2 <- sub("_10uM-atp","", substeps_df$conditions)
## substeps_df$conditions2 <- sub("10uM-atp_","", substeps_df$conditions2)

substeps_df$conditions2 <- factor(substeps_df$conditions2, levels = c("0uM-dani",
                                                     "1uM-dani",
                                                     "10uM-dani",
                                                     "100uM-dani",
                                                       "10uM-OM"))

substeps_df$conditions3 <- sub("uM-dani", " &micro;M <br> Dani", substeps_df$conditions2)
substeps_df$conditions3 <- sub("uM-OM", " &micro;M <br> OM", substeps_df$conditions3)

substeps_df$conditions3 <- factor(substeps_df$conditions3, levels = c("0 &micro;M <br> Dani",
                                                                      "1 &micro;M <br> Dani",
                                                                      "10 &micro;M <br> Dani",
                                                                      "100 &micro;M <br> Dani",
                                                                      "10 &micro;M <br> OM"))

## me[, original_filename := sub("  ", " ", original_filename)]
## me[, c("filename", "rep") := tstrsplit(original_filename, " ", fixed = TRUE)]

## me[, c("date2", "flowcell", "setup", "mogul", "exrta") := tstrsplit(filename, "_", fixed = TRUE)]


substeps_df <- substeps_df[, .(prior_unbound_position_nm = mean(prior_unbound_position_nm),
                                 bead_position_substep_1_nm = mean(bead_position_substep_1_nm),
                                 substep_1_nm = mean(substep_1_nm),
                                 bead_position_substep_2_nm = mean(bead_position_substep_2_nm),
                                 substep_2_nm = mean(substep_2_nm),
                                 after_unbound_position_nm = mean(after_unbound_position_nm)),
                             by = c("project", "conditions", "conditions2", "conditions3",
                                    "date", "obs", "event_id")
                             ]

substeps_df[, total_step := substep_1_nm + substep_2_nm ]
substeps_df[, total_step2 := bead_position_substep_2_nm - after_unbound_position_nm ]
substeps_df[, substep_2_nm2 := total_step2 - substep_1_nm ]

(gg_bottom <- plot_grid(gg_ea, gg_step_bars, rel_widths = c(1.3,1 ), labels = c("C", "")))




##stats::
con_sub <- substeps_df[conditions == "0uM-dani_10uM-atp"]
dani10_sub <- substeps_df[conditions == "10uM-dani_10uM-atp"]
dani100_sub <- substeps_df[conditions == "100uM-dani_10uM-atp"]

## kruskal.test(total_step2 ~ conditions, data = substeps_df)
## ## dani_10uM_substeps[, total := substep_1_nm+substep_2_nm]
## pairwise.wilcox.test(dani_10uM_substeps$total_step2, dani_10uM_substeps$conditions, p.adjust.method = "bonferroni")

## kruskal.test(substep_1_nm ~ conditions, data = dani_10uM_substeps)
## ## dani_10uM_substeps[, total := substep_1_nm+substep_2_nm]
## pairwise.wilcox.test(dani_10uM_substeps$substep_1_nm, dani_10uM_substeps$conditions, p.adjust.method = "bonferroni")


## kruskal.test(substep_2_nm2 ~ conditions, data = dani_10uM_substeps)
## ## dani_10uM_substeps[, total := substep_1_nm+substep_2_nm]
## pairwise.wilcox.test(dani_10uM_substeps$substep_2_nm5, dani_10uM_substeps$conditions, p.adjust.method = "bonferroni")

dani_10uM_substeps_summary <- substeps_df[, .(
  total = mean(total_step),
  total_sd = sd(total_step),
  total_se = sd(total_step)/sqrt(.N),
  total2 = mean(total_step2, na.rm = TRUE),
  total_sd2 = sd(total_step2, na.rm = TRUE),
  total_se2 = sd(total_step2, na.rm = TRUE)/sqrt(.N),
  s1 = mean(substep_1_nm),
  s1_sd = sd(substep_1_nm),
  s1_se = sd(substep_1_nm)/sqrt(.N),
  s2 = mean(substep_2_nm),
  s2_sd = sd(substep_2_nm, na.rm = TRUE),
  s2_se = sd(substep_2_nm, na.rm = TRUE)/sqrt(.N),
  s22 = mean(substep_2_nm2, na.rm = TRUE),
  s2_sd2 = sd(substep_2_nm2, na.rm = TRUE),
  s2_se2 = sd(substep_2_nm2, na.rm = TRUE)/sqrt(.N),
  n = paste0("n = ", .N)),
  by = .(project, conditions, conditions2)
  ]

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
   ggplot(data = substeps_df[conditions %in% c("0uM-dani_10uM-atp",
                                               "10uM-dani_10uM-atp",
                                               "10uM-OM_10uM-atp")])+
   stat_ecdf(aes(total_step, color = conditions2), linewidth = 1, alpha = 0.5)+
   stat_ecdf(data = rnorm_df_10uM_atp_substeps[conditions %in% c("0uM-dani_10uM-atp",
                                                                 "10uM-dani_10uM-atp",
                                                                 "10uM-OM_10uM-atp")],
               aes(total, color = conditions2),
             linetype = "dashed")+
   scale_color_manual(name = "", values = c(colorz[1], colorz[3], colorz[5]))+
   xlab("Displacement (nm)")+
   ylab("Probability")+
   ## ggtitle("Total Displacement")+
   theme(legend.position = "none",
         plot.title = element_text(hjust= 0.5))+
   coord_cartesian(xlim = c(-20, 30))
)

ggplot(data = substeps_df[conditions %in% c("0uM-dani_10uM-atp",
                                            "10uM-dani_10uM-atp",
                                            "10uM-OM_10uM-atp")])+
  stat_ecdf(aes(substep_2_nm, color = conditions2), linewidth = 1)+
  scale_color_manual(name = "", values = c(colorz[1], colorz[3], colorz[5]))+
  xlab("Displacement (nm)")+
  ylab("Cumulative Probability")+
  coord_cartesian(xlim = c(-20, 30))+
  ggtitle("Substep 1")+
  theme(legend.position = "none",
        plot.title = element_text(hjust= 0.5))

ggplot(data = substeps_df[conditions %in% c("0uM-dani_10uM-atp",
                                            "10uM-dani_10uM-atp",
                                            "10uM-OM_10uM-atp")])+
  stat_ecdf(aes(substep_2_nm, color = conditions2), linewidth = 1)+
  scale_color_manual(name = "", values = c(colorz[1], colorz[3], colorz[5]))+
  xlab("Displacement (nm)")+
  ylab("Cumulative Probability")+
  coord_cartesian(xlim = c(-20, 30))+
  ggtitle("Substep 2")+
  theme(legend.position = "none",
        plot.title = element_text(hjust= 0.5))


ggplot(data = substeps_df[conditions %in% c("0uM-dani_10uM-atp",
                                            "10uM-dani_10uM-atp",
                                            "10uM-OM_10uM-atp")])+
  geom_histogram(aes(total_step, y = after_stat(density), fill = conditions),
                 color = "black",
                 binwidth = 4)+
  scale_color_manual(name = "", values = c(colorz[1], colorz[3], colorz[5]))+
  facet_wrap(~conditions)+
  xlab("Displacement (nm)")+
  ylab("Cumulative Probability")+
  coord_cartesian(xlim = c(-20, 30))

## ggsave("~/Downloads/total-step-sizes.pdf", device = cairo_pdf)
## png("figures/step-sizes.png", width = 6.5, height = 5.35, units = "in", res = 300)
## plot_grid(fig_top, gg_bottom, nrow = 2, rel_heights = c(1, 3))
## dev.off()

ggplot(substeps_df)+
  geom_histogram(aes(x = total_step, fill = conditions2), color = "black", binwidth = 3)+
  facet_wrap(~conditions2)

## pdf("figures/for-the-boss/step-sizes.pdf", width = 6.5, height = 5.35, useDingbats=FALSE)
## plot_grid(fig_top, gg_bottom, nrow = 2, rel_heights = c(1, 3))
## dev.off()

## svg("figures/for-the-boss/step-sizes.svg", width = 6.5, height = 5.35)
## plot_grid(fig_top, gg_bottom, nrow = 2, rel_heights = c(1, 3))
## dev.off()

## ggplot()+
##   geom_errorbar(data = me_by_con,
##            aes(conditions3,
##                ymin = mean_s2-sd_s2,
##                ymax = mean_s2+sd_s2,
##                color = conditions2),
##            width = 0.2)+
##   geom_col(data = me_by_con,
##            aes(conditions3,
##                mean_s2,
##                fill = conditions2))+
##   geom_jitter(data = me_by_mol,
##               aes(conditions2,
##                   substep_2_nm,
##                   shape = as.factor(date)),
##               width = 0.2)





## ggplot()+
##   geom_boxplot(data = me_by_mol,
##                aes(conditions2,
##                    med_ton,
##                    fill = conditions2))









## ggsave("~/Downloads/ea-draft.png", bg = "white")


#######

## me <- fread("~/lasertrapr/project_dani-single-molecule/summary/2024-05-02_project_dani-single-molecule_all-measured-events.csv")
## me <- me[conditions != "1uM-dani_10uM-atp"]

## me_om <- fread("~/lasertrapr/project_10uM-OM_10uM-atp/summary/2024-04-11_project_10uM-OM_10uM-atp_all-measured-events.csv")

me <- rbind(me, me_om)

me$conditions <- factor(me$conditions, levels = c("0uM-dani_10uM-atp",
                                                  "1uM-dani_10uM-atp",
                                                  "10uM-dani_10uM-atp",
                                                  "100uM-dani_10uM-atp",
                                                  "10uM-OM_10uM-atp"))

me[, ton := (attachment_duration_bead_1_s+attachment_duration_bead_2_s)/2]


wilcox.test(x = me[conditions=="0uM-dani_10uM-atp"]$ton,
            y = me[conditions=="10uM-dani_10uM-atp"]$ton)


library(fitdistrplus)

data <- me[conditions=="0uM-dani_10uM-atp"]$ton

fit_gamma <- function(data){
  gamma_fit_real <- fitdist(data, distr = "gamma", method = "mle")
  rate <- coef(gamma_fit_real)[[2]]
  sample_and_fit_gamma <- function(data){
    new_data <- sample(data, replace = T)
    gamma_fit <- fitdist(new_data, distr = "gamma", method = "mle")
    dt <- data.table(shape = coef(gamma_fit)[1],
                     rate = coef(gamma_fit)[2])
    return(dt)
  }

  all_reps <- replicate(1000, sample_and_fit_gamma(data = data), simplify = FALSE )

  dt <- rbindlist(all_reps)
  setorder(dt, rate)

  rate_lower <- dt$rate[25]
  rate_upper <- dt$rate[975]

  return(list(rate = rate,
              rate_lower = rate_lower,
              rate_upper = rate_upper,
              boot = dt$rate))

}

ton_fit_0uM_dani <- fit_gamma(data = me[conditions=="0uM-dani_10uM-atp"]$ton)
ton_fit_10uM_dani <- fit_gamma(data = me[conditions=="10uM-dani_10uM-atp"]$ton)


fwrite(data.table(x = ton_fit_0uM_dani$boot), "~/Downloads/try-0-dani.csv", col.names = FALSE)
fwrite(data.table(x = ton_fit_10uM_dani$boot), "~/Downloads/try-10-dani.csv", col.names = FALSE)

summary(fit.gamma)

plot(fit.gamma)









me_sum <- me[, .(med_ton = median(ton),
                 med_ton_ms = round(median(ton)*1000),
                 quarter = quantile(ton, probs = 0.25)),
             by = conditions]

fwrite(me[conditions == "0uM-dani_10uM-atp", .(ton)],
       "optical-trap/0uM-10uM-atp.csv",
       col.names = F)


fwrite(me[conditions == "10uM-dani_10uM-atp", .(ton)],
       "optical-trap/10uM-10uM-atp.csv",
       col.names = F)


fwrite(me[conditions == "10uM-OM_10uM-atp", .(ton)],
       "optical-trap/10uM-OM_10uM-atp.csv",
       col.names = F)


me_sum <- me[, .(avg = mean(ton*1000),
                 sd = sd(ton*1000),
                 n = .N),
             by = conditions]

me_sum[, labels := sub("uM-dani_10uM-atp", " &micro;M Dani", conditions)]
me_sum[, labels := sub("uM-OM_10uM-atp", " &micro;M OM", labels)]

me_sum$labels <- factor(me_sum$labels, levels = c("0 &micro;M Dani",
                                                  "1 &micro;M Dani",
                                                  "10 &micro;M Dani",
                                                  "100 &micro;M Dani",
                                                  "10 &micro;M OM"))


(gg_ton <-
ggplot(me[conditions %in% c("0uM-dani_10uM-atp", "10uM-dani_10uM-atp", "10uM-OM_10uM-atp")])+
  stat_ecdf(aes(ton, color = conditions), pad = F, linewidth = 1)+
  coord_cartesian(xlim = c(0,0.8))+
  scale_color_manual(values = c(colorz[1], colorz[3], colorz[5]))+
  xlab("Attachment Duration (s)")+
  ylab("Probability")+
  ## ggtitle("Attachment Durations")+
  theme(legend.position = "none",
        plot.title = element_text(hjust= 0.5))
  )



gg_ton
## ggsave("~/Downloads/10uM-atp-attachment-durations.pdf", device = cairo_pdf)

kruskal.test(ton ~ conditions, data = me)

dunn.test::dunn.test(me$ton, me$conditions, method = "bonferroni")

stats_df <- data.frame(x = c(2, 3, 4),
                       y = c(0.08, 0.082, 0.098),
                       label = c("*", "*", "*#"))


## ggplot(me)+
##   geom_boxplot(aes(x = conditions, y = ton, fill = conditions))+
##   scale_y_log10()+
##   scale_fill_manual(values = colorz)
##   ## scale_y_log10()


## ggplot(me_sum)+
##   geom_errorbar(aes(labels, ymin = med_ton-quarter, ymax = med_ton+quarter, color = conditions), width = 0.2)+
##   geom_col(aes(x = labels, y = med_ton, fill = conditions))+
##   ## geom_text(data = stats_df, aes(x = x, y = y, label = label))+
##   scale_color_manual(values = colorz)+
##   scale_fill_manual(values = colorz)+
##   ylab("Time (s)")+
##   xlab("")+
##   scale_y_continuous(expand = expansion(c(0, 0.01), c(0, 0.01)))+
##  theme_cowplot(basesize)+
##   theme(
##     axis.text.x = element_markdown(),
##     legend.position = "none"
##   )



bottom_right <- plot_grid(gg_total_ecdf, gg_ton, nrow = 2, labels = c("D", "E"))

bottom <- plot_grid(gg_ea, bottom_right, labels = c("C", ""), rel_widths = c(1, 0.75))

png("figures/optical-trap.png", width = 6.5, height = 5.5, units = "in", res = 300)
plot_grid(fig_top, bottom, nrow = 2, rel_heights = c(0.3, 1))
dev.off()


############
# 1mM  ATP
########

me_1mM <- fread("~/lasertrapr/project_dani-1mM-atp/summary/2024-05-28_project_dani-1mM-atp_all-measured-events.csv")
## me <- me[conditions != "1uM-dani_10uM-atp"]

me_om_1mM <- fread("~/lasertrapr/project_10uM-OM_1mM-atp/summary/2024-04-11_project_10uM-OM_1mM-atp_all-measured-events.csv")

me_1mM <- rbind(me_1mM, me_om_1mM, fill = TRUE)

me_1mM[, ton := (attachment_duration_bead_1_s+attachment_duration_bead_2_s)/2]

me_1mM_sum <- me_1mM[, .(avg = mean(ton*1000),
                 sd = sd(ton*1000),
                 n = .N),
             by = conditions]


me_1mM[, labels := sub("uM-dani_1mM-atp", " &micro;M Dani", conditions)]
me_1mM[, labels := sub("uM-OM_1mM-atp", " &micro;M OM", labels)]
me_1mM[, labels := sub("control", "0 &micro;M Drug", labels)]

me_sum_1mM <- me_1mM[, .(med_ton = median(ton),
                 quarter = quantile(ton, probs = 0.25)),
             by = conditions]

me_sum_1mM[, labels := sub("uM-dani_1mM-atp", " &micro;M Dani", conditions)]
me_sum_1mM[, labels := sub("uM-OM_1mM-atp", " &micro;M OM", labels)]
me_sum_1mM[, labels := sub("control", "0 &micro;M Dani", labels)]

me_sum_1mM$labels <- factor(me_sum_1mM$labels, levels = c("0 &micro;M Dani",
                                                  "10 &micro;M Dani",
                                                  ## "100 &micro;M Dani",
                                                  "10 &micro;M OM"))

ggplot(me_1mM)+
  stat_ecdf(aes(ton, color = labels), pad = F, linewidth = 1, alpha = 0.8)+
  coord_cartesian(xlim = c(0,0.4))+
  scale_color_manual(values = c(colorz[1], colorz[3], colorz[5]))+
  ylab("Cumulative Probability")+
  xlab("Time (s)")+
  ggtitle("1 mM ATP")+
  ## theme_cowplot(12)+
  theme(
    legend.text = element_markdown(),
    plot.title = element_text(hjust = 0.5)
  )

## extrafont::loadfonts()
ggsave("~/Downloads/1mM-atp-trap-attachment-durations.pdf", device = cairo_pdf)

kruskal.test(ton ~ conditions, data = me_1mM)

dunn.test::dunn.test(me_1mM$ton, me_1mM$conditions, method = "bonferroni")

stats_df <- data.frame(x = c(3),
                       y = c(0.06),
                       label = c("*#"))


ggplot(me_sum_1mM)+
  geom_errorbar(aes(labels, ymin = med_ton-quarter, ymax = med_ton+quarter, color = labels), width = 0.2)+
  geom_col(aes(x = labels, y = med_ton, fill = labels))+
  ## geom_text(data = stats_df, aes(x = x, y = y, label = label))+
  scale_color_manual(values = c(colorz[1], colorz[3], colorz[5]))+
  scale_fill_manual(values = c(colorz[1], colorz[3], colorz[5]))+
  ## scale_fill_manual(values = colorz)+
  ylab("Time (s)")+
  xlab("")+
  scale_y_continuous(expand = expansion(c(0, 0.01), c(0, 0.01)))+
 theme_cowplot(basesize)+
  theme(
    axis.text.x = element_markdown(),
    legend.position = "none"
  )
