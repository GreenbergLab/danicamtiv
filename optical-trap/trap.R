library(data.table)
library(ggplot2)
library(cowplot)
library(ggtext)

options_paths <- list.files("~/lasertrapr/project_dani-single-molecule",
                            pattern = "options.csv",
                            full.names = TRUE,
                            recursive = TRUE)

options_data <- rbindlist(lapply(options_paths, fread), fill = TRUE)
options_data <- options_data[include == TRUE & review == TRUE & report == "success"]

me <- vector("list")
for(o in 1:nrow(options_data)){
  print(o)
  project <- options_data$project[[o]]
  conditions <- options_data$conditions[[o]]
  date <- options_data$date[[o]]
  obs <- options_data$obs[[o]]
  original_filename <- options_data$original_filename[[o]]
  path_to_read <- file.path("~", "lasertrapr", project, conditions, date, obs, "measured-events.csv")
  data_in <- fread(path_to_read)
  data_in[, original_filename := original_filename]
  me[[o]] <- data_in
}


me <- rbindlist(me, fill = TRUE)

me$conditions <- factor(me$conditions, levels = c("0uM-dani_10uM-atp", "1uM-dani_10uM-atp", "10uM-dani_10uM-atp", "100uM-dani_10uM-atp"))
me$conditions2 <- sub("_10uM-atp","", me$conditions)
me$conditions2 <- sub("10uM-atp_","", me$conditions2)

me$conditions2 <- factor(me$conditions2, levels = c("0uM-dani",
                                                     "1uM-dani",
                                                     "10uM-dani",
                                                     "100uM-dani"))

me[, original_filename := sub("  ", " ", original_filename)]
me[, c("filename", "rep") := tstrsplit(original_filename, " ", fixed = TRUE)]

me[, c("date2", "flowcell", "setup", "mogul", "exrta") := tstrsplit(filename, "_", fixed = TRUE)]

me[, ton := (attachment_duration_bead_1_s+attachment_duration_bead_2_s)/2]
me[, s1 := (substep_1_bead_1_nm+substep_1_bead_2_nm)/2]
me[, s2 := (substep_2_bead_1_nm+substep_2_bead_2_nm)/2]
me[, total_step := s1+s2]


me <- me[event_user_excluded == FALSE]




me_by_mol <- me[, .(med_ton = median(ton, na.rm=TRUE),
                    mean_ton = mean(ton, na.rm = TRUE),
                    mean_step = mean(total_step),
                    mean_s1 = mean(s1),
                    mean_s2 = mean(s2),
                    N = .N),
                    by = c("project",
                           "conditions",
                           "conditions2",
                           "date",
                           ## "obs",
                           "flowcell",
                           ## "setup",
                           "mogul"
                           )]
g1 <-
ggplot()+
  geom_jitter(data = me_by_mol,
              aes(conditions2,
                  mean_s1,
                  color = as.factor(date),
                  size = N),
                  ## shape = as.factor(date)),
              width = 0.2)


g2 <-
ggplot()+
  geom_jitter(data = me_by_mol,
              aes(conditions2,
                  mean_s2,
                  color = as.factor(date),
                  size = N),
                  ## shape = as.factor(date)),
              width = 0.2)

plot_grid(g1, g2)


me_by_mol <- me_by_mol[N >= 20]


me_by_con <- me_by_mol[, .(mean_med_ton = mean(med_ton),
                           sd_med_ton = sd(med_ton),
                    mean_mean_ton = mean(mean_ton),
                    med_med_ton = median(med_ton)),
                    by = c("project",
                           "conditions",
                           "conditions2")]

trap_10mM_pval <- t.test(me_by_mol[conditions2 == "0uM-dani"]$med_ton,
                        me_by_mol[conditions2 == "10uM-dani"]$med_ton,
                        )

kruskal.test(ton ~ conditions2, data = me[conditions2 %in% c("0uM-dani", "10uM-dani")])

ggplot()+
  geom_errorbar(data = me_by_con,
           aes(conditions2,
               ymin = mean_med_ton-sd_med_ton,
               ymax = mean_med_ton+sd_med_ton,
               color = conditions2),
           width = 0.2)+
  geom_col(data = me_by_con,
           aes(conditions2,
               mean_med_ton
               ## fill = conditions2
               ))+
  geom_jitter(data = me_by_mol,
              aes(conditions2,
                  med_ton,
                  color = as.factor(date),
                  size = N),
                  ## shape = as.factor(date)),
              width = 0.2)


ggplot()+
  geom_boxplot(data = me_by_mol,
               aes(conditions2,
                   med_ton,
                   fill = conditions2))


ggplot(me)+
  stat_ecdf(aes(ton, color = conditions))+
  coord_cartesian(c(0, 0.5))


############################################################
######## 1mM ATP
############################################################

options_paths <- list.files("~/lasertrapr/project_dani-1mM-atp",
                            pattern = "options.csv",
                            full.names = TRUE,
                            recursive = TRUE)

options_data <- rbindlist(lapply(options_paths, fread), fill = TRUE)
options_data <- options_data[include == TRUE & review == TRUE & report == "success"]

me <- vector("list")
for(o in 1:nrow(options_data)){
  print(o)
  project <- options_data$project[[o]]
  conditions <- options_data$conditions[[o]]
  date <- options_data$date[[o]]
  obs <- options_data$obs[[o]]
  original_filename <- options_data$original_filename[[o]]
  path_to_read <- file.path("~", "lasertrapr", project, conditions, date, obs, "measured-events.csv")
  data_in <- fread(path_to_read)
  data_in[, original_filename := original_filename]
  me[[o]] <- data_in
}


me <- rbindlist(me, fill = TRUE)

## me$conditions <- factor(me$conditions, levels = c("0uM-dani_10uM-atp", "1uM-dani_10uM-atp", "10uM-dani_10uM-atp", "100uM-dani_10uM-atp"))
me$conditions2 <- sub("control","0 &micro;M", me$conditions)
me$conditions2 <- sub("_1mM-atp","", me$conditions2)
me$conditions2 <- sub("10uM-dani","10 &micro;M", me$conditions2)

me$conditions2 <- factor(me$conditions2, levels = c("0 &micro;M", "10 &micro;M"))

me[, original_filename := sub("  ", " ", original_filename)]
me[, c("filename", "rep") := tstrsplit(original_filename, " ", fixed = TRUE)]

me[, c("date2", "flowcell", "setup", "mogul", "exrta") := tstrsplit(filename, "_", fixed = TRUE)]

me[, ton := (attachment_duration_bead_1_s+attachment_duration_bead_2_s)/2]




me_by_mol <- me[, .(med_ton = median(ton, na.rm=TRUE),
                    mean_ton = mean(ton, na.rm = TRUE),
                    N = .N),
                    by = c("project",
                           "conditions",
                           "conditions2",
                           "date",
                           ## "obs",
                           "flowcell",
                           "setup",
                           "mogul"
                           )]

me_by_mol <- me_by_mol[N >= 15]

me_by_con <- me_by_mol[, .(mean_med_ton = mean(med_ton),
                           sd_med_ton = sd(med_ton),
                    mean_mean_ton = mean(mean_ton),
                    med_med_ton = median(med_ton)),
                    by = c("project",
                           "conditions",
                           "conditions2"
                           )]



trap_1mM_pval <- t.test(me_by_mol[conditions == "control"]$med_ton,
                        me_by_mol[conditions == "10uM-dani_1mM-atp"]$med_ton,
                        )


ggplot()+
  geom_errorbar(data = me_by_con,
           aes(conditions2,
               ymin = mean_med_ton-sd_med_ton,
               ymax = mean_med_ton+sd_med_ton,
               color = conditions2),
           width = 0.2)+
  geom_col(data = me_by_con,
           aes(conditions2,
               mean_med_ton,
               fill = conditions2))+
  geom_jitter(data = me_by_mol,
              aes(conditions2,
                  med_ton,
                  shape = as.factor(date)),
              width = 0.2)+
  draw_line(c(1, 2), y = 0.04)+
  annotate("text", x = 1.5, y = 0.042, label = paste0("p = ", round(trap_1mM_pval$p.value, 2)))+
  scale_color_manual(values = c("#666666", "#e7298a"))+
  scale_fill_manual(values = c("#666666", "#e7298a"))+
  scale_y_continuous(expand = expansion(c(0, 0.01), c(0, 0.006)))+
  xlab("[Dani]")+
  ylab("Attachment Duration (s)")+
  ggtitle("1 mM ATP")+
  theme_cowplot(11)+
  theme(
   axis.text.x = element_markdown(),
   plot.title = element_text(hjust = 0.5),
   legend.position = "none")


ggplot()+
  geom_boxplot(data = me_by_mol,
               aes(conditions2,
                   med_ton,
                   fill = conditions2))


shapiro.test(me_by_mol[conditions == "10uM-dani_1mM-atp"]$med_ton)


########################################################3
########## 50 NM ATP
########################################################3

options_paths <- list.files("~/lasertrapr/project_dani-50nM-atp",
                            pattern = "options.csv",
                            full.names = TRUE,
                            recursive = TRUE)

options_data <- rbindlist(lapply(options_paths, fread), fill = TRUE)
options_data <- options_data[include == TRUE & review == TRUE & report == "success"]

me <- vector("list")
for(o in 1:nrow(options_data)){
  print(o)
  project <- options_data$project[[o]]
  conditions <- options_data$conditions[[o]]
  date <- options_data$date[[o]]
  obs <- options_data$obs[[o]]
  original_filename <- options_data$original_filename[[o]]
  path_to_read <- file.path("~", "lasertrapr", project, conditions, date, obs, "measured-events.csv")
  data_in <- fread(path_to_read)
  data_in[, original_filename := original_filename]
  me[[o]] <- data_in
}


me <- rbindlist(me, fill = TRUE)

## me$conditions <- factor(me$conditions, levels = c("0uM-dani_10uM-atp", "1uM-dani_10uM-atp", "10uM-dani_10uM-atp", "100uM-dani_10uM-atp"))
## me$conditions2 <- sub("control","0 &micro;M", me$conditions)
me$conditions2 <- sub("0uM-dani","0 &micro;M", me$conditions)
me$conditions2 <- factor(me$conditions2, levels = c("0 &micro;M", "100 &micro;M"))

me[, original_filename := sub("  ", " ", original_filename)]
me[, c("filename", "rep") := tstrsplit(original_filename, " ", fixed = TRUE)]

me[, c("date2", "flowcell", "setup", "mogul") := tstrsplit(filename, "_", fixed = TRUE)]

me[, ton := (attachment_duration_bead_1_s+attachment_duration_bead_2_s)/2]


me_by_mol <- me[, .(med_ton = median(ton, na.rm=TRUE),
                    mean_ton = mean(ton, na.rm = TRUE),
                    N = .N),
                    by = c("project",
                           "conditions",
                           "conditions2",
                           "date",
                           ## "obs",
                           "flowcell",
                           "setup",
                           "mogul"
                           )]

me_by_mol <- me_by_mol[N >= 15]

me_by_con <- me_by_mol[, .(mean_med_ton = mean(med_ton),
                           sd_med_ton = sd(med_ton),
                    mean_mean_ton = mean(mean_ton),
                    med_med_ton = median(med_ton),
                    N = .N ),
                    by = c("project",
                           "conditions",
                           "conditions2"
                           )]


trap_50nm_pval <- t.test(me_by_mol[conditions == "0uM-dani"]$med_ton,
                        me_by_mol[conditions == "100uM-dani"]$med_ton,
)

ggplot(me)+
  geom_point(aes(displacement_bead_1_nm, attachment_duration_bead_1_s),
             alpha = 0.2)+
  facet_wrap(~conditions2)


ggplot(me)+
  stat_ecdf(aes(ton, color = conditions2), pad = FALSE)



 kruskal.test(ton ~ conditions2, data = me)
ggplot(me)+
  geom_boxplot(aes(conditions2, ton))

ggplot()+
  geom_errorbar(data = me_by_con,
           aes(conditions2,
               ymin = mean_med_ton-sd_med_ton,
               ymax = mean_med_ton+sd_med_ton,
               color = conditions2),
           width = 0.2)+
  geom_col(data = me_by_con,
           aes(conditions2,
               mean_med_ton,
               fill = conditions2))+
  geom_jitter(data = me_by_mol,
              aes(conditions2,
                  med_ton,
                  shape = as.factor(date)),
              width = 0.2)+
  draw_line(c(1, 2), y = 0.04)+
  annotate("text", x = 1.5, y = 0.042, label = paste0("p = ", round(trap_50nm_pval$p.value, 2)))+
  scale_color_manual(values = c("#666666", "#e7298a"))+
  scale_fill_manual(values = c("#666666", "#e7298a"))+
  scale_y_continuous(expand = expansion(c(0, 0.01), c(0, 0.006)))+
  xlab("[Dani]")+
  ylab("Attachment Duration (s)")+
  ggtitle("50 nM ATP")+
  theme_cowplot(11)+
  theme(
   axis.text.x = element_markdown(),
   plot.title = element_text(hjust = 0.5),
   legend.position = "none")


ggplot()+
  geom_boxplot(data = me_by_mol,
               aes(conditions2,
                   med_ton,
                   fill = conditions2))


shapiro.test(me_by_mol[conditions == "10uM-dani_1mM-atp"]$med_ton)



########################################################3
########## SUBSTEPS
########################################################3

options_paths <- list.files("~/lasertrapr/project_dani-single-molecule",
                            pattern = "options.csv",
                            full.names = TRUE,
                            recursive = TRUE)

options_data <- rbindlist(lapply(options_paths, fread), fill = TRUE)
options_data <- options_data[include == TRUE & review == TRUE & report == "success"]

me <- vector("list")
for(o in 1:nrow(options_data)){
  print(o)
  project <- options_data$project[[o]]
  conditions <- options_data$conditions[[o]]
  date <- options_data$date[[o]]
  obs <- options_data$obs[[o]]
  original_filename <- options_data$original_filename[[o]]
  path_to_read <- file.path("~", "lasertrapr", project, conditions, date, obs, "substeps.csv")
  data_in <- fread(path_to_read)
  data_in[, original_filename := original_filename]
  me[[o]] <- data_in
}


me <- rbindlist(me, fill = TRUE)

me$conditions <- factor(me$conditions, levels = c("0uM-dani_10uM-atp", "1uM-dani_10uM-atp", "10uM-dani_10uM-atp", "100uM-dani_10uM-atp"))
me$conditions2 <- sub("_10uM-atp","", me$conditions)
me$conditions2 <- sub("10uM-atp_","", me$conditions2)

me$conditions2 <- factor(me$conditions2, levels = c("0uM-dani",
                                                     "1uM-dani",
                                                     "10uM-dani",
                                                     "100uM-dani"))

me[, original_filename := sub("  ", " ", original_filename)]
me[, c("filename", "rep") := tstrsplit(original_filename, " ", fixed = TRUE)]

me[, c("date2", "flowcell", "setup", "mogul", "exrta") := tstrsplit(filename, "_", fixed = TRUE)]

me <- me[, .(prior_unbound_position_nm = mean(prior_unbound_position_nm),
                                 bead_position_substep_1_nm = mean(bead_position_substep_1_nm),
                                 substep_1_nm = mean(substep_1_nm),
                                 bead_position_substep_2_nm = mean(bead_position_substep_2_nm),
                                 substep_2_nm = mean(substep_2_nm),
                                 after_unbound_position_nm = mean(after_unbound_position_nm)),
                             by = c("project", "conditions", "conditions2", "date", "obs", "event_id",
                                    "flowcell", "setup", "mogul")
                             ]

me[, total_step := substep_1_nm + substep_2_nm ]
me[, total_step2 := bead_position_substep_2_nm - after_unbound_position_nm ]
me[, substep_2_nm2 := total_step2 - substep_1_nm ]

## #average beads together
## me <- me[, .(substep_1_nm = mean(substep_1_nm),
##                     substep_2_nm = mean(substep_2_nm),
##              substep_2_nm2 = mean(substep_2_nm2),
##              total_step = mean(total_step),
##              total_step2 = mean(total_step)),
##                 by = c("project",
##                        "conditions",
##                        "conditions2",
##                        "date",
##                        ## "obs",
##                        "flowcell",
##                        "setup",
##                        "mogul",
##                        "event_id"
##                        )]

#average by molecule
me_by_mol <- me[, .(substep_1_nm = mean(substep_1_nm),
                    substep_2_nm = mean(substep_2_nm),
                    substep_2_nm2 = mean(substep_2_nm2, na.rm = TRUE),
                    total_step = mean(total_step),
                    total_step2 = mean(total_step2, na.rm = TRUE),
                    N = .N),
                by = c("project",
                       "conditions",
                       "conditions2",
                       "date",
                       ## "obs",
                       "flowcell",
                       ## "setup",
                       "mogul"
                       )]

me_by_mol <- me_by_mol[N >= 20]

me_by_con <- me_by_mol[, .(
  mean_s1 = mean(substep_1_nm),
  sd_s1 = sd(substep_1_nm),
  mean_s2 = mean(substep_2_nm),
  sd_s2 = sd(substep_2_nm),
  mean_s22 = mean(substep_2_nm2),
  sd_s22 = sd(substep_2_nm2),
  mean_total_step = mean(total_step),
  sd_total_step = sd(total_step),
  mean_total_step2 = mean(total_step2),
  sd_total_step2 = sd(total_step2)),
  by = c("project",
         "conditions",
         "conditions2")]



ggplot()+
  geom_errorbar(data = me_by_con,
           aes(conditions2,
               ymin = mean_s1-sd_s1,
               ymax = mean_s1+sd_s1,
               color = conditions2),
           width = 0.2)+
  geom_col(data = me_by_con,
           aes(conditions2,
               mean_s1,
               fill = conditions2))+
  geom_jitter(data = me_by_mol,
              aes(conditions2,
                  substep_1_nm,
                  shape = as.factor(date)),
              width = 0.2)





ggplot()+
  geom_errorbar(data = me_by_con,
           aes(conditions2,
               ymin = mean_s2-sd_s2,
               ymax = mean_s2+sd_s2,
               color = conditions2),
           width = 0.2)+
  geom_col(data = me_by_con,
           aes(conditions2,
               mean_s2,
               fill = conditions2))+
  geom_jitter(data = me_by_mol,
              aes(conditions2,
                  substep_2_nm,
                  shape = as.factor(date)),
              width = 0.2)





ggplot()+
  geom_boxplot(data = me_by_mol,
               aes(conditions2,
                   med_ton,
                   fill = conditions2))
