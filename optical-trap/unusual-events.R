library(data.table)
library(ggplot2)
library(cowplot)

## rename project and conditions
## trap_data_files <- list.files("~/lasertrapr/project_dani-single-molecule-filtered-in-unusual-events",
##                               pattern = "trap-data.csv",
##                               full.names = TRUE,
##                               recursive = TRUE)



## options_files <- list.files("~/lasertrapr/project_dani-single-molecule-filtered-in-unusual-events",
##                               pattern = "options.csv",
##                               full.names = TRUE,
##                               recursive = TRUE)


## me_files <- list.files("~/lasertrapr/project_dani-single-molecule-filtered-in-unusual-events",
##                               pattern = "measured-events.csv",
##                               full.names = TRUE,
##                               recursive = TRUE)


## substeps_files <- list.files("~/lasertrapr/project_dani-single-molecule-filtered-in-unusual-events",
##                               pattern = "substeps.csv",
##                               full.names = TRUE,
##                               recursive = TRUE)


## file_list <- c(trap_data_files, options_files, me_files, substeps_files)

## for(i in seq_along(file_list)){
##   print(i)
##   fname <- file_list[[i]]
##   f <- data.table::fread(fname)
##   ## f[,
##   ##   project = "project_dani-single-molecule_1uM-atp",
##   ##   ## conditions = "25_ug-ml",
##   ##   ## date  := "2022-02-22"
##   ## ]
##   ## f$date <- "2022-02-22"
##   f$project <- "project_dani-single-molecule-filtered-in-unusual-events"
##   ## f$conditions <- "10uM-dani_1mM-atp"
##   data.table::fwrite(f, fname)
## }

substeps_files <- list.files("~/lasertrapr/project_dani-single-molecule-filtered-in-unusual-events",
                             pattern = "substeps.csv",
                             full.names = TRUE,
                             recursive = TRUE)

filter_out_unusual_events <- function(file){
  ## file <- substeps_files[[1]]
  print(file)
  substep_data <- fread(file)

  substep_data <- substep_data[, .(prior_unbound_position_nm = mean(prior_unbound_position_nm),
                                   bead_position_substep_1_nm = mean(bead_position_substep_1_nm),
                                   substep_1_nm = mean(substep_1_nm),
                                   bead_position_substep_2_nm = mean(bead_position_substep_2_nm),
                                   substep_2_nm = mean(substep_2_nm),
                                   after_unbound_position_nm = mean(after_unbound_position_nm)),
                               by = c("project", "conditions", "date", "obs", "event_id")
                               ]

  substep_data[, total_step2 := bead_position_substep_2_nm - after_unbound_position_nm ]
  substep_data[, substep_2_nm2 := total_step2 - substep_1_nm ]
  substep_data <- substep_data[substep_2_nm2 <= 6]
  me_file <- sub("substeps.csv", "measured-events.csv", file)
  me_data <- fread(me_file)
  me_data$event_id <- 1:nrow(me_data)

  merge_data <- me_data[substep_data, on = .(project, conditions, date, obs, event_id)]
  ## merge_data <- dplyr::left_join(substep_data, me_data, by = c("project", "conditions", "date", "obs", "event_id"))

  fwrite(merge_data, file = me_file)

}

lapply(substeps_files, filter_out_unusual_events)






filter_in_unusual_events <- function(file){
  ## file <- substeps_files[[1]]
  print(file)
  substep_data <- fread(file)

  substep_data <- substep_data[, .(prior_unbound_position_nm = mean(prior_unbound_position_nm),
                                   bead_position_substep_1_nm = mean(bead_position_substep_1_nm),
                                   substep_1_nm = mean(substep_1_nm),
                                   bead_position_substep_2_nm = mean(bead_position_substep_2_nm),
                                   substep_2_nm = mean(substep_2_nm),
                                   after_unbound_position_nm = mean(after_unbound_position_nm)),
                               by = c("project", "conditions", "date", "obs", "event_id")
                               ]

  substep_data[, total_step2 := bead_position_substep_2_nm - after_unbound_position_nm ]
  substep_data[, substep_2_nm2 := total_step2 - substep_1_nm ]
  substep_data <- substep_data[substep_2_nm2 >= 6]
  me_file <- sub("substeps.csv", "measured-events.csv", file)
  me_data <- fread(me_file)
  me_data$event_id <- 1:nrow(me_data)

  merge_data <- me_data[substep_data, on = .(project, conditions, date, obs, event_id)]
  ## merge_data <- dplyr::left_join(substep_data, me_data, by = c("project", "conditions", "date", "obs", "event_id"))

  fwrite(merge_data, file = me_file)

}

lapply(substeps_files, filter_in_unusual_events)
################################################
################################################
################################################



ea <- fread("~/lasertrapr/project_dani-single-molecule/summary/ensemble-averages.csv")
ea$conditions <- factor(ea$conditions, levels = c("0uM-dani_10uM-atp", "10uM-dani_10uM-atp", "100uM-dani_10uM-atp"))


ea2 <- fread("~/lasertrapr/project_0.1uM-dani/summary/ensemble-averages.csv")

ea3 <- rbind(ea, ea2)

ea3$conditions2 <- sub("_10uM-atp","", ea3$conditions)
ea3$conditions2 <- sub("10uM-atp_","", ea3$conditions2)

ea3$conditions2 <- factor(ea3$conditions2, levels = c("0uM-dani",
                                                     "0.1uM-dani",
                                                     "10uM-dani",
                                                     "100uM-dani"))

ggplot(ea3[forward_backward_index>=-50])+
  geom_line(aes(x = forward_backward_index/20000,
                 y = avg,
                 color = conditions2),
             size = 0.5)+
  xlab("Time (s)")+
  ylab("Position (nm)")+
  ## geom_point(aes(x = forward_backward_index/20000,
  ##                y = avg,
  ##                color = conditions),
  ##            size = 0.5)+
  ## geom_hline(aes(yintercept=0))+
  ## geom_vline(aes(xintercept=0))+
  coord_cartesian(xlim=c(0,0.15))+
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2"), name="")+
  ## facet_wrap(~conditions)+
  theme_bw()


## ggsave("~/Downloads/ea-dani-overlayed.png", bg = "white")




ea <- fread("~/lasertrapr/project_dani-single-molecule-filtered-out-unusual-events/summary/ensemble-averages.csv")


ggplot(ea)+
  geom_line(aes(x = forward_backward_index/20000,
                 y = avg,
                 color = conditions),
             size = 0.5)+
  xlab("Time (s)")+
  ylab("Position (nm)")+
  ## geom_point(aes(x = forward_backward_index/20000,
  ##                y = avg,
  ##                color = conditions),
  ##            size = 0.5)+
  ## geom_hline(aes(yintercept=0))+
  ## geom_vline(aes(xintercept=0))+
  coord_cartesian(xlim=c(0,0.5))+
  ## facet_wrap(~conditions)+
  theme_bw()



ea <- fread("~/lasertrapr/project_dani-single-molecule-filtered-in-unusual-events/summary/ensemble-averages.csv")


ggplot(ea)+
  geom_line(aes(x = forward_backward_index/20000,
                 y = avg,
                 color = conditions),
             size = 0.5)+
  xlab("Time (s)")+
  ylab("Position (nm)")+
  ## geom_point(aes(x = forward_backward_index/20000,
  ##                y = avg,
  ##                color = conditions),
  ##            size = 0.5)+
  ## geom_hline(aes(yintercept=0))+
  ## geom_vline(aes(xintercept=0))+
  coord_cartesian(xlim=c(0,0.1))+
  ## facet_wrap(~conditions)+
  theme_bw()






#############################################

## me <- fread("~/lasertrapr/project_dani-single-molecule/summary/2024-01-10_project_dani-single-molecule_all-measured-events.csv")

## me[, displacement_nm := (displacement_bead_1_nm + displacement_bead_2_nm)/2 ]
## me[, time_on_s := (attachment_duration_bead_1_s + attachment_duration_bead_2_s)/2 ]
## me[, step1 := (substep_1_bead_1_nm + substep_1_bead_2_nm)/2 ]
## me[, step2 := (substep_2_bead_1_nm + substep_2_bead_2_nm)/2 ]

substeps_files <- list.files("~/lasertrapr/project_dani-single-molecule-filtered-in-unusual-events",
                             pattern = "substeps.csv",
                             full.names = TRUE,
                             recursive = TRUE)

substep_data <- rbindlist(lapply(substeps_files, fread))

substep_data <- substep_data[, .(prior_unbound_position_nm = mean(prior_unbound_position_nm),
                                 bead_position_substep_1_nm = mean(bead_position_substep_1_nm),
                                 substep_1_nm = mean(substep_1_nm),
                                 bead_position_substep_2_nm = mean(bead_position_substep_2_nm),
                                 substep_2_nm = mean(substep_2_nm),
                                 total_step = mean(total_step),
                                 after_unbound_position_nm = mean(after_unbound_position_nm)),
                             by = c("project", "conditions", "date", "obs", "event_id")
                             ]

substep_data[, total_step2 := bead_position_substep_2_nm - after_unbound_position_nm ]
substep_data[, substep_2_nm2 := total_step2 - substep_1_nm ]

substep_data$conditions <- factor(substep_data$conditions, levels = c("0uM-dani_10uM-atp",
                                                  "10uM-dani_10uM-atp",
                                                  "100uM-dani_10uM-atp"))

ggplot(substep_data)+
  geom_hex(aes(substep_1_nm, substep_2_nm2, fill = after_stat(density)))+
  facet_wrap(~conditions, ncol = 1)+
  coord_cartesian(xlim = c(-35, 35))+
  ## coord_cartesian(c(xlim = c(-10, 10)))+
  theme_linedraw()

ggplot(me)+
  geom_histogram(aes(step1, fill = conditions), bins = 30)+
  facet_wrap(~conditions)+
  theme_linedraw()






me_1mm <- fread("/home/brent/lasertrapr/project_dani-1mM-atp/summary/2024-04-02_project_dani-1mM-atp_all-measured-events.csv")

me_1mm$conditions <- factor(me_1mm$conditions, levels = c("control", "10uM-dani_1mM-atp"))
ggplot(me_1mm)+
  stat_ecdf(aes(attachment_duration_bead_1_s, color = conditions), pad = F, show.legend = F)+
  ylab("Cumulative Probability")+
  xlab("Time (s)")+
  scale_color_manual(values = c("#666666", "#e7298a"))+
  coord_cartesian(xlim = c(0, 0.5))+
  theme_cowplot()
