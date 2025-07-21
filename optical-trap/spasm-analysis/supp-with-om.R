library(ggplot2)
library(ggtext)
library(cowplot)
library(tidyr)
library(minpack.lm)
library(data.table)
library(magick)
devtools::load_all("~/git/spasmtools")

colorz <- c("#666666", "#e7298a", "#d95f02")


spasm_files <- list("Control" = "optical-trap/spasm-analysis/control/spasm-included.xlsx",
                    "10 &micro;M Dani" = "optical-trap/spasm-analysis/10uM-dani/spasm-included-manual.xlsx",
                    "10 &micro;M OM" = "optical-trap/spasm-analysis/10uM-om/spasm-included.xlsx"
                    )

spasm_data <- spasm_read(spasm_files)

## spasm_data <- spasm_data |> separate(id, c("isoform", "construct"), sep = "-", remove = FALSE)

con <- spasm_ensemble_average(spasm_file = "optical-trap/spasm-analysis/control/combinedEnsembleAxesData.xlsx",
                              forward_time_filter = 0.3,
                              reverse_time_filter = 0.5,
                              step_estimation_method = "fit",
                              color = colorz[1],
                              textsize = 8,
                              ylim = c(-2, 9),
                              x_shift = 0.9,
                              title = "0 &micro;M"
                              )


dani <- spasm_ensemble_average(spasm_file = "optical-trap/spasm-analysis/10uM-dani/manual-combinedEnsembleAxesData.xlsx",
                              forward_time_filter = 0.3,
                              reverse_time_filter = 0.5,
                              step_estimation_method = "fit",
                              color = colorz[2],
                              textsize = 8,
                              ylim = c(-2, 9),
                              x_shift = 0.9,
                              title = "10 &micro;M"
                              )


om <- spasm_ensemble_average(spasm_file = "optical-trap/spasm-analysis/10uM-om/combinedEnsembleAxesData.xlsx",
                              forward_time_filter = 0.3,
                              reverse_time_filter = 0.5,
                              step_estimation_method = "fit",
                              color = colorz[3],
                              textsize = 8,
                              ylim = c(-2, 9),
                              x_shift = 0.9,
                              title = "10 &micro;M"
                              )



con_ea_f <- con$ea_f
con_ea_f$id <- "Control"

dani_ea_f <- dani$ea_f
dani_ea_f$id <- "10 &micro;M Dani"

om_ea_f <- om$ea_f
om_ea_f$id <- "10 &micro;M OM"

ea_dat <- rbind(con_ea_f, dani_ea_f, om_ea_f)

ea_dat$id <- factor(ea_dat$id, levels = c("Control", "10 &micro;M Dani", "10 &micro;M OM"))

con_f_pred <- con$f_pred
dani_f_pred <- dani$f_pred
om_f_pred <- om$f_pred

(
ggea <-
ggplot(data = dplyr::filter(ea_dat, seconds > -0.02))+
    geom_line(
              aes(x = seconds,
                  y = nanometers,
                  color = id),
              ## alpha = 0.8,
              linewidth = 0.6,
              show.legend = FALSE)+
    annotate("text",
             x = 0.19,
             y = min(ea_dat$nanometers),
             label = "0.25 s",
             vjust = 1.25,
             hjust = 0.5,
             color = "black",
             size = 12/.pt)+
    geom_line(data = con_f_pred,
              aes(x = seconds,
                  y = y_fit),
              linewidth = 0.4)+
    geom_line(data = dani_f_pred,
              aes(x = seconds,
                  y = y_fit),
              linewidth = 0.4)+
    geom_line(data = om_f_pred,
              aes(x = seconds,
                  y = y_fit),
              linewidth = 0.4)+
  draw_line(c(0.1, 0.26), c(-0.5))+
  coord_cartesian(ylim = c(-1, 5.75))+
  scale_y_continuous(breaks = 0:5)+
  scale_color_manual(values = colorz)+
  ylab("Displacement (nm)")+
  theme_cowplot(11)+
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
)




ton <- spasm_fit_ton(spasm_data = spasm_data, tmin = 0.01, colorz = colorz)

ggton <-
ton$plot+
coord_cartesian(xlim = c(0, 1))+
ggtitle("")+
ylab("Probability")+
xlab("Attachment Duration (s)")+
theme_cowplot(11)+
  theme(legend.position = "none")



## con_to_write <-dplyr::filter(ton$boot_data, id == "0 &micro;M")$boot_df[[1]]

## fwrite(con_to_write, file = "optical-trap/spasm-analysis/control/ton-boot-control.csv", col.names = F)


## dani_to_write <-dplyr::filter(ton$boot_data, id == "10 &micro;M")$boot_df[[1]]
## fwrite(dani_to_write, file = "optical-trap/spasm-analysis/10uM-dani/ton-boot-dani.csv", col.names = F)


spasm_summary <-
  spasm_data |>
  dplyr::group_by(id)|>
  dplyr::summarize("Avg Total Step" = mean(Displacements),
                   "SD Total Step" = sd(Displacements, na.rm = T),
                   "Avg Substep 1" = mean(`Substep 1`),
                   "SD Substep 1" = sd(`Substep 1`, na.rm = T),
                   "Avg Substep 2" = mean(`Substep 2`),
                   "SD Substep 2" = sd(`Substep 2`, na.rm = T),
                   n = n()
                   )

rnorm_list <- vector("list")
  for(i in 1:nrow(spasm_summary)){
    rnorm_list[[i]] <-
      data.table(total = rnorm(10000,
                                       mean = spasm_summary$`Avg Total Step`[[i]],
                                       sd = spasm_summary$`SD Total Step`[[i]]),
                 id = spasm_summary$id[[i]]
                 )
  }
rnorm_df <- do.call("rbind", rnorm_list)

ggstep1 <-
ggplot(spasm_data)+
  stat_ecdf(aes(x = Displacements, color = id), linewidth = 1, alpha = 0.7)+
  stat_ecdf(data = rnorm_df, aes(x = total, color = id), linewidth = 0.5, linetype = "dashed")+
  scale_color_manual(values = colorz)+
  coord_cartesian(xlim = c(-15, 20))+
  ylab("Probability")+
  xlab("Total Displacement (nm)")+
  ggtitle("")+
  theme_cowplot(11)+
  theme(
  legend.position = "none"
  )

ggstep2 <- ggstep1+coord_cartesian(xlim = c(-10, 15))+theme_cowplot(11)+theme(legend.position = "none")

ggstep3 <- ggdraw(ggstep2)+draw_plot(ggstep1, 0.6, 0.15, 0.4, 0.4)

## d0_facet <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = Displacements, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   facet_wrap(~isoform)+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Total Step")+
##   theme_cowplot()

## d1 <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = `Substep 1`, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Substep 1")+
##   theme_cowplot()+
##   theme(
##   legend.position = "none"
##   )



## d1_facet <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = `Substep 1`, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   facet_wrap(~isoform)+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Substep 1")+
##   theme_cowplot()

## d2 <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = `Substep 2`, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Substep 2")+
##   theme_cowplot()

## d2_facet <-
## ggplot(spasm_data)+
##   stat_ecdf(aes(x = `Substep 2`, color = id), linewidth = 0.8)+
##   scale_color_manual(values = colorz)+
##   coord_cartesian(xlim = c(-20, 25))+
##   facet_wrap(~isoform)+
##   ylab("Probability")+
##   xlab("Displacement (nm)")+
##   ggtitle("Substep 2")+
##   theme_cowplot()


png("supplemental-figures/optical-trap-supplment-om.png", width = 6.5, height = 2.5, units = "in", res = 1000)
plot_grid(ggea, ggstep1, ggton, nrow = 1, labels = c("A", "B", "C"))
dev.off()






ton_df <- dplyr::select(ton$boot_data, id, html_label)

spasm_summary <-
  spasm_data |>
  dplyr::group_by(id)|>
  dplyr::summarize("Avg Total Step" = mean(Displacements),
                   "SD Total Step" = sd(Displacements, na.rm = T),
                   "Avg Substep 1" = mean(`Substep 1`),
                   "SD Substep 1" = sd(`Substep 1`, na.rm = T),
                   "Avg Substep 2" = mean(`Substep 2`),
                   "SD Substep 2" = sd(`Substep 2`, na.rm = T),
                   n = n()
                   )


dat <- left_join(spasm_summary, ton_df)


write.csv(dat, "supplementary-summary-trap-table.csv")






con$mod_f



###########################
############
# stats
##########################
###########################


atest <- aov(Displacements ~ id, data = spasm_data)

tuks <- TukeyHSD(atest)

sink(file = "~/Downloads/test.txt")
    print(atest)
    print(tuks)
sink()



 tuks$id[,"p adj"]

str(tuks)

print(tuks, digits = 22)


kstest <- kruskal.test(`Attachment Times`~id, data = spasm_data)

dunn.test::dunn.test(x = spasm_data$`Attachment Times`, spasm_data$id, method = "sidak")
