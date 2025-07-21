library(ggplot2)
library(ggtext)
library(cowplot)
library(tidyr)
library(minpack.lm)
library(data.table)
library(magick)
library(ggbeeswarm)

devtools::load_all("~/git/spasmtools")

theme_set(theme_cowplot(10))

colorz <- c("#666666",
            "#f7b8d8",
            "#ef72b2",
            "#e7298a")


spasm_files <- list("0 &micro;M" = "optical-trap/spasm-analysis/control/spasm-included.xlsx",
                    "1 &micro;M" = "optical-trap/spasm-analysis/1uM-dani/spasm_1uM-dani_10uM-atp.xlsx",
                    "10 &micro;M" = "optical-trap/spasm-analysis/10uM-dani/spasm-included-manual.xlsx",
                    "100 &micro;M" = "optical-trap/spasm-analysis/100uM-dani/spasm_100uM-dani_10uM-atp.xlsx"
                    )

spasm_data <- spasm_read(spasm_files)

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
ggplot(dplyr::filter(spasm_data, id %in% c("0 &micro;M", "10 &micro;M")))+
  stat_ecdf(aes(x = Displacements, color = id), linewidth = 1, alpha = 0.7)+
  stat_ecdf(data = dplyr::filter(rnorm_df, id %in% c("0 &micro;M", "10 &micro;M")),
                                 aes(x = total, color = id), linewidth = 0.5, linetype = "dashed")+
  scale_color_manual(values = c(colorz[1], colorz[3]))+
  coord_cartesian(xlim = c(-15, 20))+
  ylab("Probability")+
  xlab("Total Displacement (nm)")+
  ggtitle("")+
  theme_cowplot(11)+
  theme(
  legend.position = "none"
  )

png("supplemental-figures/step-size-distributions.png", width = 3, height = 3, res = 300, units = "in")
ggstep1
dev.off()


pdf("supplemental-figures/step-size-distributions.pdf", width = 3, height = 3)
ggstep1
dev.off()

## ggstep2 <- ggstep1+coord_cartesian(xlim = c(-10, 15))+theme_cowplot(11)+theme(legend.position = "none")

## ggstep3 <- ggdraw(ggstep2)+draw_plot(ggstep1, 0.6, 0.15, 0.4, 0.4)

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
