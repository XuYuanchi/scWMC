setwd("D:/MyWorkWorld/Imputation/scWMC/DEG/")
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

color = readRDS("D:/MyWorkWorld/Imputation/scWMC/color.rds")
methods_name = c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
                 "scImpute", "scTSSR", "scVI", "VIPER")


# 50
dat = readRDS("Results/blk50.rds")

dat = as.data.frame(dat[1:6, ])

dat$x = c('50', '100', '150', '200', '250', '300')

dat = gather(dat, Methods, y_value, Raw, scWMC, ALRA, DCA, DrImpute, EnImpute, MAGIC, PBLR, SAVER, 
             scImpute, scTSSR, scVI, VIPER)

dat$Methods = factor(dat$Methods, levels = methods_name)
dat$x = factor(dat$x, levels = c('50', '100', '150', '200', '250', '300'))

p1 = ggplot(data = dat, aes(x=x, y=y_value,group=Methods, color=Methods)) +
          # geom_area() +
          scale_color_manual(values = color) + 
          geom_line(size=1) +
          geom_point() +
          xlab("number of selected  differentially expressed genes") +
          ylab("overlap") +
          theme_bw() +
          theme( panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank())



# 100
dat = readRDS("Results/blk100.rds")

dat = as.data.frame(dat[1:6, ])

dat$x = c('50', '100', '150', '200', '250', '300')

dat = gather(dat, Methods, y_value, Raw, scWMC, ALRA, DCA, DrImpute, EnImpute, MAGIC, PBLR, SAVER, 
             scImpute, scTSSR, scVI, VIPER)

dat$Methods = factor(dat$Methods, levels = methods_name)
dat$x = factor(dat$x, levels = c('50', '100', '150', '200', '250', '300'))

p2 = ggplot(data = dat, aes(x=x, y=y_value,group=Methods, color=Methods)) +
  # geom_area() +
  scale_color_manual(values = color) + 
  geom_line(size=1) +
  geom_point() +
  xlab("number of selected  differentially expressed genes") +
  ylab("overlap") +
  theme_bw() +
  theme( panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank())

p = p1 + p2 + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')

ggsave(plot = p, filename = "Results/DEG_sdeg.pdf", width = 8, height = 5)
