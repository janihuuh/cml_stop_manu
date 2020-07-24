
me=system("whoami", intern = T)
setwd(paste0("/Users/", me, "/Dropbox/cml_stop_clean/"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
library(data.table)
library(gridExtra)
library(RColorBrewer)

theme_set(theme_classic(base_size = 12))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))

add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))
facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))


## Run all fun_* codes
for(code in list.files("src/R/", "fun", full.names = T, recursive = T)){

  message(code)
  source(code)

}

