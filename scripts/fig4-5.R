## Show cell type change over time for blood and nose and correlate FC with IFN stimulated signaure

library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape)

vst <- read.delim("./output/processed_data/Nasal_2018_vst_normalized.tsv")
vst$geneID <- row.names(vst)
markers <- read.delim("data/nasal_cell_markers.tsv", header = F)

