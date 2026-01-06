library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(corrplot)
library(tidyr)
library(ggpubr)
library(ggrepel)
library(rstatix)
library(enrichR)
library(GO.db)
library(circlize)
set.seed(123)

####
## Read in data
####

# Normalised enrhichment scores (NES)
NESlist <- list()
NESlist$bld17 <- read.delim("./output/processed_data/gsea/NES_Blood_2017.tsv")
NESlist$bld18 <- read.delim("./output/processed_data/gsea/NES_Blood_2018.tsv")
NESlist$nsl17 <- read.delim("./output/processed_data/gsea/NES_Nasal_2017.tsv")
NESlist$nsl18 <- read.delim("./output/processed_data/gsea/NES_Nasal_2018.tsv")
NESlist$BLD <- full_join(NESlist$bld17, NESlist$bld18, by="X")
NESlist$NSL <- full_join(NESlist$nsl17, NESlist$nsl18, by="X")

## Select genesets. Present in at least 50% of participants
Bloodpathways <- read.delim("./data/NES-Blood_pathwaySelection_less50Missing.tsv")
Bloodpathways <- Bloodpathways$Pathway
Nasalpathways <- read.delim("./data/NES-Nasal_pathwaySelection_less50Missing.tsv")
Nasalpathways <- Nasalpathways$Pathway

## Blood NES filtered by genesets
BLDdata <- NESlist$BLD %>% filter(X %in% Bloodpathways)
row.names(BLDdata) <- BLDdata$X
BLDdata$X <- NULL
BLDdata <- as.matrix(BLDdata)

## Nasal NES filtered by genesets
NSLdata <- NESlist$NSL %>% filter(X %in% Nasalpathways)
row.names(NSLdata) <- NSLdata$X
NSLdata$X <- NULL
NSLdata <- as.matrix(NSLdata)

####
## Fig 2A - Nasal heatmap
####
col_fun <- colorRamp2(c(-2,0,2), c("blue","white" ,"red"))
NSLdata <- NSLdata[intersect(Nasalpathways, row.names(NSLdata)),]
new_names <- substr(rownames(NSLdata), 0, nchar(rownames(NSLdata)) - 11)
new_names <- str_to_title(new_names)
rownames(NSLdata) <- new_names

svg("./output/figs/fig2a.svg", 20,10)
nasal_heatmap <- Heatmap(NSLdata,
                         name = "Heatmap\nName",
                         cluster_rows = TRUE, 
                         cluster_columns = TRUE,
                         col = col_fun, 
                         show_heatmap_legend = F,
                         column_km_repeats = 100,
                         row_km_repeats = 100,
                         heatmap_height = unit(15, "cm"),
                         heatmap_width = unit(30, "cm"),
                         show_column_names = F,
                         show_row_dend = T,
                         show_column_dend = F,
                         row_names_gp = gpar(fontsize = 14),
                         row_split = 2,
                         row_gap = unit(2, "mm"), border = TRUE
)
nasal_heatmap

lgd = Legend(col_fun = col_fun, title = "Normalised enrichment score", direction = "horizontal")
draw(lgd, x = unit(0.5, "npc"), y = unit(0.075, "npc"))
dev.off()

####
## Fig 2B - Nose NES correlation
####
nsl_genesets <- t(NSLdata)

custom_col <- colorRampPalette(c("blue", "white", "red"))

svg("./output/figs/fig2b.svg", 20,10)
corrplot(cor(nsl_genesets), 
         p.mat = cor.mtest(nsl_genesets, method = "spearman", conf.level = 0.95)$p, 
         sig.level = 0.05, 
         insig = "blank", 
         method = "circle",
         col = custom_col(200),
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 1.2,
         cl.cex = 1.2,
         order = "hclust")

dev.off()