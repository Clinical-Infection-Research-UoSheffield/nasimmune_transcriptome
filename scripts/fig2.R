library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(corrplot)
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
Nasalpathways <- read.delim("./data/NES-Nasal_pathwaySelection_less50Missing24ThushanFiltred.tsv")
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
## Fig 2A - Blood heatmap
####
col_fun <- colorRamp2(c(-2,0,2), c("blue","white" ,"red"))
BLDdata <- BLDdata[intersect(Bloodpathways, row.names(BLDdata)),]
new_names <- substr(rownames(BLDdata), 0, nchar(rownames(BLDdata)) - 11)
new_names <- str_to_title(new_names)
new_names[17] <- "Positive Regulation Of Transcription Factor Import"
new_names[13] <- "Negative Regulation Of Apoptotic Signaling"
new_names[2] <- "Regulation Of Cysteine-Type Endopeptidase Activity"
new_names[24] <- "Protein Targeting To ER" 

rownames(BLDdata) <- new_names

svg("./output/figs/fig2a.svg", 20,10)
blood_heatmap <- Heatmap(BLDdata,
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
                       show_column_dend = F,
                       show_row_dend = T,
                       row_names_gp = gpar(fontsize = 14),
                       row_split = 2,
                       row_gap = unit(2, "mm"), border = TRUE
)

blood_heatmap
lgd = Legend(col_fun = col_fun, title = "Normalised enrichment score", direction = "horizontal")
draw(lgd, x = unit(0.5, "npc"), y = unit(0.075, "npc"))
dev.off()

####
## Fig 2B - Nasal heatmap
####
NSLdata <- NSLdata[intersect(Nasalpathways, row.names(NSLdata)),]
new_names <- substr(rownames(NSLdata), 0, nchar(rownames(NSLdata)) - 11)
new_names <- str_to_title(new_names)
rownames(NSLdata) <- new_names

## 2  horizontal split heatmap
svg("./output/figs/fig2b.svg", 20,10)
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
dev.off()

####
## Fig 2C - Blood and nose NES correlation
####
shared_pids <- colnames(NSLdata)[colnames(NSLdata) %in% colnames(BLDdata)]
shared_pathways <- row.names(NSLdata)[row.names(NSLdata) %in% row.names(BLDdata)]
shared_nsl <- NSLdata[shared_pathways, shared_pids]
shared_bld <- BLDdata[shared_pathways, shared_pids]

nes_correlation <- 
    lapply(shared_pathways, function(x){
      
        lapply(shared_pathways, function(y){
            pathway_df <- bind_rows(shared_bld[x,], shared_nsl[y,]) 
            pathway_cor <- cor.test(unlist(as.vector(pathway_df[1,])), unlist(as.vector(pathway_df[2,])), method = "spearman")
            tibble(pathway_bld = x, 
                   pathway_nsl = y,
                   cor = pathway_cor$estimate, 
                   p = pathway_cor$p.value)
        }) %>% bind_rows()
        
    }) %>% bind_rows()

cor_rho_mat <- 
    lapply(shared_pathways, function(x){
        row <- 
            nes_correlation %>% 
            filter(pathway_bld == x) %>%
            select(cor) %>% t()
        tibble(row)
    }) %>% bind_rows() %>%
    as.matrix()

row.names(cor_rho_mat) <- shared_pathways    
colnames(cor_rho_mat) <- shared_pathways    

cor_p_mat <- 
  lapply(shared_pathways, function(x){
    row <- 
      nes_correlation %>% 
      filter(pathway_bld == x) %>%
      select(p) %>% t()
    tibble(row)
  }) %>% bind_rows() %>%
  as.matrix()

row.names(cor_p_mat) <- shared_pathways    
colnames(cor_p_mat) <- shared_pathways 

custom_col <- colorRampPalette(c("blue", "white", "red"))


svg("./output/figs/fig2c.svg", 20,10)
corrplot(cor_rho_mat, p.mat = cor_p_mat, 
         sig.level = 0.05, 
         insig = "blank", 
         method = "circle",
         col = custom_col(200),
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 1.2,
         cl.cex = 1.2)

n <- ncol(cor_rho_mat)
coords <- seq(from = n, to = 1)  # because corrplot plots from top to bottom

for (i in 1:n) {
  rect(xleft = i - 0.5, 
       ybottom = coords[i] - 0.5,
       xright = i + 0.5, 
       ytop = coords[i] + 0.5,
       border = "black", 
       lwd = 3)
}

rect(
  xleft = 0.5, 
  ybottom = 0.5, 
  xright = 3.5, 
  ytop = n + 0.5,
  border = "blue", 
  lwd = 3
)

dev.off()


