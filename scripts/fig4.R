library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(corrplot)
set.seed(123)

####
## Functions
####
reorder_infile <- function(fc_file){
  output_df <- read.delim(fc_file)
  output_df$gene <- row.names(output_df)
  row.names(output_df) <- NULL
  output_df <- output_df %>% dplyr::select("gene", everything())
  return(output_df)
}

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

# Gene foldchange data
FClist <- list()
bld17 <- reorder_infile("./output/processed_data/Blood_2017_rlog_log2FC.tsv")
bld18 <- reorder_infile("./output/processed_data/Blood_2018_rlog_log2FC.tsv")
nsl17 <- reorder_infile("./output/processed_data/Nasal_2017_FiltredToSharedGenes_log2FC_Ranks.tsv")
nsl18 <- reorder_infile("./output/processed_data/Nasal_2018_rlog_log2FC.tsv")
FClist$BLD <- full_join(bld17, bld18, by = "gene")
FClist$NSL <- full_join(nsl17, nsl18, by="gene")

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
NSLdata <- NSLdata[intersect(Nasalpathways, row.names(NSLdata)),]
new_names <- substr(rownames(NSLdata), 0, nchar(rownames(NSLdata)) - 11)
new_names <- str_to_title(new_names)
rownames(NSLdata) <- new_names

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

svg("./output/figs/fig4a.svg", 20,10)
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
## Fig 4B - Blood and nose NES correlation
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


svg("./output/figs/fig4b.svg", 20,10)
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

####
## Fig 3c
####
nose_correlated_genes <- c("MX1", "IL1B")
blood_correlated_genes <- c("MX1", "IL1B", "CSF1")

nasal_fc <- 
  FClist$NSL %>% 
  dplyr::filter(gene %in% nose_correlated_genes) %>%
  t() %>%
  as.data.frame()

colnames(nasal_fc) <- nasal_fc[1,]
nasal_fc$subid <- rownames(nasal_fc)
nasal_fc <- nasal_fc[-1,]
rownames(nasal_fc) <- NULL
nasal_fc <- nasal_fc %>% 
    dplyr::select("subid", everything()) %>%
    dplyr::filter(subid %in% shared_pids) 
colnames(nasal_fc) <- paste0("nose_", colnames(nasal_fc))
nasal_fc <- rename(nasal_fc, "subid" = nose_subid)

blood_fc <- 
  FClist$BLD %>% 
  dplyr::filter(gene %in% blood_correlated_genes) %>%
  t() %>%
  as.data.frame() 

colnames(blood_fc) <- blood_fc[1,]
blood_fc$subid <- rownames(blood_fc)
blood_fc <- blood_fc[-1,]
rownames(blood_fc) <- NULL
blood_fc <- blood_fc %>% 
    dplyr::select("subid", everything()) %>%
    dplyr::filter(subid %in% shared_pids) 
colnames(blood_fc) <- paste0("blood_", colnames(blood_fc))
blood_fc <- rename(blood_fc, "subid" = blood_subid)

cor_df <- left_join(nasal_fc, blood_fc)
cor_df <- cor_df %>% mutate(across(2:6, as.numeric))

ggplot(cor_df, aes(x = nose_MX1, y = blood_MX1)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) + 
    stat_cor(label.y = 4, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    theme_bw()

ggplot(cor_df, aes(x = nose_MX1, y = blood_IL1B)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  stat_cor(label.y = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()

ggplot(cor_df, aes(x = nose_MX1, y = blood_CSF1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  stat_cor(label.y = 1.5, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()

ggplot(cor_df, aes(x = nose_IL1B, y = blood_MX1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  stat_cor(label.y = 2, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()

ggplot(cor_df, aes(x = nose_IL1B, y = blood_CSF1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  stat_cor(label.y = 1, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()

ggplot(cor_df, aes(x = nose_IL1B, y = blood_IL6)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  stat_cor(label.y = 1, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw()

