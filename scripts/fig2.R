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
## Functions
####
reorder_infile <- function(fc_file){
  output_df <- read.delim(fc_file)
  output_df$gene <- row.names(output_df)
  row.names(output_df) <- NULL
  output_df <- output_df %>% dplyr::select("gene", everything())
  return(output_df)
}

get_all_children <- function(go_id, mapping_list) {
  # Get direct children in ontology
  children <- mapping_list[[go_id]]
  if (is.null(children)) {
    return(character(0))
  }
  
  # Recursively get children of children
  all_children <- unlist(lapply(children, get_all_children, mapping_list = mapping_list))
  return(unique(c(children, all_children)))
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

# Gene foldchange data
FClist <- list()
bld17 <- reorder_infile("./output/processed_data/Blood_2017_rlog_log2FC.tsv")
bld18 <- reorder_infile("./output/processed_data/Blood_2018_rlog_log2FC.tsv")
nsl17 <- reorder_infile("./output/processed_data/Nasal_2017_FiltredToSharedGenes_log2FC_Ranks.tsv")
nsl18 <- reorder_infile("./output/processed_data/Nasal_2018_rlog_log2FC.tsv")
FClist$BLD <- full_join(bld17, bld18, by = "gene")
FClist$NSL <- full_join(nsl17, nsl18, by="gene")

metadata <- read.csv("./data/metadata.csv", header = T)

####
## Clean metadata
####
metadata <- 
  metadata %>%
  dplyr::select("subid1", contains("v2_shed")) %>%
  dplyr::rename("subid" = subid1)

metadata$n_shed <- rowSums(metadata[,c("h1_v2_shed", "h3_v2_shed", "b_v2_shed")])


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

####
## Fig2c
####
shed_degs <- read.csv("./output/processed_data/degs/meta/nasal_shedsignif_genes.csv")

combined_nsl_degs <- read.csv("./output/processed_data/degs/meta/nasal_shedall_genes.csv")

unique_degs <- shed_degs$geneID
combined_nsl_degs$unique_deg <- combined_nsl_degs$geneID %in% unique_degs
combined_nsl_degs$p_adj_random <- as.numeric(combined_nsl_degs$p_adj_random)

unique_degs_df <-
  combined_nsl_degs %>%
  dplyr::filter(unique_deg == T) %>%
  dplyr::select("geneID", "eff_meta_random", "p_adj_random") 
names(unique_degs_df) <- c("geneID", "avg_fc", "avg_p")

combined_nsl_degs$minuslog10_p <- -log10(combined_nsl_degs$p_adj_random)

# obtained from enrichcr (see fig 4b below)
enriched_degs <- c(
  "ISG20", "IFITM1", "OAS1", "OAS3", "MX1", "IFIT5",
  "ISG15", "IFIT1", "APOBEC3A", "STAT2", "IFIT3", "IFIT2",
  "IL1B", "IRF7", "LCP1", "ZFP36", "SBNO2", "CXCL11",
  "LILRB2", "CSF1", "GPR183", "PRKCB", "LPIN1", "FFAR2",
  "TNFAIP6", "SLC7A5", "CXCR4", "MX2", "NFKBIA", "IFI6",
  "CCL20", "IRF8", "TNFSF13B", "RASGRP3"
)

combined_nsl_degs$unique_deg <- combined_nsl_degs$geneID %in% enriched_degs

combined_nsl_degs$select_label = ifelse(combined_nsl_degs$unique_deg == F, "", combined_nsl_degs$geneID)

ggplot() +
  geom_point(data = dplyr::filter(combined_nsl_degs, geneID %in% unique_degs_df$geneID), 
             aes(x = eff_meta_random, y = minuslog10_p), 
             color = "red",
             size = 0.1) +
  geom_point(data = dplyr::filter(combined_nsl_degs, !(geneID %in% unique_degs_df$geneID)), 
             aes(x = eff_meta_random, y = minuslog10_p), 
             color = "black",
             size = 0.06,
             alpha = 0.25) +
  geom_vline(xintercept = 0.322, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = -0.322, linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.3) +
  annotate("text", x = 1, y = 2.5, size = 4, colour = "red", label = "Up = 67", fontface = "bold") +
  annotate("text", x = -0.55, y = 2.5, size = 4, colour = "red", label = "Down = 2", fontface = "bold") +
  geom_text_repel(data = dplyr::filter(combined_nsl_degs, geneID %in% c("AADAT", "DNAJC28")), 
                  aes(x = eff_meta_random, y = minuslog10_p, label = geneID),
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  size = 3,
                  color = "red",
                  fontface = "bold") +
  theme_classic() +
  theme(legend.position = "none") +
  xlab(expression(Log[2]~Fold~Change)) +
  ylab(expression(-log[10](padj))) +
  xlim(-0.75, 1.3)

ggsave("./output/figs/fig4a.svg", height = 4, width = 4.5)
ggsave("./output/figs/fig4a.png", height = 4, width = 4.5)

ggplot() +
  geom_point(data = dplyr::filter(combined_nsl_degs, geneID %in% unique_degs_df$geneID), 
             aes(x = eff_meta_random, y = minuslog10_p), 
             color = "red",
             size = 0.06) +
  geom_point(data = dplyr::filter(combined_nsl_degs, !(geneID %in% unique_degs_df$geneID)), 
             aes(x = eff_meta_random, y = minuslog10_p), 
             color = "black",
             size = 0.06,
             alpha = 0.25) +
  geom_text_repel(data = dplyr::filter(combined_nsl_degs, unique_deg), 
                  aes(x = eff_meta_random, y = minuslog10_p, label = geneID),
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  size = 2.75,
                  color = "red",
                  fontface = "bold",
                  text.padding = 1) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.title = element_blank()) +
  xlim(0.322, 1) +
  ylim(1, 3.75)

ggsave("./output/figs/fig4a_zoom.svg", height = 4, width = 5)
ggsave("./output/figs/fig4a_zoom.png", height = 4, width = 5)

####
## Fig 4b
####
nasal_unique_up <- unique_degs_df$geneID[unique_degs_df$geneID != "AADAT"]

enriched <- enrichr(nasal_unique_up, "GO_Biological_Process_2023")

## Filter GO to include the lowest node
signif_go <- enriched$GO_Biological_Process_2023$Term[enriched$GO_Biological_Process_2023$Adjusted.P.value < 0.05]

go_ids <- toupper(gsub(".*(GO:\\d+).*", "\\1", signif_go))

bp_children_list <- as.list(GOBPCHILDREN)

go_links <- 
  lapply(go_ids, 
         function(x){
           all_descendants <- get_all_children(x, bp_children_list)
           descendent_n <- all_descendants[all_descendants %in% go_ids]
           tibble(GO = x, descendent_n = length(descendent_n))
         }) %>%
  bind_rows()

terminal_go <- go_links$GO[go_links$descendent_n == 0]

combined_pattern <- paste(substr(terminal_go, 4, 10), collapse = "|")
terminal_names <- enriched$GO_Biological_Process_2023$Term[grepl(combined_pattern, enriched$GO_Biological_Process_2023$Term)]

## filter enrichr output
enriched_filtered <- enriched$GO_Biological_Process_2023[enriched$GO_Biological_Process_2023$Term %in% terminal_names,]
row.names(enriched_filtered) <- seq(1, nrow(enriched_filtered))

enriched_filtered$Term <- substr(enriched_filtered$Term, 0, nchar(enriched_filtered$Term) - 12)
enriched_filtered[[1]][3] <- "Regulation Of Monocyte Chemotactic Protein-1 Production"
enriched_filtered[[1]][4] <- "Regulation Of CXCL2 Production"

plotEnrich(enriched_filtered, 
           showTerms = 7, 
           numChar = 70, 
           y = "Count", 
           orderBy = "FDR",
           title = "",
           x = "") +
  guides(fill=guide_colorbar(title="P (Adjusted)", reverse = T))

ggsave("./output/figs/fig4b.svg", height = 4, width = 8)
ggsave("./output/figs/fig4b.png", height = 4, width = 8)

all_enriched_genes <- unique(unlist(strsplit(enriched_filtered$Genes, ";")))

