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
## Fig 3A - Shedding
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

# obtained from enrichcr (see fig 3b below)
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

ggsave("./output/figs/fig3a.svg", height = 4, width = 4.8)

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

ggsave("./output/figs/fig3a_zoom.svg", height = 4, width = 4.8)


####
## Fig 3b
####
nasal_unique_up <- unique_degs_df$geneID[!(unique_degs_df$geneID %in% c("AADAT", "DNAJC28"))]

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

plotEnrich(enriched_filtered, 
           showTerms = 15, 
           numChar = 70, 
           y = "Count", 
           orderBy = "FDR",
           title = "",
           x = "") +
  guides(fill=guide_colorbar(title="P (Adjusted)", reverse = T)) +
  scale_y_continuous(breaks = seq(0,9))

ggsave("./output/figs/fig3b.svg", height = 4, width = 8)

all_enriched_genes <- unique(unlist(strsplit(enriched_filtered$Genes, ";")))

####
## Fig3c
nsl_hai_genes <- c("IFIT1", "IRF4", "HERC5", "ZBP1", "DHX58", "MX1", "IFI6", "APOBEC3A", "IFIT1")
nasal_fc <- 
  FClist$NSL %>% 
  dplyr::filter(gene %in% nsl_hai_genes) %>%
  t() %>%
  as.data.frame()

colnames(nasal_fc) <- nasal_fc[1,]
nasal_fc$subid <- rownames(nasal_fc)
nasal_fc <- nasal_fc[-1,]
rownames(nasal_fc) <- NULL
nasal_fc <- nasal_fc %>% dplyr::select("subid", everything())

nasal_shedding <- 
  metadata %>%
  dplyr::filter(subid %in% nasal_fc$subid) %>%
  dplyr::select("subid", "n_shed") %>%
  left_join(nasal_fc, by= "subid")

nasal_shedding$n_shed <- as.factor(nasal_shedding$n_shed)
nasal_shedding$shedding <- nasal_shedding$n_shed != 0
nasal_shedding$MX1 <- as.numeric(nasal_shedding$MX1)
nasal_shedding$IFIT1 <- as.numeric(nasal_shedding$IFIT1)
nasal_shedding$HERC5 <- as.numeric(nasal_shedding$HERC5)

my_comparisons <- list( c("0", "1"), c("0", "2"), c("0", "3") )

stat_test <- nasal_shedding %>% 
  wilcox_test(IFIT1 ~ n_shed, comparisons = my_comparisons) %>%
  add_xy_position(x = "n_shed") %>%
  mutate(myformatted.p = paste0("p = ", p))

lm_test <- lm(IFIT1 ~ as.numeric(n_shed), data = nasal_shedding)
summary(lm_test)

global_p <- summary(lm_test)$coefficients[2, 4] 
formatted_global_p <- paste0("Linear trend p = ", signif(global_p, 3))

nasal_shedding %>%
  ggplot(aes(x = n_shed, y = IFIT1)) +
  geom_boxplot(outlier.shape = NA, colour = "darkblue") +
  geom_jitter(width = 0.1, size = 0.75, colour = "darkblue") +
  annotate("text", x = 1.75, y = 5.75, label = formatted_global_p, size = 4, hjust = 0) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  xlab("Number of strains shed") +
  ylab(expression(italic(IFIT1) ~ "expression (" * Log[2] * " fold change)"))

ggsave("./output/figs/fig3c.svg", width = 4.8, height = 4)


####
## Blood
####
nsl_hai_genes <- c("IFIT1", "IRF4", "HERC5", "ZBP1", "DHX58", "MX1", "IFI6", "APOBEC3A", "IFIT1")
blood_fc <- 
  FClist$BLD %>% 
  dplyr::filter(gene %in% nsl_hai_genes) %>%
  t() %>%
  as.data.frame()

colnames(blood_fc) <- blood_fc[1,]
blood_fc$subid <- rownames(blood_fc)
blood_fc <- blood_fc[-1,]
rownames(blood_fc) <- NULL
blood_fc <- blood_fc %>% dplyr::select("subid", everything())

blood_shedding <- 
  metadata %>%
  dplyr::filter(subid %in% blood_fc$subid) %>%
  dplyr::select("subid", "n_shed") %>%
  left_join(blood_fc, by= "subid")

blood_shedding$n_shed <- as.factor(blood_shedding$n_shed)
blood_shedding$shedding <- blood_shedding$n_shed != 0
blood_shedding$MX1 <- as.numeric(blood_shedding$MX1)
blood_shedding$IFIT1 <- as.numeric(blood_shedding$IFIT1)
blood_shedding$HERC5 <- as.numeric(blood_shedding$HERC5)

my_comparisons <- list( c("0", "1"), c("0", "2"), c("0", "3") )

stat_test <- nblood_shedding %>% 
  wilcox_test(IFIT1 ~ n_shed, comparisons = my_comparisons) %>%
  add_xy_position(x = "n_shed") %>%
  mutate(myformatted.p = paste0("p = ", p))

lm_test <- lm(IFIT1 ~ as.numeric(n_shed), data = blood_shedding)
summary(lm_test)

global_p <- summary(lm_test)$coefficients[2, 4] 
formatted_global_p <- paste0("Linear trend p = ", signif(global_p, 3))

blood_shedding %>%
  ggplot(aes(x = n_shed, y = IFIT1)) +
  geom_boxplot(outlier.shape = NA, colour = "darkblue") +
  geom_jitter(width = 0.1, size = 0.75, colour = "darkblue") +
  annotate("text", x = 1.75, y = 5.75, label = formatted_global_p, size = 4, hjust = 0) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  xlab("Number of strains shed") +
  ylab(expression(italic(IFIT1) ~ "expression (" * Log[2] * " fold change)"))

ggsave("./output/figs/fig3d.svg", width = 4.8, height = 4)


