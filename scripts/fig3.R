library(dplyr)
library(ggplot2)
library(ggvenn)
library(ggrepel)
library(ggpubr)
library(enrichR)
library(stringr)
library(GO.db)
library(scales)

####
## Functions
####
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
## read in DEGS
####
blood_hai <- readRDS("./output/processed_data/degs/unique/blood_hai_degs.rds")
blood_hai_responder_up <- blood_hai$unique_responders$geneID[blood_hai$unique_responders$up_down == "Up"]
blood_hai_responder_down <- blood_hai$unique_responders$geneID[blood_hai$unique_responders$up_down == "Down"]


blood_iga <- readRDS("./output/processed_data/degs/unique/blood_iga_degs.rds")
blood_iga_responder_up <- blood_iga$unique_responders$geneID[blood_iga$unique_responders$up_down == "Up"]
blood_iga_responder_down <- blood_iga$unique_responders$geneID[blood_iga$unique_responders$up_down == "Down"]

blood_cd4 <- readRDS("./output/processed_data/degs/unique/blood_cd4_degs.rds")
blood_cd4_responder_up <- blood_cd4$unique_responders$geneID[blood_cd4$unique_responders$up_down == "Up"]
blood_cd4_responder_down <- blood_cd4$unique_responders$geneID[blood_cd4$unique_responders$up_down == "Down"]

blood_cd8 <- readRDS("./output/processed_data/degs/unique/blood_cd8_degs.rds")
blood_cd8_responder_up <- blood_cd8$unique_responders$geneID[blood_cd8$unique_responders$up_down == "Up"]
blood_cd8_responder_down <- blood_cd8$unique_responders$geneID[blood_cd8$unique_responders$up_down == "Down"]

blood_hai <- readRDS("./output/processed_data/degs/unique/blood_hai_degs.rds")
blood_hai_responder_up <- blood_hai$unique_responders$geneID[blood_hai$unique_responders$up_down == "Up"]
blood_hai_responder_down <- blood_hai$unique_responders$geneID[blood_hai$unique_responders$up_down == "Down"]

nasal_hai <- readRDS("./output/processed_data/degs/unique/nasal_hai_degs.rds")
nasal_hai_responder_up <- nasal_hai$unique_responders$geneID[nasal_hai$unique_responders$up_down == "Up"]
nasal_hai_responder_down <- nasal_hai$unique_responders$geneID[nasal_hai$unique_responders$up_down == "Down"]

nasal_iga <- readRDS("./output/processed_data/degs/unique/nasal_iga_degs.rds")
nasal_iga_responder_up <- nasal_iga$unique_responders$geneID[nasal_iga$unique_responders$up_down == "Up"]
nasal_iga_responder_down <- nasal_iga$unique_responders$geneID[nasal_iga$unique_responders$up_down == "Down"]

nasal_cd4 <- readRDS("./output/processed_data/degs/unique/nasal_cd4_degs.rds")
nasal_cd4_responder_up <- nasal_cd4$unique_responders$geneID[nasal_cd4$unique_responders$up_down == "Up"]
nasal_cd4_responder_down <- nasal_cd4$unique_responders$geneID[nasal_cd4$unique_responders$up_down == "Down"]

nasal_cd8 <- readRDS("./output/processed_data/degs/unique/nasal_cd8_degs.rds")
nasal_cd8_responder_up <- nasal_cd8$unique_responders$geneID[nasal_cd8$unique_responders$up_down == "Up"]
nasal_cd8_responder_down <- nasal_cd8$unique_responders$geneID[nasal_cd8$unique_responders$up_down == "Down"]

####
## Read in rlog values
####
blood_2017_rlog <- read.table("./output/processed_data/Blood_2017_Rlog_normalized.tsv", header = T)
blood_2018_rlog <- read.table("./output/processed_data/Blood_2018_Rlog_normalized.tsv", header = T)
nasal_2018_rlog <- read.table("./output/processed_data/Nasal_2018_Rlog_normalized.tsv", header = T)

####
## Read in metadata
####
metadata <- readRDS("./output/processed_data/metadata_cleaned.rds")

####
## Read in DESEQ
####
blood_hai_responder <- read.table("./output/processed_data/degs/Blood2018_HIGH-HR_4x_max_gmfr_DESEq2.tsv", header = T)

####
## Fig 3b
####
response_venn <- ggvenn(
  list("HAI" = blood_hai_responder_up, 
       "IgA" = blood_iga_responder_up,
       "CD4" = blood_cd4_responder_up,
       "CD8" = blood_cd8_responder_up),
  show_percentage = F
) 

response_venn
response_venn$layers[[3]]$data$x <- c(-1.5, -1.3, 1.3, 1.5)
response_venn$layers[[3]]$data$y <- c(-0.9, 0.95, 0.95, -0.9)

ggsave("./output/figs/fig3b1.svg")
ggsave("./output/figs/fig3b1.png")

blood_non_hai_degs <- c(blood_iga_responder_up, blood_cd4_responder_up, blood_cd8_responder_up)
blood_unique_up <- blood_hai_responder_up[!(blood_hai_responder_up %in% blood_non_hai_degs)]

#save file to use in cytoscape
write.csv(blood_unique_up, "./output/processed_data/degs/unique/blood_hai_unique_up.csv", row.names = F)

enriched <- enrichr(blood_unique_up, "GO_Biological_Process_2023")

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
enriched_filtered$Term[4] <- "Negative Regulation of Type I Interferon" 

plotEnrich(enriched_filtered, 
           showTerms = 8, 
           numChar = 60, 
           y = "Count", 
           orderBy = "FDR",
           title = "",
           x = "") +
           guides(fill=guide_colorbar(title="P (Adjusted)", reverse = T))

ggsave("./output/figs/fig3b2.svg", width = 8, height = 3)
ggsave("./output/figs/fig3b2.png", width = 8, height = 3)

####
## Fig 3c
####

response_venn <- ggvenn(
  list("HAI" = nasal_hai_responder_up, 
       "IgA" = nasal_iga_responder_up,
       "CD4" = nasal_cd4_responder_up,
       "CD8" = nasal_cd8_responder_up),
  show_percentage = F
) 

response_venn
response_venn$layers[[3]]$data$x <- c(-1.5, -1.3, 1.3, 1.5)
response_venn$layers[[3]]$data$y <- c(-0.9, 0.95, 0.95, -0.9)

ggsave("./output/figs/fig3c1.svg")
ggsave("./output/figs/fig3c1.png")

# enrichr
nasal_non_cd4hai_degs <- c(nasal_iga_responder_up, nasal_cd8_responder_up)
nasal_unique_up <- nasal_cd4_responder_up[!(nasal_cd4_responder_up %in% nasal_non_cd4hai_degs)]

#save file to use in cytoscape
write.csv(nasal_unique_up, "./output/processed_data/degs/unique/nasal_haicd4_unique_up.csv", row.names = F)

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
           showTerms = 6, 
           numChar = 60, 
           y = "Count", 
           orderBy = "FDR",
           title = "",
           x = "") +
  guides(fill=guide_colorbar(title="P (Adjusted)", reverse = T))

ggsave("./output/figs/fig3c2.svg", width = 6, height = 3)

####
## Fig 3f
####

nasal_non_cd4hai_degs <- c(nasal_iga_responder_up, nasal_cd8_responder_up)
nasal_unique_up <- nasal_cd4_responder_up[!(nasal_cd4_responder_up %in% nasal_non_cd4hai_degs)]

#save file to use in cytoscape
write.csv(nasal_unique_up, "./output/processed_data/degs/unique/nasal_haicd4_unique_up.csv", row.names = F)


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
           showTerms = 6, 
           numChar = 60, 
           y = "Count", 
           orderBy = "FDR",
           title = "",
           x = "") +
  guides(fill=guide_colorbar(title="P (Adjusted)", reverse = T))

ggsave("./output/figs/fig3c2.svg", width = 6, height = 3)




####
## Fig 3d
####
blood_rlog <- 
    bind_cols(
        blood_2017_rlog[c("IFIT1","HERC5"),],
        blood_2018_rlog[c("IFIT1","HERC5"),]
    )

blood_rlog <- as.data.frame(t(blood_rlog))
blood_rlog$subid1 <- substr(rownames(blood_rlog), 1, 5)
blood_rlog$visit <- substr(rownames(blood_rlog), 7, 8)
blood_rlog <- 
    left_join(blood_rlog, metadata[,c("subid1", "HR_4x_max_gmfr")])

blood_rlog <- blood_rlog %>%
    mutate(HR_4x_max_gmfr = case_when(
    HR_4x_max_gmfr == 0 ~ "Non-Responders",
    HR_4x_max_gmfr == 1 ~ "Responders"
  )) %>%
    mutate(visit = case_when(
      visit == "V0" ~ "0",
      visit == "V2" ~ "2"
    ))
  
blood_rlog$HR_4x_max_gmfr <- 
    factor(blood_rlog$HR_4x_max_gmfr, levels = c("Non-Responders", "Responders"))

blood_rlog %>%
    ggplot(aes(x = visit, y = IFIT1)) +
    geom_boxplot(aes(fill = HR_4x_max_gmfr), outliers = F) +
    geom_jitter(width = 0.1, size = 0.6, alpha = 0.5) +
    stat_compare_means(method = "t.test", paired = T, label = "p.format", 
                       label.x = 1.25, label.y = 15.1) +
    facet_wrap(. ~ HR_4x_max_gmfr) +
    ylab(expression(italic(IFIT1) ~ "Expression (rlog normalised)")) +
    xlab("Days Post Vaccine") +
    scale_fill_manual(values = c("Non-Responders" = "#0C0C86" , 
                               "Responders" = "#0C860C")) +
    theme_classic2() +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 10),
          axis.title.x = element_text(size = 10))

ggsave("./output/figs/fig3d.svg", width = 4)

####
## Fig 3e
####
blood_rlog %>%
  ggplot(aes(x = visit, y = HERC5)) +
  geom_boxplot(aes(fill = HR_4x_max_gmfr), outliers = F) +
  geom_jitter(width = 0.1, size = 0.6, alpha = 0.5) +
  stat_compare_means(method = "t.test", paired = T, label = "p.format", 
                     label.x = 1.25, label.y = 13.8) +
  facet_wrap(. ~ HR_4x_max_gmfr) +
  ylab(expression(italic(HERC5) ~ "Expression (rlog normalised)")) +  
  xlab("Days Post Vaccine") +
  scale_fill_manual(values = c("Non-Responders" = "#0C0C86" , 
                               "Responders" = "#0C860C")) +
  theme_classic2() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10))
ggsave("./output/figs/fig3e.svg", width = 4)

####
## Fig 3k
####
nasal_rlog <- nasal_2018_rlog[c("IFIT1","HELZ2"),]

nasal_rlog <- as.data.frame(t(nasal_rlog))
nasal_rlog$subid1 <- substr(rownames(nasal_rlog), 1, 5)
nasal_rlog$visit <- substr(rownames(nasal_rlog), 7, 8)
nasal_rlog <- 
  left_join(nasal_rlog, metadata[,c("subid1", "HR_4x_max_gmfr")])

nasal_rlog <- nasal_rlog %>%
  mutate(HR_4x_max_gmfr = case_when(
    HR_4x_max_gmfr == 0 ~ "Non-Responders",
    HR_4x_max_gmfr == 1 ~ "Responders"
  )) %>%
  mutate(visit = case_when(
    visit == "V0" ~ "0",
    visit == "V2" ~ "2"
  ))

nasal_rlog$HR_4x_max_gmfr <- 
  factor(nasal_rlog$HR_4x_max_gmfr, levels = c("Non-Responders", "Responders"))

nasal_rlog %>%
  ggplot(aes(x = visit, y = HELZ2)) +
  geom_boxplot(aes(fill = HR_4x_max_gmfr), outliers = F) +
  geom_jitter(width = 0.1, size = 0.6, alpha = 0.5) +
  stat_compare_means(method = "t.test", paired = T, label = "p.format", 
                     label.x = 1.25, label.y = 10.8) +
  facet_wrap(. ~ HR_4x_max_gmfr) +
  ylab(expression(italic(MX1) ~ "Expression (rlog normalised)")) + 
  xlab("Days Post Vaccine") +
  scale_fill_manual(values = c("Non-Responders" = "#0C0C86" , 
                               "Responders" = "#0C860C")) +
  theme_classic2() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10))
ggsave("./output/figs/fig3f.svg", width = 4)


####
## Fig 3k
####
nasal_rlog %>%
  ggplot(aes(x = visit, y = HERC5)) +
  geom_boxplot(aes(fill = HR_4x_max_gmfr), outliers = F) +
  geom_jitter(width = 0.1, size = 0.6, alpha = 0.5) +
  stat_compare_means(method = "t.test", paired = T, label = "p.format", 
                     label.x = 1.25, label.y = 9.8) +
  facet_wrap(. ~ HR_4x_max_gmfr) +
  ylab(expression(italic(HERC5) ~ "Expression (rlog normalised)")) + 
  xlab("Days Post Vaccine") +
  scale_fill_manual(values = c("Non-Responders" = "#0C0C86" , 
                               "Responders" = "#0C860C")) +
  theme_classic2() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10))
ggsave("./output/figs/fig3g.svg", width = 4)



