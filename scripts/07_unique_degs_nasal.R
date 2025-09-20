library(dplyr)

####
## Functions
####
compare_years_nasal <- 
  function(file2017, file2018, sig_level, fc_lower, fc_upper, output_filename){
    
    if(grepl("Nasal2017", file2017) == T){
      degs2017 <- 
        read.delim(file2017) %>%
        filter(P.Value < sig_level,
               logFC > fc_upper | logFC > fc_lower)
      degs2017 <-
        tibble(
          "geneID" = row.names(degs2017),
          "log2FC_2017" = degs2017$logFC,
          "p_2017"= degs2017$P.Value,
          "padj_2017" = degs2017$adj.P.Val
        )
    }
    
    if(grepl("Nasal2017", file2017) != T){
      degs2017 <- 
        read.delim(file2017) %>%
        filter(pvalue < sig_level,
               log2FoldChange > fc_upper | log2FoldChange < fc_lower) %>%
        select(geneID,
               log2FoldChange,
               pvalue,
               padj)
      
      names(degs2017) <- c("geneID", "log2FC_2017", "p_2017", "padj_2017")
    }
    
    
    print(paste0("Number of 2017 DEGs: ", nrow(degs2017)))
    
    degs2018 <- 
      read.delim(file2018) %>%
      filter(pvalue < sig_level,
             log2FoldChange > fc_upper | log2FoldChange < fc_lower) %>%
      select(geneID,
             log2FoldChange,
             pvalue,
             padj)
    
    names(degs2018) <- c("geneID", "log2FC_2018", "p_2018", "padj_2018")
    print(paste0("Number of 2018 DEGs: ", nrow(degs2018)))
    
    shared_degs <- 
      degs2017$geneID[degs2017$geneID %in% degs2018$geneID]
    
    shared_df <- 
      left_join(
        filter(degs2017, geneID %in% shared_degs),
        filter(degs2018, geneID %in% shared_degs),
        by = "geneID"
      ) %>%
      select(
        "geneID", 
        "log2FC_2017",
        "log2FC_2018",
        "p_2017",
        "p_2018", 
        "padj_2017",
        "padj_2018"
      ) %>%
      mutate(
        avg_log2FC = (log2FC_2017 + log2FC_2018)/2,
        avg_padj= (p_2017 + p_2018) /2,
        up_down_reg = ifelse(avg_log2FC > 0, "up", "down")
      )
    
    
    print(paste0("Number of shared DEGs: ", nrow(shared_df)))
    write.csv(shared_df, 
              output_filename, 
              quote = F, row.names = F)
    shared_df
  }


get_unique_degs <-
  function(nonresponder_degs, responder_degs){
    n_shared_degs <-table(nonresponder_degs$geneID %in% responder_degs$geneID)["TRUE"]
    n_unique_nonresponder <- nrow(nonresponder_degs) - n_shared_degs
    n_unique_responder <- nrow(responder_degs) - n_shared_degs
    shared_degs <- nonresponder_degs$geneID[nonresponder_degs$geneID %in% responder_degs$geneID]
    unique_nonresponder <- nonresponder_degs[!(nonresponder_degs$geneID %in% responder_degs$geneID),]
    unique_responder <-responder_degs[!(responder_degs$geneID %in% nonresponder_degs$geneID),]
    
    
    print(paste0("Number of shared DEGS: ", n_shared_degs))
    print(paste0("Number of unique non-responder DEGS: ", n_unique_nonresponder))
    print(paste0("Number of unique non-responder Up-DEGS: ", table(unique_nonresponder$up_down_reg)["up"]))
    print(paste0("Number of unique non-responder Down-DEGS: ", table(unique_nonresponder$up_down_reg)["down"]))
    print(paste0("Number of unique responder DEGS: ", n_unique_responder))
    print(paste0("Number of unique responder Up-DEGS: ", table(unique_responder$up_down_reg)["up"]))
    print(paste0("Number of unique responder Down-DEGS: ", table(unique_responder$up_down_reg)["down"]))
    
    list(
      shared_degs = shared_degs,
      unique_nonresponders = unique_nonresponder,
      unique_responders = unique_responder
    )
  }


####
## Get DEG shared between years
####

##Blood HAI
nasal_hai <-
  compare_years_nasal("./output/processed_data/degs/Nasal2017_HIGH-HR_4x_max_gmfr_Limma.tsv",
                      "./output/processed_data/degs/Nasal2018_HIGH-HR_4x_max_gmfr_DESEq2.tsv",
                      0.05,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/nasal_hai_shared.csv")

nasal_nohai <-
  compare_years_nasal("./output/processed_data/degs/Nasal2017_LOW-HR_4x_max_gmfr_Limma.tsv",
                      "./output/processed_data/degs/Nasal2018_LOW-HR_4x_max_gmfr_DESEq2.tsv",
                      0.05,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/nasal_nohai_shared.csv")

##Blood IgA
nasal_iga <-
  compare_years_nasal("./output/processed_data/degs/Nasal2017_HIGH-HR_2x_max_iga_fc_Limma.tsv",
                      "./output/processed_data/degs/Nasal2018_HIGH-HR_2x_max_iga_fc_DESEq2.tsv",
                      0.05,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/nasal_iga_shared.csv")

nasal_noiga <-
  compare_years_nasal("./output/processed_data/degs/Nasal2017_LOW-HR_2x_max_iga_fc_Limma.tsv",
                      "./output/processed_data/degs/Nasal2018_LOW-HR_2x_max_iga_fc_DESEq2.tsv",
                      0.05,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/nasal_noiga_shared.csv")

##CD4
nasal_cd4 <-
  compare_years_nasal("./output/processed_data/degs/Nasal2017_HIGH-HR_2x_max_mnp_cd4_fc_Limma.tsv",
                      "./output/processed_data/degs/Nasal2018_HIGH-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                      0.05,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/nasal_cd4_shared.csv")

nasal_nocd4 <-
  compare_years_nasal("./output/processed_data/degs/Nasal2017_LOW-HR_2x_max_mnp_cd4_fc_Limma.tsv",
                      "./output/processed_data/degs/Nasal2018_LOW-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                      0.05,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/nasal_nocd4_shared.csv")

## CD8
nasal_cd8 <-
  compare_years_nasal("./output/processed_data/degs/Nasal2017_HIGH-HR_2x_max_mnp_cd8_fc_Limma.tsv",
                      "./output/processed_data/degs/Nasal2018_HIGH-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                      0.05,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/nasal_cd8_shared.csv")

nasal_nocd8 <-
  compare_years_nasal("./output/processed_data/degs/Nasal2017_LOW-HR_2x_max_mnp_cd8_fc_Limma.tsv",
                      "./output/processed_data/degs/Nasal2018_LOW-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                      0.05,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/nasal_nocd8_shared.csv")

####
## Unique DEGS
####

##Blood
nasal_hai_degs <- get_unique_degs(nasal_nohai, nasal_hai)
saveRDS(nasal_hai_degs, file = "./output/processed_data/degs/unique/nasal_hai_degs.rds")
nasal_iga_degs <- get_unique_degs(nasal_noiga, nasal_iga)
saveRDS(nasal_iga_degs, file = "./output/processed_data/degs/unique/nasal_iga_degs.rds")
nasal_cd4_degs <- get_unique_degs(nasal_nocd4, nasal_cd4)
saveRDS(nasal_cd4_degs, file = "./output/processed_data/degs/unique/nasal_cd4_degs.rds")
nasal_cd8_degs <- get_unique_degs(nasal_nocd8, nasal_cd8)
saveRDS(nasal_cd8_degs, file = "./output/processed_data/degs/unique/nasal_cd8_degs.rds")
