library(dplyr)

####
## Functions
####
compare_years_blood <- 
  function(file2017, file2018, sig_level, fc_lower, fc_upper, output_filename){
    
    degs2017 <- 
      read.delim(file2017, header = T) %>%
      dplyr::filter(padj < sig_level,
             log2FoldChange > fc_upper | log2FoldChange < fc_lower) %>%
      dplyr::select(geneID,
             log2FoldChange,
             pvalue,
             padj)
    
    names(degs2017) <- c("geneID", "log2FC_2017", "p_2017", "padj_2017")
    
    
    print(paste0("Number of 2017 DEGs: ", nrow(degs2017)))
    
    degs2018 <- 
      read.delim(file2018) %>%
      dplyr::filter(padj < sig_level,
             log2FoldChange > fc_upper | log2FoldChange < fc_lower) %>%
      dplyr::select(geneID,
             log2FoldChange,
             pvalue,
             padj)
    
    names(degs2018) <- c("geneID", "log2FC_2018", "p_2018", "padj_2018")
    print(paste0("Number of 2018 DEGs: ", nrow(degs2018)))
    
    shared_degs <- 
      degs2017$geneID[degs2017$geneID %in% degs2018$geneID]
    
    shared_df <- 
      left_join(
        dplyr::filter(degs2017, geneID %in% shared_degs),
        dplyr::filter(degs2018, geneID %in% shared_degs),
        by = "geneID"
      ) %>%
      dplyr::select(
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
          avg_padj= (padj_2017 + padj_2018) /2,
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
blood_hai <-
  compare_years_blood("./output/processed_data/degs/Blood2017_HIGH-HR_4x_max_gmfr_DESEq2.tsv",
                      "./output/processed_data/degs/Blood2018_HIGH-HR_4x_max_gmfr_DESEq2.tsv",
                      0.1,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/blood_hai_shared.csv")

blood_nohai <-
  compare_years_blood("./output/processed_data/degs/Blood2017_LOW-HR_4x_max_gmfr_DESEq2.tsv",
                      "./output/processed_data/degs/Blood2018_LOW-HR_4x_max_gmfr_DESEq2.tsv",
                      0.1,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/blood_nohai_shared.csv")

ordered_hai_responder <- 
  blood_hai %>%
  mutate(avg_fc = (log2FC_2017 + log2FC_2018)/2) %>%
  filter(!(geneID %in% blood_nohai$geneID)) %>%
  arrange(desc(avg_fc))

ordered_hai_nonresponder <- 
  blood_nohai %>%
  mutate(avg_fc = (log2FC_2017 + log2FC_2018)/2) %>%
  filter(!(geneID %in% blood_hai$geneID)) %>%
  arrange(desc(avg_fc))


##Blood IgA
blood_iga <-
  compare_years_blood("./output/processed_data/degs/Blood2017_HIGH-HR_2x_max_iga_fc_DESEq2.tsv",
                      "./output/processed_data/degs/Blood2018_HIGH-HR_2x_max_iga_fc_DESEq2.tsv",
                      0.1,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/blood_iga_shared.csv")

blood_noiga <-
  compare_years_blood("./output/processed_data/degs/Blood2017_LOW-HR_2x_max_iga_fc_DESEq2.tsv",
                      "./output/processed_data/degs/Blood2018_LOW-HR_2x_max_iga_fc_DESEq2.tsv",
                      0.1,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/blood_noiga_shared.csv")

##CD4
blood_cd4 <-
  compare_years_blood("./output/processed_data/degs/Blood2017_HIGH-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                      "./output/processed_data/degs/Blood2018_HIGH-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                      0.1,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/blood_cd4_shared.csv")

blood_nocd4 <-
  compare_years_blood("./output/processed_data/degs/Blood2017_LOW-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                      "./output/processed_data/degs/Blood2018_LOW-HR_2x_max_mnp_cd4_fc_DESEq2.tsv",
                      0.1,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/blood_nocd4_shared.csv")

## CD8
blood_cd8 <-
  compare_years_blood("./output/processed_data/degs/Blood2017_HIGH-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                      "./output/processed_data/degs/Blood2018_HIGH-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                      0.1,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/blood_cd8_shared.csv")

blood_nocd8 <-
  compare_years_blood("./output/processed_data/degs/Blood2017_LOW-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                      "./output/processed_data/degs/Blood2018_LOW-HR_2x_max_mnp_cd8_fc_DESEq2.tsv",
                      0.1,
                      -0.322,
                      0.322,
                      "./output/processed_data/degs/unique/blood_nocd8_shared.csv")

####
## Unique DEGS
####

##Blood
blood_hai_degs <- get_unique_degs(blood_nohai, blood_hai)
saveRDS(blood_hai_degs, file = "./output/processed_data/degs/unique/blood_hai_degs.rds")
blood_iga_degs <- get_unique_degs(blood_noiga, blood_iga)
saveRDS(blood_iga_degs, file = "./output/processed_data/degs/unique/blood_iga_degs.rds")
blood_cd4_degs <- get_unique_degs(blood_nocd4, blood_cd4)
saveRDS(blood_cd4_degs, file = "./output/processed_data/degs/unique/blood_cd4_degs.rds")
blood_cd8_degs <- get_unique_degs(blood_nocd8, blood_cd8)
saveRDS(blood_cd8_degs, file = "./output/processed_data/degs/unique/blood_cd8_degs.rds")
