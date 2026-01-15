library(tidyr)
library(dplyr)

####
## Get unique HAI DEGs
####
blood_hai <- readRDS("./output/processed_data/degs/unique/blood_hai_degs.rds")
genes <- blood_hai$unique_responders$geneID[blood_hai$unique_responders$up_down == "Up"]


####
## Read in expression data
####

## RNAseq data
genetokeep <- read.delim("data/2017_nasal_rma_normalised.txt")
genetokeep <- genetokeep$geneid

rawdata <- list() 
rawdata$Nasal2018 <- read.delim("data/2018_nasal_counts.tsv", sep = "\t")
sharedgenes <- intersect(row.names(rawdata$Nasal2018),genetokeep)
rawdata$Nasal2018 <- rawdata$Nasal2018[sharedgenes,]
rawdata$Blood2017 <- read.delim("data/2017_blood_counts.tsv", sep = "\t")
sharedgenes <- intersect(row.names(rawdata$Blood2017),genetokeep )
rawdata$Blood2017 <- rawdata$Blood2017[sharedgenes,]
rawdata$Blood2018 <- read.delim("data/2018_blood_counts.tsv", sep = "\t")
sharedgenes <- intersect(row.names(rawdata$Blood2018),genetokeep )
rawdata$Blood2018 <- rawdata$Blood2018[sharedgenes,]

## Metadata
samples <- read.delim("output/processed_data/mastercoldata.tsv")
samples <- samples$samples_id

responsetable <- readRDS("output/processed_data/metadata_cleaned.rds")
responsetable$sheddD2_2orMore<- ifelse(responsetable$SheddingD2Sum >= 2, 1, 0)
responsetable$sheddD2_1orMore<- ifelse(responsetable$SheddingD2Sum >= 1, 1, 0)
responsetable$sheddD7_2orMore<- ifelse(responsetable$SheddingD7Sum >= 2, 1, 0)
responsetable$sheddD7_1orMore<- ifelse(responsetable$SheddingD7Sum >= 1, 1, 0)
responsetable$sheddD2_Noshedd<- ifelse(responsetable$SheddingD2Sum == 0, 1, 0)
responsetable$sheddD7_Noshedd<- ifelse(responsetable$SheddingD7Sum == 0, 1, 0)

responsetable <- 
  responsetable %>% 
  dplyr::select(subid1, year, HR_2x_max_iga_fc,
         HR_4x_max_gmfr, HR_2x_max_mnp_cd4_fc, HR_2x_max_mnp_cd8_fc,
         sheddD2_Noshedd, sheddD2_1orMore,
         sheddD2_2orMore,sheddD7_Noshedd,
         sheddD7_1orMore, sheddD7_2orMore)

responsetableV0 <- responsetable
responsetableV0$sample_id <- paste0(responsetableV0$subid1, "_V0")
responsetableV2 <- responsetable
responsetableV2$sample_id <- paste0(responsetableV0$subid1, "_V2")


mastercoldata <- rbind(responsetableV0, responsetableV2) %>% arrange(subid1)
mastercoldata <- mastercoldata %>% dplyr::filter(sample_id  %in% samples)
mastercoldata$timepoint <- gsub("^[A-Z][0-9][0-9][0-9][A-Z]_", "", mastercoldata$sample_id)
colnames(mastercoldata)[1] <- "id"
mastercoldata <- 
  mastercoldata %>% 
  dplyr::select(id, year, sample_id, timepoint,  
         HR_2x_max_iga_fc, HR_4x_max_gmfr, HR_2x_max_mnp_cd4_fc,
         HR_2x_max_mnp_cd8_fc,
         sheddD2_Noshedd, sheddD2_1orMore, sheddD2_2orMore,
         sheddD7_Noshedd, 
         sheddD7_1orMore, sheddD7_2orMore
  )

subjects <- 
  lapply(1:length(rawdata),
         function(cohort){
           unique(substr(colnames(rawdata[[cohort]]), 1, 5))
         })

subjects <-  unlist(subjects)

mastercoldata <- 
  mastercoldata %>% filter(id %in% subjects)

# New dataset
blood <- bind_cols(rawdata$Blood2017, rawdata$Blood2018)

df_long <- blood %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "counts") %>%
  mutate(
    subid = substr(sample, 1, 4),
    visit = substr(sample, nchar(sample)-1, nchar(sample))
  ) %>%
  dplyr::select(subid, gene, visit, counts)

# Make V0 and V2 columns
df_wide_visits <- df_long %>%
  pivot_wider(names_from = visit, values_from = counts) %>%
  arrange(subid) %>%
  dplyr:: filter(gene %in% genes) 

names(df_wide_visits) <- c("subid", "gene", "V0_count", "V2_count")

#same with fold change
fc2017 <- read.table("./output/processed_data/Blood_2017FiltredToRMAgenes_rlog_log2FC_Ranks.tsv", sep = "\t")
fc2017 <- fc2017[genes,]

fc2018 <- read.table("./output/processed_data/Blood_2018FiltredToRMAgenes_rlog_log2FC_Ranks.tsv", sep = "\t")
fc2018 <- fc2018[genes,]

fc <- bind_cols(fc2017, fc2018)
df_t <- as.data.frame(t(fc))
colnames(df_t) <- rownames(fc)
df_t$subid <- rownames(df_t)
df_t <- df_t[, c("subid", setdiff(names(df_t), "subid"))]

#metadata to merge
metadata <- responsetable %>% dplyr::select("subid1", "HR_4x_max_gmfr")
colnames(metadata) <- c("subid", "seroconversion")
df_t <- left_join(df_t, metadata)
df_t <- dplyr::select(df_t, -c(subid))
write.csv(df_t, "./output/processed_data/hai_degs_for_gene_score.csv", row.names = F, quote = F)
