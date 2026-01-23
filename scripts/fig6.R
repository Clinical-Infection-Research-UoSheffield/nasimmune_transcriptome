library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)

####
## Make Fig 6C. Fig 6A-C make in matplotlib in a seperate script
####

genes <- c("KIF4A", "DIAPH3", "IL17D", "ARHGAP23", "KIR2DS4", "SHROOM3", "GBP1", "SASH1", "MID1", "LTBP1")

####
## Read in fold change
####

#same with fold change
fc2017 <- read.table("./output/processed_data/Blood_2017FiltredToRMAgenes_rlog_log2FC_Ranks.tsv", sep = "\t")
fc2017 <- fc2017[genes,]

fc2018 <- read.table("./output/processed_data/Blood_2018FiltredToRMAgenes_rlog_log2FC_Ranks.tsv", sep = "\t")
fc2018 <- fc2018[genes,]

fc <- bind_cols(fc2017, fc2018)

#metadata to merge
## Metadata
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
metadata <- responsetable %>% dplyr::select("subid1", "HR_4x_max_gmfr")
metadata <- metadata %>% dplyr::filter(subid1 %in% colnames(fc))
rm(fc2017, fc2018, responsetable, responsetableV0, responsetableV2)

## 
metadata$fc_sum <- colSums(fc)

ggplot(metadata, aes(x = as.factor(HR_4x_max_gmfr), y = fc_sum)) +
    geom_boxplot() 
