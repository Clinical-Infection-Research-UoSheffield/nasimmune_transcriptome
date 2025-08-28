library(dplyr)

# read in metadata
metadata <- read.csv("./data/metadata.csv")

# create binary variables for each immune response

metadata$max_iga_fc <- apply(metadata[,c("h1_iga_fc", "h3_iga_fc", "b_iga_fc")],1,max)
metadata$HR_2x_max_iga_fc <- ifelse(log10(metadata$max_iga_fc) >= log10(2), 1, 0)

metadata$max_gmfr <- apply(metadata[,c("h1_gmfr", "h3_gmfr", "b_gmfr")],1,max)
metadata$HR_4x_max_gmfr <- ifelse(log2(metadata$max_gmfr) >= log2(4), 1, 0)

metadata$max_mnp_cd4_fc <- apply(metadata[,c("mnp_cd4_ifng_fc", "mnp_cd4_il2_fc")],1,max)
metadata$HR_2x_max_mnp_cd4_fc <- ifelse(log10(metadata$max_mnp_cd4_fc) >= log10(2), 1, 0)

metadata$max_mnp_cd8_fc <- apply(metadata[,c("mnp_cd8_ifng_fc", "mnp_cd8_il2_fc")],1,max)
metadata$HR_2x_max_mnp_cd8_fc <- ifelse(log10(metadata$max_mnp_cd8_fc) >= log10(2), 1, 0)

metadata$SheddingD2Sum <- rowSums(metadata[,c("h1_v2_shed", "h3_v2_shed", "b_v2_shed")])
metadata$SheddingD7Sum <- rowSums(metadata[,c("h1_v7_shed", "h3_v7_shed", "b_v7_shed")])
metadata$N_responses <-  rowSums(metadata[,c("HR_2x_max_iga_fc","HR_4x_max_gmfr","HR_2x_max_mnp_cd4_fc", "HR_2x_max_mnp_cd8_fc")])

# output cleaned metadata
saveRDS(metadata, "./output/processed_data/metadata_cleaned.rds")
