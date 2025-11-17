library(reshape)
library(plyr)
library(DESeq2)
library(dplyr)
library(corrplot)
library(ggplot2)
library(vsn)
library(BiocParallel)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(limma)
library(fastDummies)
options(stringsAsFactors = F)
register(MulticoreParam(4))

####
## Read in data
####

## RNAseq data
genetokeep <- read.delim("data/2017_nasal_rma_normalised.txt")
genetokeep <-genetokeep$geneid

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
  mastercoldata %>% dplyr::filter(id %in% subjects)


####
## DESeq2
####

for (cohort in 1:length(rawdata)) {
  for (immune in 5:ncol(mastercoldata)) {
    variable <- colnames(mastercoldata)[immune]
    filenameLow <- paste0("output/processed_data/degs/", names(rawdata)[cohort],"_LOW-",variable,"_DESEq2.tsv" )
    filenameHigh <- paste0("output/processed_data/degs/", names(rawdata)[cohort],"_HIGH-",variable,"_DESEq2.tsv" )
    data <- rawdata[[cohort]]
    keep <- rowSums(data) >= 30
    data <- data[keep,]
    data_ids <- gsub("_V.", "", colnames(data))
    coldata <- mastercoldata
    coldata <- coldata %>%
      dplyr::select(id, sample_id, timepoint, all_of(variable)) %>%
      dplyr::filter(sample_id %in% colnames(data))
    coldata <- coldata[complete.cases(coldata),]
    data <- data %>% dplyr::select(all_of(coldata$sample_id))
    row.names(coldata) <- coldata$sample_id
    coldata$sample_id <- NULL
    coldata$id <- factor(coldata$id)
    coldata$timepoint <- factor(coldata$timepoint, levels=c("V0", "V2"))
    colnames(coldata)[3] <- "group"

    coldataLow <- coldata %>% filter(group == 0) %>% dplyr::select(id, timepoint)
    coldataLow$id <- factor(coldataLow$id)
    dataLow  <- data %>% dplyr::select(all_of(row.names(coldataLow)))
    
    coldataHigh <- coldata %>% filter(group == 1) %>% dplyr::select(id, timepoint)
    coldataHigh$id <- factor(coldataHigh$id)
    dataHigh  <- data %>% dplyr::select(all_of(row.names(coldataHigh)))
    
    if(!(all(rownames(coldata) %in% colnames(data)))) {break}
    if(!(all(rownames(coldata) == colnames(data)))) {break}
    if(!(all(rownames(coldataLow) %in% colnames(dataLow)))) {break}
    if(!(all(rownames(coldataLow) == colnames(dataLow)))) {break}  
    if(!(all(rownames(coldataHigh) %in% colnames(dataHigh)))) {break}
    if(!(all(rownames(coldataHigh) == colnames(dataHigh)))) {break}
    if(length(rownames(coldataHigh)) < 6) {next}
        
    ddsLow <- DESeqDataSetFromMatrix(countData = dataLow, colData = coldataLow, design = ~ id + timepoint)
    ddsLow <- DESeq(ddsLow, parallel = TRUE)
    resLow <- results(ddsLow, parallel = TRUE)
    resOrderedLow<- resLow[order(resLow$log2FoldChange),]
    DEGsLow <- as.data.frame(resOrderedLow)
    DEGsLow$geneID <- rownames(DEGsLow)
    DEGsLow <- dplyr::select(DEGsLow, c("geneID", "baseMean", "log2FoldChange","lfcSE", "stat","pvalue","padj"))
    write.table(DEGsLow, filenameLow, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    ddsHigh <- DESeqDataSetFromMatrix(countData = dataHigh, colData = coldataHigh, design = ~ id + timepoint)
    ddsHigh <- DESeq(ddsHigh, parallel = TRUE)
    resHigh <- results(ddsHigh, parallel = TRUE)
    resOrderedHigh<- resHigh[order(resHigh$log2FoldChange),]
    DEGsHigh <- as.data.frame(resOrderedHigh)
    DEGsHigh$geneID <- rownames(DEGsHigh)
    DEGsHigh <- dplyr::select(DEGsHigh, c("geneID", "baseMean", "log2FoldChange","lfcSE", "stat","pvalue","padj"))
    write.table(DEGsHigh, filenameHigh, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    }
 }

########################################## BEGIN LIMMA BLOCK ##########################################

rawdata <- list() 
rawdata$Nasal2017 <- read.delim("data/2017_nasal_rma_normalised.txt")
row.names(rawdata$Nasal2017) <- rawdata$Nasal2017$geneid
rawdata$Nasal2017$geneid <- NULL

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
  mastercoldata %>% dplyr::filter(id %in% subjects)

for (cohort in 1:length(rawdata)) {
  for (immune in 5:ncol(mastercoldata)) {
    variable <- colnames(mastercoldata)[immune]
    filenameLow <- paste0("output/processed_data/degs/", names(rawdata)[cohort],"_LOW-",variable,"_Limma.tsv" )
    filenameHigh <- paste0("output/processed_data/degs/", names(rawdata)[cohort],"_HIGH-",variable,"_Limma.tsv" )
    data <- rawdata[[cohort]]
    mean <- as.data.frame(rowMeans(data))
    colnames(mean) <- "mean"
    mean$gene <- row.names(mean)
    mean <- mean %>% arrange(dplyr::desc(mean))
    keep <- mean[c(1:(nrow(mean)*0.8)),2]
    data <- data[keep,]
    data_ids <- gsub("_V.", "", colnames(data))
    coldata <- mastercoldata
    coldata <- coldata %>%
      dplyr::select(id, sample_id, timepoint, all_of(variable)) %>%
      dplyr::filter(sample_id %in% colnames(data))
    coldata <- coldata[complete.cases(coldata),]
    data <- data %>% dplyr::select(all_of(coldata$sample_id))
    row.names(coldata) <- coldata$sample_id
    coldata$sample_id <- NULL
    coldata$id <- factor(coldata$id)
    coldata$timepoint <- factor(coldata$timepoint, levels=c("V0", "V2"))
    colnames(coldata)[3] <- "group"
    
    coldataLow <- coldata %>% dplyr::filter(group == 0) %>% dplyr::select(id, timepoint)
    coldataLow$id <- factor(coldataLow$id)
    dataLow  <- data %>% dplyr::select(all_of(row.names(coldataLow)))
    coldataHigh <- coldata %>% dplyr::filter(group == 1) %>% dplyr::select(id, timepoint)
    coldataHigh$id <- factor(coldataHigh$id)
    dataHigh  <- data %>% dplyr::select(all_of(row.names(coldataHigh)))
    
    if(!(all(rownames(coldata) %in% colnames(data)))) {break}
    if(!(all(rownames(coldata) == colnames(data)))) {break}
    if(!(all(rownames(coldataLow) %in% colnames(dataLow)))) {break}
    if(!(all(rownames(coldataLow) == colnames(dataLow)))) {break}  
    if(!(all(rownames(coldataHigh) %in% colnames(dataHigh)))) {break}
    if(!(all(rownames(coldataHigh) == colnames(dataHigh)))) {break}
    if(length(rownames(coldataHigh)) < 6) {next}
    if(length(rownames(coldataLow)) < 6) {next}
    
    # Define model 
    design <- model.matrix(~0 + coldataLow$timepoint + coldataLow$id)
    colnames(design) <- gsub('coldataLow\\$id|coldataLow\\$timepoint', '', colnames(design))
    
    fit <- lmFit(dataLow, design)
    contrast.matrix <- makeContrasts(contrasts = "V2-V0", levels = design)
    
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)

    coef_name <- colnames(fit2$coefficients)[1]  
    # add in moderated SE for downstream meta-analysis: sqrt(s2.post) * stdev.unscaled
    se_mod <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
    se_mod_vec <- se_mod[, coef_name]
    
    res_paired <- topTable(fit2, coef = coef_name, n = Inf, adjust = "BH")
    res_paired$SE_mod <- se_mod_vec[rownames(res_paired)]
    res_paired$geneID <- rownames(res_paired)
    
    write.table(res_paired, filenameLow, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    rm(fit,fit2,res_paired)
    
    design <- model.matrix(~0 + coldataHigh$timepoint + coldataHigh$id)
    colnames(design) <- gsub('coldataHigh\\$id|coldataHigh\\$timepoint', '', colnames(design))
    fit <- lmFit(dataHigh, design)
    contrast.matrix <- makeContrasts(contrasts = "V2-V0", levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    coef_name <- colnames(fit2$coefficients)[1]
    se_mod <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
    se_mod_vec <- se_mod[, coef_name]
    
    res_paired <- topTable(fit2, coef = coef_name, n = Inf, adjust = "BH")
    res_paired$SE_mod <- se_mod_vec[rownames(res_paired)]
    res_paired$geneID <- rownames(res_paired)

    write.table(res_paired, filenameHigh, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
    
  }}
