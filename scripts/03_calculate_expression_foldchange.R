library(DESeq2)
library(dplyr)
library(edgeR)
library(BiocParallel)
register(MulticoreParam(3))

####
## Takes raw microarray and RNAseq data as input, log transforms and generates
## individual gene expression fold change between baseline and day 2 post vaccine
####

####
## Read data
####
rawdata <- list() 
rawdata$Nasl17 <- read.delim("./data/2017_nasal_rma_normalised.txt")
row.names(rawdata$Nasl17) <- rawdata$Nasl17$geneid
rawdata$Nasl17$geneid <- NULL
rawdata$Nasl18 <- read.delim("./data/2018_nasal_counts.tsv")
rawdata$Bld17 <- read.delim("./data/2017_blood_counts.tsv")
rawdata$Bld18 <- read.delim("./data/2018_blood_counts.tsv")

####
## Filter data
####
# Select genes according to Microarray data because it has fewer genes and it could bias 
# the GSEA and some downstream analysis while comparing datasets 
marray_genes <- row.names(rawdata$Nasl17)
blood_genes <- intersect(row.names(rawdata$Bld17),row.names(rawdata$Bld18))
blood_genes <- intersect(blood_genes,marray_genes)
nasal_genes <- intersect(marray_genes, row.names(rawdata$Nasl18))

genelist <- list()
genelist$Nasl17 <- marray_genes
genelist$Nasl18 <- marray_genes
genelist$Bld17 <- marray_genes
genelist$Bld18 <- marray_genes

####
## Log2 scale and normalise with respect to library size
####

# Generate baseline files for further analysis (if necessary)
baseline <- list()
baseline$Nasl17 <- rawdata$Nasl17 %>% select(all_of(grep("_V0", colnames(rawdata$Nasl17),value = T)))
baseline$Nasl17 <- baseline$Nasl17[nasal_genes,]
baseline$Nasl18 <- rawdata$Nasl18 %>% select(all_of(grep("_V0", colnames(rawdata$Nasl18),value = T)))
baseline$Nasl18 <- baseline$Nasl18[nasal_genes,]
baseline$Bld17 <- rawdata$Bld17 %>% select(all_of(grep("_V0", colnames(rawdata$Bld17),value = T))) 
baseline$Bld17 <- baseline$Bld17[blood_genes,]
baseline$Bld18 <- rawdata$Bld18 %>% select(all_of(grep("_V0", colnames(rawdata$Bld18),value = T))) 
baseline$Bld18 <- baseline$Bld18[blood_genes,]

for (df in 2:length(baseline))  { 
  data <- as.matrix(baseline[[df]])
  rlogdata <- rlog(data, blind=TRUE)
  rlogdata[which(is.infinite(rlogdata))] <- 0
  rlogdata <- as.data.frame(rlogdata)
  baseline[[df]] <- rlogdata
}

# Write baseline files
write.table(baseline$Bld17,"./output/processed_data/blood_2017_rlog_baseline.tsv", quote = F, sep = "\t")
write.table(baseline$Bld18,"./output/processed_data/blood_2018_rlog_baseline.tsv", quote = F, sep = "\t")
write.table(baseline$Nasl17,"./output/processed_data/nasal_2017_rma_baseline.tsv", quote = F, sep = "\t")
write.table(baseline$Nasl18,"./output/processed_data/nasal_2018_rlog_baseline.tsv", quote = F, sep = "\t")

####
## Generate individual participant foldchange data
####

CalculateFC <- function(eset, sample_name) {
  v0 <- eset[, paste0(sample_name, "_V0")]
  v2 <- eset[, paste0(sample_name, "_V2")]
  v2-v0
}

data_foldChange <- list()
data_foldChangeFiltred <- list()

# Processing the MicroArray data 
datanames <- list()
datanames$Nasl17 <- "Nasal_2017"
datanames$Nasl18 <- "Nasal_2018"
datanames$Bld17 <- "Blood_2017"
datanames$Bld18 <- "Blood_2018"

data <- rawdata$Nasl17
mean <- as.data.frame(rowMeans(data))
colnames(mean) <- "mean"
mean$gene <- row.names(mean)
mean <- mean %>% arrange(dplyr::desc(mean))
keep <- mean[c(1:(nrow(mean)*0.8)),2]
data <- data[keep,]
rawdata$Nasl17 <- data

micro_ids <- unique(gsub("_V.", "", colnames(rawdata$Nasl17)))
micro_fc <- data.frame(geneid = rownames(rawdata$Nasl17))

for(id in micro_ids) {
  fc <- CalculateFC(rawdata$Nasl17, id)
  micro_fc[, id] <- fc
}
micro_fc[, -1] <- do.call(data.frame, lapply(micro_fc[, -1], function(x) {
  replace(x, is.infinite(x) | is.na(x), 0)
})
)

filename <- paste0("./output/processed_data/", datanames[[1]], "_FiltredToSharedGenes_log2FC_Ranks.tsv" )
row.names(micro_fc) <-  micro_fc$geneid
micro_fc$geneid <- NULL
write.table(micro_fc, filename, sep = "\t", quote = F)
data_foldChange$Nasl17 <- micro_fc
data_foldChangeFiltred$Nasl17 <- micro_fc

Samplelist <- colnames(rawdata$Nasl18) 
Samplelist <- c(Samplelist, colnames(rawdata$Nasl17))
Samplelist <- c(Samplelist, colnames(rawdata$Bld18))
Samplelist <- c(Samplelist, colnames(rawdata$Bld17))
Samplelist <- unique(Samplelist) 

ids <- gsub("_V.", "", Samplelist)
timepoint <- gsub("^[A-Z][0-9][0-9][0-9][A-Z]_", "", Samplelist)
mastercoldata <- data.frame(samples_id=Samplelist, id=ids,timepoint=timepoint ) 

write.table(mastercoldata, "./output/processed_data/mastercoldata.tsv", quote = F, row.names = F, sep = "\t")

#Processing the RNAseq data
for (df in 2:length(rawdata))  { #Skipping the micro array data frame
  filename <- paste0("./output/processed_data/", datanames[[df]], "_TMM_log2FC.tsv" )
  filenameFCrlog <- paste0("./output/processed_data/", datanames[[df]], "_rlog_log2FC.tsv" )
  filenamefiltred <- paste0("./output/processed_data/", datanames[[df]], "FiltredToRMAgenes_TMM_log2FC_Ranks.tsv" )
  filenamefiltredrlog <- paste0("./output/processed_data/", datanames[[df]], "FiltredToRMAgenes_rlog_log2FC_Ranks.tsv" )
  filenamerlog <- paste0("./output/processed_data/", datanames[[df]], "_Rlog_normalized.tsv" )
  filenameTMM <- paste0("./output/processed_data/", datanames[[df]], "_TMM_normalized.tsv" )
  filenamefiltredvst <- paste0("./output/processed_data/", datanames[[df]], "FiltredToRMAgenes_vst_log2FC_Ranks.tsv" )
  filenamervst <- paste0("./output/processed_data/", datanames[[df]], "_vst_normalized.tsv" )
  
  dataraw <- rawdata[[df]]
  cohort <- datanames[[df]]
  # Get IDs
  data_ids <- unique(gsub("_V.", "", colnames(data)))
  
  # Normalizing raw counts
  keep <- rowSums(dataraw) >= 30
  dataraw <- dataraw[keep,]
  data_ids <- gsub("_V.", "", colnames(dataraw))
  data_groups <- gsub("....._", "", colnames(dataraw))
  exp <- DGEList(counts=dataraw,group=factor(data_groups))
  exp <- calcNormFactors(exp, method = "TMM")
  data <-  cpm(exp, log=TRUE)
  
  
  coldata <- mastercoldata
  coldata <- coldata %>%
    select(id, samples_id, timepoint) %>%
    filter(samples_id %in% colnames(dataraw))
  row.names(coldata) <- coldata$samples_id
  coldata$samples_id <- NULL
  coldata$id <- factor(coldata$id)
  coldata$timepoint <- factor(coldata$timepoint, levels=c("V0", "V2"))
  
  coldata <-  coldata[colnames(dataraw),]
  
  if(!(all(rownames(coldata) %in% colnames(dataraw)))) {break}
  if(!(all(rownames(coldata) == colnames(dataraw)))) {break}
  
  dds <- DESeqDataSetFromMatrix(countData = dataraw, colData = coldata, design = ~ id + timepoint)
  
  rlogdata <- rlog(dds, blind=FALSE)
  vstgdata <- vst(dds, blind=FALSE)
  write.table(assay(vstgdata), filenamervst, sep = "\t", quote = F)
  write.table(assay(rlogdata), filenamerlog, sep = "\t", quote = F)
  write.table(data, filenameTMM, sep = "\t", quote = F)
  rlogtable <- as.data.frame(assay(rlogdata))
  
  # Calculate FC
  data_fc <- data.frame(geneid = rownames(data))
  for(id in data_ids) {
    fc <- CalculateFC(data, id)
    data_fc[, id] <- fc
  }
  data_fcfiltred <- data_fc %>% filter(geneid %in% genelist[[df]])
  row.names(data_fc) <- data_fc$geneid
  data_fc$geneid <- NULL
  row.names(data_fcfiltred) <- data_fcfiltred$geneid
  data_fcfiltred$geneid <- NULL
  data_foldChangeFiltred[[df]] <- data_fcfiltred
  data_foldChange[[df]] <- as.data.frame(data_fc)
  write.table(data_fcfiltred, filenamefiltred, sep = "\t", quote = F)
  write.table(data_fc, filename, sep = "\t", quote = F)
  
  
  data_fc <- data.frame(geneid = rownames(rlogtable)) 
  for(id in data_ids) {
    fc <- CalculateFC(rlogtable, id)
    data_fc[, id] <- fc
  }
  data_fcfiltred <- data_fc %>% filter(geneid %in% genelist[[df]])
  row.names(data_fc) <- data_fc$geneid
  data_fc$geneid <- NULL
  row.names(data_fcfiltred) <- data_fcfiltred$geneid
  
  data_fcfiltred$geneid <- NULL
  data_foldChangeFiltred[[df]] <- data_fcfiltred
  data_foldChange[[df]] <- as.data.frame(data_fc)
  write.table(data_fcfiltred, filenamefiltredrlog, sep = "\t", quote = F)
  write.table(data_fc, filenameFCrlog, sep = "\t", quote = F)
  
}
