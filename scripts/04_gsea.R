library(data.table)
library(clusterProfiler)
library(fgsea)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

####
## Performs geneset enrichment on individual gene expression fold change
####
options(stringsAsFactors = F)
####
## Function block
####
read.gmt <- function(fname){
  res <- list(genes=list(), 
              desc=list())
  gmt <- file(fname)
  gmt.lines <- readLines(gmt)
  close(gmt)
  gmt.list <- lapply(gmt.lines, 
                     function(x) unlist(strsplit(x, split="\t")))
  gmt.names <- sapply(gmt.list, '[', 1)
  gmt.desc <- lapply(gmt.list, '[', 2)
  gmt.genes <- lapply(gmt.list,
                      function(x){x[3:length(x)]})
  names(gmt.desc) <- names(gmt.genes) <- gmt.names
  return(gmt.genes)
}
ssGSEA <- function(gmtfile=gmtfile,fileranks=fileranks,Ptype=Ptype,pval_cutoff=pval_cutoff){
  
  gmt <- read.gmt(gmtfile)
  ranks_df  <- read.delim(fileranks)
  names <- colnames(ranks_df)
  ranks_df$geneID <- row.names(ranks_df)
  ranks_df <- ranks_df %>% select(geneID, all_of(names)) 
  ranks_df <- ranks_df[complete.cases(ranks_df), ]
  ranks_ch <- colnames(ranks_df)[-1]
  ranks_ch <- setNames(ranks_ch, ranks_ch)
  
  #run fastGSEA
  tmp_ranks <- lapply(ranks_ch, function(rankname){
    tmpdf <- ranks_df[,c('geneID', rankname)]
    tmpdf <- tmpdf[complete.cases(tmpdf),]
    tmpranks <- tmpdf[[rankname]]
    names(tmpranks) <- tmpdf$geneID
    tmpranks <- tmpranks[order(tmpranks)]
    # print(tmpranks)
    fgseaRes <- fgsea(pathways = gmt, stats = tmpranks, 
                      minSize = 15, maxSize = 2000, nperm = 1000)
    fgseaRes <- as.data.frame(fgseaRes)
    fgseaRes <- fgseaRes[,c('pathway', 'pval', 'padj', 'NES','size','leadingEdge')]
    fgseaRes
  })
  
  
  # Remove ranks without enrichment
  tmp_ranks <- Filter(function(x) nrow(x) > 1, tmp_ranks)
  
  # Write output NES
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"NES"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
  }
  df <- df[,-1]
  
  fileranks <- str_remove(fileranks, "./output/processed_data/")
  fileranks <- str_remove(fileranks, "FiltredToRMAgenes_rlog_log2FC_Ranks")
  fileranks <- str_remove(fileranks, "_FiltredToSharedGenes_log2FC_Ranks")
  
  rank_nameout <- paste0('./output/processed_data/gsea/NES_', fileranks,sep="")
  
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  # Write output AdjP
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"padj"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
  }
  df <- df[,-1]
  
  rank_nameout <- paste0('./output/processed_data/gsea/padj_', fileranks,sep="")
  
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  # 
  
  # Write output pval
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"pval"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
  }
  df <- df[,-1]
  rank_nameout <- paste0('./output/processed_data/gsea/pval_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  # Write output Leading Edge genes
  rm(df)
  for(name_rank in names(tmp_ranks)){
    rank_out <- tmp_ranks[[name_rank]] 
    xxx   <- rank_out[,"leadingEdge"]
    names(xxx) <- rank_out[,"pathway"]
    df <- cbind(df,xxx)
    colnames(df)[ncol(df)] <- name_rank
  }
  df <- df[,-1]
  rank_nameout <- paste0('./output/processed_data/gsea/LE_', fileranks,sep="")
  write.table(df, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
  
  
  nes <- data.table::fread(paste0('./output/processed_data/gsea/NES_', fileranks,sep=""))
  colnames(nes)[1] <- "pathway"
  pval <- data.table::fread(paste0('./output/processed_data/gsea/',Ptype,'_', fileranks,sep=""))
  colnames(pval)[1] <- "pathway"
  
  nes_melt <- nes %>%
    gather(sample, nes, -pathway)
  
  pval_melt <- pval %>%
    gather(sample, pval, -pathway)
  
  result <- full_join(nes_melt, pval_melt, by=c("pathway", "sample")) %>%
    filter(pval <= pval_cutoff) %>%
    select(pathway, sample, nes) %>%
    spread(sample, nes)
  
  rank_nameout <- paste0('./output/processed_data/gsea/NES_',Ptype,pval_cutoff,"_", fileranks,sep="")
  write.table(result, file = rank_nameout, sep = '\t', quote = F, row.names = T,col.names=NA)
}

####
## Run GSEA on each dataset using expression fold change
####

###  Blood 2017
gmtfile <- "./data/GO_Biological_Process_2018.gmt"
fileranks <- "./output/processed_data/Blood_2017FiltredToRMAgenes_rlog_log2FC_Ranks.tsv"
Ptype <- "padj"
pval_cutoff <- 0.1

#Run ssGSEA
ssGSEA(gmtfile = gmtfile, 
       fileranks = fileranks,
       Ptype = Ptype,
       pval_cutoff = pval_cutoff)

###  Blood 2018
gmtfile <- "./data/GO_Biological_Process_2018.gmt"
fileranks <- "./output/processed_data/Blood_2018FiltredToRMAgenes_rlog_log2FC_Ranks.tsv"
Ptype <- "padj"
pval_cutoff <- 0.1

#Run ssGSEA
ssGSEA(gmtfile = gmtfile, 
       fileranks = fileranks,
       Ptype = Ptype,
       pval_cutoff = pval_cutoff)

###  Nasal 2017
gmtfile <- "./data/GO_Biological_Process_2018.gmt"
fileranks <- "./output/processed_data/Nasal_2017_FiltredToSharedGenes_log2FC_Ranks.tsv"
Ptype <- "padj"
pval_cutoff <- 0.1

ssGSEA(gmtfile = gmtfile, 
       fileranks = fileranks,
       Ptype = Ptype,
       pval_cutoff = pval_cutoff)


###  Nasal 2018
gmtfile <- "./data/GO_Biological_Process_2018.gmt"
fileranks <- "./output/processed_data/Nasal_2018FiltredToRMAgenes_rlog_log2FC_Ranks.tsv"
Ptype <- "padj"
pval_cutoff <- 0.1

ssGSEA(gmtfile = gmtfile, 
       fileranks = fileranks,
       Ptype = Ptype,
       pval_cutoff = pval_cutoff)
