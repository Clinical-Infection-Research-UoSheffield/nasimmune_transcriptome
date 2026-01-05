library(dplyr)

metadata <- readRDS("./output/processed_data/metadata_cleaned.rds")

# Get pids with transcriptome available
rawdata <- list() 
rawdata$Nasl17 <- read.delim("./data/2017_nasal_rma_normalised.txt")
row.names(rawdata$Nasl17) <- rawdata$Nasl17$geneid
rawdata$Nasl17$geneid <- NULL
rawdata$Nasl18 <- read.delim("./data/2018_nasal_counts.tsv")
rawdata$Bld17 <- read.delim("./data/2017_blood_counts.tsv")
rawdata$Bld18 <- read.delim("./data/2018_blood_counts.tsv")

bld_pids <- c(names(rawdata$Bld17), names(rawdata$Bld18))
bld_pids <- unique(substr(bld_pids, 1, 5))

nsl_pids <- c(names(rawdata$Nasl17), names(rawdata$Nasl18))
nsl_pids <- unique(substr(nsl_pids, 1, 5))

# label metadata
metadata$bld_transcriptome <- metadata$subid1 %in% bld_pids
metadata$nsl_transcriptome <- metadata$subid1 %in% nsl_pids
metadata$both_transcriptome <- metadata$bld_transcriptome * metadata$nsl_transcriptome

#sex
table(metadata$sex[metadata$bld_transcriptome==1])
prop.table(table(metadata$sex[metadata$bld_transcriptome==1]))
table(metadata$sex[metadata$bld_transcriptome==0])
prop.table(table(metadata$sex[metadata$bld_transcriptome==0]))
chisq.test(metadata$sex, metadata$bld_transcriptome)

table(metadata$sex[metadata$nsl_transcriptome==1])
prop.table(table(metadata$sex[metadata$nsl_transcriptome==1]))
table(metadata$sex[metadata$nsl_transcriptome==0])
prop.table(table(metadata$sex[metadata$nsl_transcriptome==0]))
chisq.test(metadata$sex, metadata$nsl_transcriptome)

#age 
summary(metadata$age[metadata$bld_transcriptome==1])
summary(metadata$age[metadata$bld_transcriptome==0])
wilcox.test(metadata$age[metadata$bld_transcriptome==1], metadata$age[metadata$bld_transcriptome==0])

summary(metadata$age[metadata$nsl_transcriptome==1])
summary(metadata$age[metadata$nsl_transcriptome==0])
wilcox.test(metadata$age[metadata$nsl_transcriptome==1], metadata$age[metadata$nsl_transcriptome==0])

#seropos
table(metadata$h1_v0_gmt[metadata$bld_transcriptome==1]==5)
prop.table(table(metadata$h1_v0_gmt[metadata$bld_transcriptome==1]==5))
table(metadata$h1_v0_gmt[metadata$bld_transcriptome==0]==5)
prop.table(table(metadata$h1_v0_gmt[metadata$bld_transcriptome==0]==5))
chisq.test(metadata$h1_v0_gmt==5, metadata$bld_transcriptome)

table(metadata$h3_v0_gmt[metadata$bld_transcriptome==1]==5)
prop.table(table(metadata$h3_v0_gmt[metadata$bld_transcriptome==1]==5))
table(metadata$h3_v0_gmt[metadata$bld_transcriptome==0]==5)
prop.table(table(metadata$h3_v0_gmt[metadata$bld_transcriptome==0]==5))
chisq.test(metadata$h3_v0_gmt==5, metadata$bld_transcriptome)

table(metadata$b_v0_gmt[metadata$bld_transcriptome==1]==5)
prop.table(table(metadata$b_v0_gmt[metadata$bld_transcriptome==1]==5))
table(metadata$b_v0_gmt[metadata$bld_transcriptome==0]==5)
prop.table(table(metadata$b_v0_gmt[metadata$bld_transcriptome==0]==5))
chisq.test(metadata$b_v0_gmt==5, metadata$bld_transcriptome)

table(metadata$h1_v0_gmt[metadata$nsl_transcriptome==1]==5)
prop.table(table(metadata$h1_v0_gmt[metadata$nsl_transcriptome==1]==5))
table(metadata$h1_v0_gmt[metadata$nsl_transcriptome==0]==5)
prop.table(table(metadata$h1_v0_gmt[metadata$nsl_transcriptome==0]==5))
chisq.test(metadata$h1_v0_gmt==5, metadata$nsl_transcriptome)

table(metadata$h3_v0_gmt[metadata$nsl_transcriptome==1]==5)
prop.table(table(metadata$h3_v0_gmt[metadata$nsl_transcriptome==1]==5))
table(metadata$h3_v0_gmt[metadata$nsl_transcriptome==0]==5)
prop.table(table(metadata$h3_v0_gmt[metadata$nsl_transcriptome==0]==5))
chisq.test(metadata$h3_v0_gmt==5, metadata$nsl_transcriptome)

table(metadata$b_v0_gmt[metadata$nsl_transcriptome==1]==5)
prop.table(table(metadata$b_v0_gmt[metadata$nsl_transcriptome==1]==5))
table(metadata$b_v0_gmt[metadata$nsl_transcriptome==0]==5)
prop.table(table(metadata$b_v0_gmt[metadata$nsl_transcriptome==0]==5))
chisq.test(metadata$b_v0_gmt==5, metadata$nsl_transcriptome)
