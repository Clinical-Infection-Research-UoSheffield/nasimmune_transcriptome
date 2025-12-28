library(dplyr)

metadata <- readRDS("./output/processed_data/metadata_cleaned.rds")
#metadata <- metadata[!is.na(metadata$N_responses),]

# Sex
sex_2017 <- table(metadata$sex[metadata$year==2017])
sex_2017
prop.table(sex_2017)

sex_2018 <- table(metadata$sex[metadata$year==2018])
sex_2018
prop.table(sex_2018)

chisq.test(rbind(sex_2017, sex_2018))

sex_total <- table(metadata$sex)
sex_total
prop.table(sex_total)

# age
summary(metadata$age[metadata$year==2017])

summary(metadata$age[metadata$year==2018])

wilcox.test(metadata$age[metadata$year==2017], metadata$age[metadata$year==2018])

summary(metadata$age)

#seropositive
h1_2017_seropos <- table(metadata$h1_v0_gmt[metadata$year==2017]!=5)
h1_2017_seropos
prop.table(h1_2017_seropos)

h1_2018_seropos <- table(metadata$h1_v0_gmt[metadata$year==2018]!=5)
h1_2018_seropos
prop.table(h1_2018_seropos)

chisq.test(rbind(h1_2017_seropos, h1_2018_seropos))

h3_2017_seropos <- table(metadata$h3_v0_gmt[metadata$year==2017]!=5)
h3_2017_seropos
prop.table(h3_2017_seropos)

h3_2018_seropos <- table(metadata$h3_v0_gmt[metadata$year==2018]!=5)
h3_2018_seropos
prop.table(h3_2018_seropos)

chisq.test(rbind(h3_2017_seropos, h3_2018_seropos))

b_2017_seropos <- table(metadata$b_v0_gmt[metadata$year==2017]!=5)
b_2017_seropos
prop.table(b_2017_seropos)

b_2018_seropos <- table(metadata$b_v0_gmt[metadata$year==2018]!=5)
b_2018_seropos
prop.table(b_2018_seropos)

chisq.test(rbind(b_2017_seropos, b_2018_seropos))


h1_total_seropos <- table(metadata$h1_v0_gmt!=5)
h1_total_seropos
prop.table(h1_total_seropos)

h3_total_seropos <- table(metadata$h3_v0_gmt!=5)
h3_total_seropos
prop.table(h3_total_seropos)

b_total_seropos <- table(metadata$b_v0_gmt!=5)
b_total_seropos
prop.table(b_total_seropos)
